import cv2 as cv
import numpy as np
import os
import pandas as pd
import time
import warnings

import array_analyzer.extract.background_estimator as background_estimator
import array_analyzer.extract.image_parser as image_parser
import array_analyzer.extract.img_processing as img_processing
import array_analyzer.extract.txt_parser as txt_parser
from array_analyzer.extract.metadata import MetaData
import array_analyzer.extract.constants as constants
import array_analyzer.load.debug_plots as debug_plots
import array_analyzer.transform.point_registration as registration
import array_analyzer.transform.array_generation as array_gen
import array_analyzer.utils.io_utils as io_utils
import array_analyzer.load.report as report


def point_registration(input_dir, output_dir):
    """
    For each image in input directory, detect spots using particle filtering
    to register fiducial spots to blobs detected in the image.

    :param str input_dir: Input directory containing images and an xml file
        with parameters
    :param str output_dir: Directory where output is written to
    """

    MetaData(input_dir, output_dir)

    xlwriter_od_well = pd.ExcelWriter(
        os.path.join(constants.RUN_PATH, 'median_ODs_per_well.xlsx'),
    )
    pdantigen = pd.DataFrame(constants.ANTIGEN_ARRAY)
    pdantigen.to_excel(xlwriter_od_well, sheet_name='antigens')

    # Initialize background estimator
    bg_estimator = background_estimator.BackgroundEstimator2D(
        block_size=128,
        order=2,
        normalize=False,
    )

    # Get grid rows and columns from params
    nbr_grid_rows = constants.params['rows']
    nbr_grid_cols = constants.params['columns']
    fiducials_idx = constants.FIDUCIALS_IDX
    # Create spot detector instance
    spot_detector = img_processing.SpotDetector(
        imaging_params=constants.params,
    )

    # ================
    # loop over well images
    # ================
    well_images = io_utils.get_image_paths(input_dir)

    for well_name, im_path in well_images.items():
        start_time = time.time()
        image = io_utils.read_gray_im(im_path)
        # Get max intensity
        max_intensity = np.iinfo(image.dtype).max
        # Crop image to well only
        try:
            well_center, well_radi, well_mask = image_parser.find_well_border(
                image,
                detmethod='region',
                segmethod='otsu',
            )
            im_well, _ = img_processing.crop_image_at_center(
                image,
                well_center,
                2 * well_radi,
                2 * well_radi,
            )
        except IndexError:
            warnings.warn("Couldn't find well in {}".format(well_name))
            im_well = image
        # Find spot center coordinates
        spot_coords = spot_detector.get_spot_coords(im_well)

        # Initial estimate of spot center
        center_point = tuple((im_well.shape[0] / 2, im_well.shape[1] / 2))
        grid_coords = registration.create_reference_grid(
            center_point=center_point,
            nbr_grid_rows=nbr_grid_rows,
            nbr_grid_cols=nbr_grid_cols,
            spot_dist=constants.SPOT_DIST_PIX,
        )
        fiducial_coords = grid_coords[fiducials_idx, :]

        # Use particle filter to register fiducials to detected spots
        reg_ok = True
        particles = registration.create_gaussian_particles(
            stds=np.array(constants.STDS),
            nbr_particles=4000,
        )
        t_matrix, min_dist = registration.particle_filter(
            fiducial_coords=fiducial_coords,
            spot_coords=spot_coords,
            particles=particles,
            stds=np.array(constants.STDS),
            debug=constants.DEBUG,
        )
        if min_dist > constants.REG_DIST_THRESH:
            warnings.warn("Registration failed, repeat with outlier removal")
            t_matrix, min_dist = registration.particle_filter(
                fiducial_coords=fiducial_coords,
                spot_coords=spot_coords,
                particles=particles,
                stds=np.array(constants.STDS),
                remove_outlier=True,
                debug=constants.DEBUG,
            )
            # Warn if fit is still bad
            if min_dist > constants.REG_DIST_THRESH:
                reg_ok = False

        if not reg_ok:
            warnings.warn("Final registration failed,"
                          "will not write OD for {}".format(well_name))
            if constants.DEBUG:
                reg_coords = np.squeeze(cv.transform(np.array([grid_coords]), t_matrix))
                debug_plots.plot_registration(
                    im_well,
                    spot_coords,
                    grid_coords[fiducials_idx, :],
                    reg_coords,
                    os.path.join(constants.RUN_PATH, well_name + '_failed'),
                    max_intensity=max_intensity,
                )
            continue

        # Transform grid coordinates
        reg_coords = np.squeeze(cv.transform(np.array([grid_coords]), t_matrix))
        # Crop image
        im_crop, crop_coords = img_processing.crop_image_from_coords(
            im=im_well,
            grid_coords=reg_coords,
        )
        im_crop = im_crop / max_intensity
        # Estimate background
        background = bg_estimator.get_background(im_crop)
        # Find spots near grid locations
        props_placed_by_loc, bgprops_by_loc = array_gen.get_spot_intensity(
            coords=crop_coords,
            im_int=im_crop,
            background=background,
            params=constants.params,
        )
        # Create arrays and assign properties
        props_array = txt_parser.create_array(
            constants.params['rows'],
            constants.params['columns'],
            dtype=object,
        )
        bgprops_array = txt_parser.create_array(
            constants.params['rows'],
            constants.params['columns'],
            dtype=object,
        )
        props_array_placed = image_parser.assign_props_to_array(
            props_array,
            props_placed_by_loc,
        )
        bgprops_array = image_parser.assign_props_to_array(
            bgprops_array,
            bgprops_by_loc,
        )
        od_well, int_well, bg_well = image_parser.compute_od(
            props_array_placed,
            bgprops_array,
        )

        # populate 96-well plate constants with OD, INT, BG arrays
        report.write_od_to_plate(od_well, well_name, 'od')
        report.write_od_to_plate(int_well, well_name, 'int')
        report.write_od_to_plate(bg_well, well_name, 'bg')

        # Write ODs per well
        pd_od = pd.DataFrame(od_well)
        pd_od.to_excel(xlwriter_od_well, sheet_name=well_name)

        print("Time to register grid to {}: {:.3f} s".format(
            well_name,
            time.time() - start_time),
        )

        # ==================================
        # SAVE FOR DEBUGGING
        if constants.DEBUG:
            start_time = time.time()
            # Save spot and background intensities
            output_name = os.path.join(constants.RUN_PATH, well_name)
            # Save OD plots, composite spots and registration
            debug_plots.plot_od(
                od_well,
                int_well,
                bg_well,
                output_name,
            )
            debug_plots.save_composite_spots(
                im_crop,
                props_array_placed,
                output_name,
                from_source=True,
            )
            debug_plots.plot_background_overlay(
                im_crop,
                background,
                output_name,
            )
            debug_plots.plot_registration(
                im_well,
                spot_coords,
                grid_coords[fiducials_idx, :],
                reg_coords,
                output_name,
                max_intensity=max_intensity,
            )
            print("Time to save debug images: {:.3f} s".format(
                time.time() - start_time),
            )

    xlwriter_od = pd.ExcelWriter(
        os.path.join(constants.RUN_PATH, 'median_ODs.xlsx'),
    )
    xlwriter_int = pd.ExcelWriter(
        os.path.join(constants.RUN_PATH, 'median_intensities.xlsx'),
    )
    xlwriter_bg = pd.ExcelWriter(
        os.path.join(constants.RUN_PATH, 'median_backgrounds.xlsx'),
    )

    report.write_antigen_report(xlwriter_od, 'od')
    report.write_antigen_report(xlwriter_int, 'int')
    report.write_antigen_report(xlwriter_bg, 'bg')

    xlwriter_od.close()
    xlwriter_int.close()
    xlwriter_bg.close()
    xlwriter_od_well.close()
