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

    # ================
    # loop over well images
    # ================
    well_images = io_utils.get_image_paths(input_dir)
    well_failed = []
    # well_names = ['A6', 'A11', 'D12', 'E1', 'E6', 'F1', 'F6']
    # well_images = {well: well_images[well] for well in well_names}
    for well_name, im_path in well_images.items():
        start_time = time.time()
        image = io_utils.read_gray_im(im_path)

        nbr_expected_spots = nbr_grid_cols * nbr_grid_rows
        spot_coords, im_filtered = img_processing.get_spot_coords(
            image,
            min_area=500,
            min_thresh=100,
            max_thresh=255,
            min_circularity=0.1,
            min_convexity=0.5,
            min_dist_between_blobs=10,
            min_repeatability=2,
            nbr_expected_spots=nbr_expected_spots,
        )

        # Initial estimate of spot center
        mean_point = tuple(np.mean(spot_coords, axis=0))
        grid_coords = registration.create_reference_grid(
            mean_point=mean_point,
            nbr_grid_rows=nbr_grid_rows,
            nbr_grid_cols=nbr_grid_cols,
            spot_dist=constants.SPOT_DIST_PIX,
        )
        fiducial_coords = grid_coords[fiducials_idx, :]

        particles = registration.create_gaussian_particles(
            mean_point=(0, 0),
            stds=np.array(constants.STDS),
            scale_mean=1.,
            angle_mean=0.,
            nbr_particles=4000,
        )
        # Optimize estimated coordinates with iterative closest point
        t_matrix, min_dist = registration.particle_filter(
            fiducial_coords=fiducial_coords,
            spot_coords=spot_coords,
            particles=particles,
            stds=np.array(constants.STDS),
            stop_criteria=0.1,
            debug=constants.DEBUG,
        )

        # Warn if fit is still bad
        if min_dist > constants.REG_DIST_THRESH:
            spot_coords, im_filtered = img_processing.get_spot_coords(
                image,
                min_area=200,
                min_thresh=100,
                max_thresh=255,
                min_circularity=0.1,
                min_convexity=0.5,
                min_dist_between_blobs=10,
                min_repeatability=2,
                sigma_gauss=10,
                nbr_expected_spots=nbr_expected_spots,
            )
            t_matrix, min_dist = registration.particle_filter(
                fiducial_coords=fiducial_coords,
                spot_coords=spot_coords,
                particles=particles,
                stds=np.array(constants.STDS),
                stop_criteria=0.1,
                debug=constants.DEBUG,
            )
            if min_dist > constants.REG_DIST_THRESH:
                warnings.warn("Registration failed, repeat with outlier removal")
                t_matrix, min_dist = registration.particle_filter(
                    fiducial_coords=fiducial_coords,
                    spot_coords=spot_coords,
                    particles=particles,
                    stds=np.array(constants.STDS),
                    stop_criteria=0.1,
                    remove_outlier=True,
                    debug=constants.DEBUG,
                )
                if min_dist > constants.REG_DIST_THRESH:
                    warnings.warn("Final registration failed, registration may be flawed")
                    well_failed.append(well_name)

        # Transform grid coordinates
        reg_coords = np.squeeze(cv.transform(np.array([grid_coords]), t_matrix))
        # Crop image
        im_crop, crop_coords = img_processing.crop_image_from_coords(
            im=image,
            grid_coords=reg_coords,
        )
        im_crop = im_crop / np.iinfo(im_crop.dtype).max
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
                image,
                spot_coords,
                grid_coords[fiducials_idx, :],
                reg_coords,
                output_name,
            )
            im_filtered_name = '_'.join([output_name, 'LoG_filtered.png'])
            cv.imwrite(im_filtered_name, im_filtered)
            print("Time to save debug images: {:.3f} s".format(
                time.time() - start_time),
            )
    print('Registration for well {} failed.'.format(', '.join(well_failed)))
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
