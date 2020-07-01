import cv2 as cv
import logging
import numpy as np
import os
import pandas as pd
import time

import array_analyzer.extract.background_estimator as background_estimator
import array_analyzer.extract.image_parser as image_parser
import array_analyzer.extract.img_processing as img_processing
import array_analyzer.extract.metadata as metadata
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
    logger = logging.getLogger(constants.LOG_NAME)

    metadata.MetaData(input_dir, output_dir)
    nbr_outliers = constants.params['nbr_outliers']

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
    nbr_fiducials = len(fiducials_idx)
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
        logger.info("Extracting well: {}".format(well_name))
        # Get max intensity
        max_intensity = np.iinfo(image.dtype).max
        logger.debug("Image max intensity: {}".format(max_intensity))
        # Crop image to well only
        try:
            well_center, well_radi, well_mask = image_parser.find_well_border(
                image,
                detmethod='region',
                segmethod='otsu',
            )
            im_well, _ = img_processing.crop_image_at_center(
                im=image,
                center=well_center,
                height=2 * well_radi,
                width=2 * well_radi,
            )
        except IndexError:
            logging.warning("Couldn't find well in {}".format(well_name))
            im_well = image

        # Find spot center coordinates
        spot_coords = spot_detector.get_spot_coords(
            im=im_well,
            max_intensity=max_intensity,
        )
        if spot_coords.shape[0] < constants.MIN_NBR_SPOTS:
            logging.warning("Not enough spots detected in {},"
                            "continuing.".format(well_name))
            continue

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
            nbr_particles=constants.NBR_PARTICLES,
        )
        t_matrix, min_dist = registration.particle_filter(
            fiducial_coords=fiducial_coords,
            spot_coords=spot_coords,
            particles=particles,
            stds=np.array(constants.STDS),
            logger=logger,
        )
        reg_dist = min_dist / nbr_fiducials
        if reg_dist > constants.REG_DIST_THRESH:
            logger.warning("Registration failed for {}, "
                           "repeat with outlier removal".format(well_name))
            t_matrix, min_dist = registration.particle_filter(
                fiducial_coords=fiducial_coords,
                spot_coords=spot_coords,
                particles=particles,
                stds=np.array(constants.STDS),
                logger=logger,
                nbr_outliers=nbr_outliers,
            )
            # Warn if fit is still bad
            reg_dist = min_dist / (nbr_fiducials - nbr_outliers)
            if reg_dist > constants.REG_DIST_THRESH:
                reg_ok = False
        logger.info("Particle filter min dist: {}".format(reg_dist))

        # Transform grid coordinates
        reg_coords = np.squeeze(cv.transform(np.array([grid_coords]), t_matrix))
        # Check that registered coordinates are inside well
        reg_ok = registration.check_reg_coords(
            reg_coords=reg_coords,
            im_shape=im_well.shape,
            reg_ok=reg_ok,
        )
        if not reg_ok:
            logger.warning("Final registration failed,"
                           "will not write OD for {}".format(well_name))
            if constants.DEBUG:
                debug_plots.plot_registration(
                    im_well,
                    spot_coords,
                    grid_coords[fiducials_idx, :],
                    reg_coords,
                    os.path.join(constants.RUN_PATH, well_name + '_failed'),
                    max_intensity=max_intensity,
                )
            continue

        # Crop image
        im_crop, crop_coords = img_processing.crop_image_from_coords(
            im=im_well,
            grid_coords=reg_coords,
        )
        im_crop = im_crop / max_intensity
        # Estimate background
        background = bg_estimator.get_background(im_crop)
        # Find spots near grid locations and compute properties
        spot_props, bg_props = array_gen.get_spot_intensity(
            coords=crop_coords,
            im=im_crop,
            background=background,
            params=constants.params,
        )
        od_well, int_well, bg_well = image_parser.compute_od(
            spot_props,
            bg_props,
        )
        # populate 96-well plate constants with OD, INT, BG arrays
        report.write_od_to_plate(od_well, well_name, constants.WELL_OD_ARRAY)
        report.write_od_to_plate(int_well, well_name, constants.WELL_INT_ARRAY)
        report.write_od_to_plate(bg_well, well_name, constants.WELL_BG_ARRAY)

        # Write ODs per well
        pd_od = pd.DataFrame(od_well)
        pd_od.to_excel(xlwriter_od_well, sheet_name=well_name)

        time_msg = "Time to extract OD in {}: {:.3f} s".format(
            well_name,
            time.time() - start_time,
        )
        print(time_msg)
        logger.info(time_msg)

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
                spot_props=spot_props,
                output_name=output_name,
                image=im_crop,
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
            logger.debug("Time to save debug images: {:.3f} s".format(
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

    report.write_antigen_report(
        xlwriter_od,
        constants.WELL_OD_ARRAY,
        array_type='od',
        logger=logger,
    )
    report.write_antigen_report(
        xlwriter_int,
        constants.WELL_INT_ARRAY,
        array_type='int',
        logger=logger,
    )
    report.write_antigen_report(
        xlwriter_bg,
        constants.WELL_BG_ARRAY,
        array_type='bg',
        logger=logger,
    )

    xlwriter_od.close()
    xlwriter_int.close()
    xlwriter_bg.close()
    xlwriter_od_well.close()
