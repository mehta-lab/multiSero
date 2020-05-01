import cv2 as cv
import numpy as np
import os
import pandas as pd
import time

import array_analyzer.extract.background_estimator as background_estimator
import array_analyzer.extract.image_parser as image_parser
import array_analyzer.extract.img_processing as img_processing
import array_analyzer.extract.txt_parser as txt_parser
from array_analyzer.extract.metadata import MetaData
import array_analyzer.extract.constants as c
import array_analyzer.load.debug_plots as debug_plots
import array_analyzer.transform.point_registration as registration
import array_analyzer.transform.array_generation as array_gen
import array_analyzer.utils.io_utils as io_utils
import array_analyzer.load.report as report


def point_registration(input_dir, output_dir, debug=False):
    """
    For each image in input directory, detect spots using particle filtering
    to register fiducial spots to blobs detected in the image.

    :param str input_dir: Input directory containing images and an xml file
        with parameters
    :param str output_dir: Directory where output is written to
    :param bool debug: For saving debug plots
    """

    MetaData(input_dir, output_dir)

    os.makedirs(c.RUN_PATH, exist_ok=True)

    xlwriter_od_well = pd.ExcelWriter(os.path.join(c.RUN_PATH, 'python_median_ODs_per_well.xlsx'))
    pdantigen = pd.DataFrame(c.ANTIGEN_ARRAY)
    pdantigen.to_excel(xlwriter_od_well, sheet_name='antigens')

    # Initialize background estimator
    bg_estimator = background_estimator.BackgroundEstimator2D(
        block_size=128,
        order=2,
        normalize=False,
    )

    # Get grid rows and columns from params
    nbr_grid_rows = c.params['rows']
    nbr_grid_cols = c.params['columns']
    fiducials_idx = c.FIDUCIALS_IDX

    # ================
    # loop over well images
    # ================
    well_images = io_utils.get_image_paths(input_dir)

    for well_name, im_path in well_images.items():
        start_time = time.time()
        image = io_utils.read_gray_im(im_path)

        spot_coords = img_processing.get_spot_coords(
            image,
            min_area=250,
            min_thresh=25,
            max_thresh=255,
            min_circularity=0,
            min_convexity=0,
            min_dist_between_blobs=10,
            min_repeatability=2,
        )

        # Initial estimate of spot center
        mean_point = tuple(np.mean(spot_coords, axis=0))
        grid_coords = registration.create_reference_grid(
            mean_point=mean_point,
            nbr_grid_rows=nbr_grid_rows,
            nbr_grid_cols=nbr_grid_cols,
            spot_dist=c.SPOT_DIST_PIX,
        )
        fiducial_coords = grid_coords[fiducials_idx, :]

        particles = registration.create_gaussian_particles(
            mean_point=(0, 0),
            stds=np.array(c.STDS),
            scale_mean=1.,
            angle_mean=0.,
            nbr_particles=1000,
        )
        # Optimize estimated coordinates with iterative closest point
        t_matrix, min_dist = registration.particle_filter(
            fiducial_coords=fiducial_coords,
            spot_coords=spot_coords,
            particles=particles,
            stds=np.array(c.STDS),
            stop_criteria=0.1,
            debug=debug,
        )
        if min_dist > c.REG_DIST_THRESH:
            print("Registration failed, repeat with outlier removal")
            t_matrix, min_dist = registration.particle_filter(
                fiducial_coords=fiducial_coords,
                spot_coords=spot_coords,
                particles=particles,
                stds=np.array(c.STDS),
                stop_criteria=0.1,
                remove_outlier=True,
                debug=debug,
            )
            # TODO: Flag bad fit if min_dist is still above threshold

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
            params=c.params,
        )
        # Create arrays and assign properties
        props_array = txt_parser.create_array(
            c.params['rows'],
            c.params['columns'],
            dtype=object,
        )
        bgprops_array = txt_parser.create_array(
            c.params['rows'],
            c.params['columns'],
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
        if debug:
            start_time = time.time()
            # Save spot and background intensities

            output_name = os.path.join(c.RUN_PATH, well_name)
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
            print("Time to save debug images: {:.3f} s".format(
                time.time() - start_time),
            )

    xlwriter_od = pd.ExcelWriter(os.path.join(c.RUN_PATH, 'python_median_ODs.xlsx'))
    xlwriter_int = pd.ExcelWriter(os.path.join(c.RUN_PATH, 'python_median_intensities.xlsx'))
    xlwriter_bg = pd.ExcelWriter(os.path.join(c.RUN_PATH, 'python_median_backgrounds.xlsx'))

    report.write_antigen_report(xlwriter_od, 'od')
    report.write_antigen_report(xlwriter_int, 'int')
    report.write_antigen_report(xlwriter_bg, 'bg')

    # loop all antigens
    # well_to_image = {v: k for k, v in c.IMAGE_TO_WELL.items()}
    # for antigen_position, antigen in np.ndenumerate(c.ANTIGEN_ARRAY):
    #     if antigen == '' or antigen is None:
    #         continue
    #     print(f"writing antigen {antigen} to excel sheets")
    #
    #     od_sheet = deepcopy(c.WELL_OUTPUT_TEMPLATE)
    #     int_sheet = deepcopy(c.WELL_OUTPUT_TEMPLATE)
    #     bg_sheet = deepcopy(c.WELL_OUTPUT_TEMPLATE)
    #
    #     # loop all wells and write OD, INT, BG of this antigen
    #     for od_position, od_well in np.ndenumerate(c.WELL_OD_ARRAY):
    #         od_val = od_well[antigen_position[0], antigen_position[1]]
    #         well_name = well_to_image[(od_position[0]+1, od_position[1]+1)]
    #         od_sheet[well_name[0]][int(well_name[1:])] = od_val
    #
    #     for int_position, int_well in np.ndenumerate(c.WELL_INT_ARRAY):
    #         int_val = int_well[antigen_position[0], antigen_position[1]]
    #         well_name = well_to_image[(int_position[0]+1, int_position[1]+1)]
    #         int_sheet[well_name[0]][int(well_name[1:])] = int_val
    #
    #     for bg_position, bg_well in np.ndenumerate(c.WELL_BG_ARRAY):
    #         bg_val = bg_well[antigen_position[0], antigen_position[1]]
    #         well_name = well_to_image[(bg_position[0]+1, bg_position[1]+1)]
    #         bg_sheet[well_name[0]][int(well_name[1:])] = bg_val
    #
    #     # write the outputs from the above three to a worksheet
    #     od_sheet_df = pd.DataFrame(od_sheet).T
    #     int_sheet_df = pd.DataFrame(int_sheet).T
    #     bg_sheet_df = pd.DataFrame(bg_sheet).T
    #
    #     od_sheet_df.to_excel(xlwriter_od, sheet_name=f'od_{antigen_position}_{antigen}')
    #     int_sheet_df.to_excel(xlwriter_int, sheet_name=f'int_{antigen_position}_{antigen}')
    #     bg_sheet_df.to_excel(xlwriter_bg, sheet_name=f'bg_{antigen_position}_{antigen}')

    xlwriter_od.close()
    xlwriter_int.close()
    xlwriter_bg.close()
