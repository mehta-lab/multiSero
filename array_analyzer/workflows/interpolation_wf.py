import time
import os
import numpy as np
import pandas as pd
import skimage.io as io

import array_analyzer.extract.image_parser as image_parser
import array_analyzer.extract.txt_parser as txt_parser
import array_analyzer.extract.img_processing as img_processing
import array_analyzer.load.debug_plots as debug_plots
import array_analyzer.extract.constants as constants
import array_analyzer.transform.array_generation as array_gen
import array_analyzer.extract.background_estimator as background_estimator
import array_analyzer.utils.io_utils as io_utils

from array_analyzer.extract.metadata import MetaData
import array_analyzer.load.report as report


def interp(input_dir, output_dir):

    MetaData(input_dir, output_dir)

    # Write an excel file that can be read into jupyter notebook with minimal parsing.
    xlwriter_od_well = pd.ExcelWriter(
        os.path.join(constants.RUN_PATH, 'median_ODs_per_well.xlsx'),
    )

    # Initialize background estimator
    bg_estimator = background_estimator.BackgroundEstimator2D(
        block_size=128,
        order=2,
        normalize=False,
    )

    # ================
    # loop over images => good place for multiproc?  careful with columns in report
    # ================
    well_images = io_utils.get_image_paths(input_dir)

    for well_name, im_path in well_images.items():
        start = time.time()
        image = io_utils.read_gray_im(im_path)

        spot_props_array = txt_parser.create_array(
            constants.params['rows'],
            constants.params['columns'],
            dtype=object,
        )
        bgprops_array = txt_parser.create_array(
            constants.params['rows'],
            constants.params['columns'],
            dtype=object,
        )

        # finding center of well and cropping
        well_center, well_radi, well_mask = image_parser.find_well_border(image, detmethod='region', segmethod='otsu')
        im_crop, _ = img_processing.crop_image_at_center(
            image,
            well_center,
            2 * well_radi,
            2 * well_radi
        )

        # find center of spots from crop
        spot_mask = img_processing.thresh_and_binarize(im_crop, method='bright_spots')
        spot_props = image_parser.generate_props(spot_mask, intensity_image_=im_crop)

        # if debug:

        crop_coords = image_parser.grid_from_centroids(
            spot_props,
            constants.params['rows'],
            constants.params['columns']
        )

        # convert to float64
        im_crop = im_crop / np.iinfo(im_crop.dtype).max
        background = bg_estimator.get_background(im_crop)
        props_by_loc, bgprops_by_loc = array_gen.get_spot_intensity(
            coords=crop_coords,
            im_int=im_crop,
            background=background,
            params=constants.params
        )
        props_array_placed = image_parser.assign_props_to_array(spot_props_array, props_by_loc)
        bgprops_array = image_parser.assign_props_to_array(bgprops_array, bgprops_by_loc)

        od_well, int_well, bg_well = image_parser.compute_od(props_array_placed, bgprops_array)

        pd_OD = pd.DataFrame(od_well)
        pd_OD.to_excel(xlwriter_od_well, sheet_name=well_name)

        # populate 96-well plate constants with OD, INT, BG arrays
        report.write_od_to_plate(od_well, well_name, 'od')
        report.write_od_to_plate(int_well, well_name, 'int')
        report.write_od_to_plate(bg_well, well_name, 'bg')

        stop = time.time()
        print(f"\ttime to process={stop-start}")

        # SAVE FOR DEBUGGING
        if constants.DEBUG:
            # Save spot and background intensities.
            output_name = os.path.join(constants.RUN_PATH, well_name)

            # # Save mask of the well, cropped grayscale image, cropped spot segmentation.
            io.imsave(output_name + "_well_mask.png",
                      (255 * well_mask).astype('uint8'))
            io.imsave(output_name + "_crop.png",
                      (255 * im_crop).astype('uint8'))
            io.imsave(output_name + "_crop_binary.png",
                      (255 * spot_mask).astype('uint8'))

            # Evaluate accuracy of background estimation with green (image), magenta (background) overlay.
            im_bg_overlay = np.stack([background, im_crop, background], axis=2)

            io.imsave(output_name + "_crop_bg_overlay.png",
                      (255 * im_bg_overlay).astype('uint8'))

            # This plot shows which spots have been assigned what index.
            debug_plots.plot_centroid_overlay(
                im_crop,
                constants.params,
                props_by_loc,
                bgprops_by_loc,
                output_name,
            )
            debug_plots.plot_od(
                od_well,
                int_well,
                bg_well,
                output_name,
            )
            # save a composite of all spots, where spots are from source or from region prop
            debug_plots.save_composite_spots(im_crop, props_array_placed, output_name, from_source=True)
            debug_plots.save_composite_spots(im_crop, props_array_placed, output_name, from_source=False)

            stop2 = time.time()
            print(f"\ttime to save debug={stop2-stop}")

    xlwriter_od_well.close()

    # create excel writers to write reports
    xlwriter_od = pd.ExcelWriter(os.path.join(constants.RUN_PATH, 'python_median_ODs.xlsx'))
    xlwriter_int = pd.ExcelWriter(os.path.join(constants.RUN_PATH, 'python_median_intensities.xlsx'))
    xlwriter_bg = pd.ExcelWriter(os.path.join(constants.RUN_PATH, 'python_median_backgrounds.xlsx'))

    report.write_antigen_report(xlwriter_od, 'od')
    report.write_antigen_report(xlwriter_int, 'int')
    report.write_antigen_report(xlwriter_bg, 'bg')

    xlwriter_od.close()
    xlwriter_int.close()
    xlwriter_bg.close()
