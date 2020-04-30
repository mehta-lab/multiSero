import re
import time
import os
import numpy as np
import pandas as pd
import skimage.io as io

import array_analyzer.extract.image_parser as image_parser
import array_analyzer.extract.txt_parser as txt_parser
import array_analyzer.extract.img_processing as img_processing
import array_analyzer.load.debug_plots as debug_plots
import array_analyzer.extract.constants as c
import array_analyzer.transform.array_generation as array_gen
import array_analyzer.extract.background_estimator as background_estimator
import array_analyzer.utils.io_utils as io_utils

from array_analyzer.extract.metadata import MetaData


def interp(input_dir, output_dir, debug=False):

    # xml = [f for f in os.listdir(input_folder_) if '.xml' in f]
    # if len(xml) > 1:
    #     raise IOError("more than one .xml file found, aborting")
    # xml_path = input_folder_+os.sep+xml[0]

    # parsing .xml
    # fiduc, spots, repl, params = create_xml_dict(xml_path)

    # creating our arrays
    # spot_ids = create_array(params['rows'], params['columns'])
    # antigen_array = create_array(params['rows'], params['columns'])

    # adding .xml info to these arrays
    # spot_ids = populate_array_id(spot_ids, spots)
    # spot_ids = populate_array_fiduc(spot_ids, fiduc)

    # antigen_array = populate_array_antigen(antigen_array, spot_ids, repl)

    # save a sub path for this processing run
    # run_path = output_folder_ + os.sep + f'{datetime.now().month}_{datetime.now().day}_{datetime.now().hour}_{datetime.now().minute}_{datetime.now().second}'

    MetaData(input_dir, output_dir)

    os.makedirs(c.RUN_PATH, exist_ok=True)

    # Write an excel file that can be read into jupyter notebook with minimal parsing.
    xlwriter_od = pd.ExcelWriter(os.path.join(c.RUN_PATH, 'python_median_ODs.xlsx'))
    if debug:
        xlwriter_int = pd.ExcelWriter(
            os.path.join(c.RUN_PATH, 'python_median_intensities.xlsx')
        )
        xlwriter_bg = pd.ExcelWriter(
            os.path.join(c.RUN_PATH, 'python_median_backgrounds.xlsx')
        )

    # Initialize background estimator
    bg_estimator = background_estimator.BackgroundEstimator2D(
        block_size=128,
        order=2,
        normalize=False,
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

        spot_props_array = txt_parser.create_array(c.params['rows'], c.params['columns'], dtype=object)
        bgprops_array = txt_parser.create_array(c.params['rows'], c.params['columns'], dtype=object)

        # finding center of well and cropping
        well_center, well_radi, well_mask = image_parser.find_well_border(image, detmethod='region', segmethod='otsu')
        im_crop, _ = \
            img_processing.crop_image_at_center(image, well_center, 2 * well_radi, 2 * well_radi)

        # find center of spots from crop
        spot_mask = img_processing.thresh_and_binarize(im_crop, method='bright_spots')
        background = bg_estimator.get_background(im_crop)

        spot_props = image_parser.generate_props(spot_mask, intensity_image_=im_crop)

        if debug:
            output_name = os.path.join(c.RUN_PATH, well_name)

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

        bg_props = image_parser.generate_props(spot_mask, intensity_image_=background)
        # eccentricities = np.array([prop.eccentricity for prop in spot_props])
        # eccent_ub = eccentricities.mean() + 2.5 * eccentricities.std()
        # spot_props = select_props(spot_props, attribute="area", condition="greater_than", condition_value=300)
        # spot_props = select_props(spot_props, attribute="eccentricity", condition="less_than",
        #                           condition_value=eccent_ub)
        spot_labels = [p.label for p in spot_props]
        bg_props = image_parser.select_props(bg_props, attribute="label", condition="is_in", condition_value=spot_labels)

        crop_coords = image_parser.grid_from_centroids(
            spot_props,
            im_crop,
            background,
            c.params['rows'],
            c.params['columns']
        )
        props_by_loc, bgprops_by_loc = array_gen.get_spot_intensity(
            coords=crop_coords,
            im_int=im_crop,
            background=background,
            params=c.params
        )
        props_array_placed = image_parser.assign_props_to_array(spot_props_array, props_by_loc)
        bgprops_array = image_parser.assign_props_to_array(bgprops_array, bgprops_by_loc)

        od_well, int_well, bg_well = image_parser.compute_od(props_array_placed, bgprops_array)

        pd_OD = pd.DataFrame(od_well)
        pd_OD.to_excel(xlwriter_od, sheet_name=well_name)

        stop = time.time()
        print(f"\ttime to process={stop-start}")

        # SAVE FOR DEBUGGING
        if debug:
            # Save spot and background intensities.
            pd_int = pd.DataFrame(int_well)
            pd_int.to_excel(xlwriter_int, sheet_name=well_name)
            pd_bg = pd.DataFrame(bg_well)
            pd_bg.to_excel(xlwriter_bg, sheet_name=well_name)

            # This plot shows which spots have been assigned what index.
            debug_plots.plot_centroid_overlay(
                im_crop,
                c.params,
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
    if debug:
        xlwriter_int.close()
        xlwriter_bg.close()
    xlwriter_od.close()
