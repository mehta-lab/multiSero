# bchhun, {2020-04-02}

import os

import array_analyzer.extract.image_parser as image_parser
import array_analyzer.extract.txt_parser as txt_parser
import array_analyzer.extract.img_processing as img_processing
import array_analyzer.load.debug_plots as debug_plots
import array_analyzer.extract.constants as c
from array_analyzer.extract.metadata import MetaData

import time
import skimage.io as io
import pandas as pd
import numpy as np
import re


def interp(input_folder_, output_folder_, method='interp', debug=False):
    """

    :param input_folder_:
    :param output_folder_:
    :param method:
    :param debug:
    :return:
    """

    MetaData(input_folder_, output_folder_)

    os.makedirs(c.RUN_PATH, exist_ok=True)

    # Write an excel file that can be read into jupyter notebook with minimal parsing.
    xl_writer_od = pd.ExcelWriter(os.path.join(c.RUN_PATH, 'ODs.xlsx'))
    pdantigen = pd.DataFrame(c.ANTIGEN_ARRAY)
    pdantigen.to_excel(xl_writer_od, sheet_name='antigens')

    if debug:
        xlwriter_int = pd.ExcelWriter(os.path.join(c.RUN_PATH, 'intensities.xlsx'))
        xlwriter_bg = pd.ExcelWriter(os.path.join(c.RUN_PATH, 'backgrounds.xlsx'))

    # ================
    # loop over images => good place for multiproc?  careful with columns in report
    # ================
    images = [file for file in os.listdir(input_folder_) if '.png' in file or '.tif' in file or '.jpg' in file]

    # remove any images that are not images of wells.
    wellimages = [file for file in images if re.match(r'[A-P][0-9]{1,2}', file)]

    # sort by letter, then by number (with '10' coming AFTER '9')
    wellimages.sort(key=lambda x: (x[0], int(x[1:-4])))

    # wellimages = ['H10.png','H11.png','H12.png']
    # wellimages = ['B8.png', 'B9.png', 'B10.png']
    # wellimages = ['A12.png', 'A11.png', 'A8.png', 'A1.png']
    # wellimages = ['A9.png']
    # wellimages = ['E5.png']
    well_path = None
    output_name = None

    for well in wellimages:
        start = time.time()
        image, image_name = image_parser.read_to_grey(input_folder_, well)
        print(image_name)

        spot_props_array = txt_parser.create_array(c.params['rows'], c.params['columns'], dtype=object)
        bgprops_array = txt_parser.create_array(c.params['rows'], c.params['columns'], dtype=object)

        # finding center of well and cropping
        cx, cy, r, well_mask = image_parser.find_well_border(image, detmethod='region', segmethod='otsu')
        im_crop = image_parser.crop_image(image, cx, cy, r, border_=0)

        # find center of spots from crop
        spot_mask = image_parser.thresh_and_binarize(im_crop, method='bright_spots')
        background = img_processing.get_background(im_crop, fit_order=2)


        spot_props = image_parser.generate_props(spot_mask, intensity_image_=im_crop)

        if method == 'fit':
            spot_props = image_parser.select_props(spot_props, attribute="area", condition="greater_than", condition_value=300)
            fiducial_locations = [(0, 0), (0, 1), (0, 5), (7, 0), (7, 5)]
            pix_size = 0.0049 # in mm
            props_by_loc = image_parser.find_fiducials_markers(spot_props,
                                                  fiducial_locations,
                                                  c.params['rows'],
                                                  c.params['columns'],
                                                  c.params['v_pitch'],
                                                  c.params['h_pitch'],
                                                  im_crop.shape,
                                                  pix_size)

            spot_props_array = image_parser.assign_props_to_array_2(spot_props_array, props_by_loc)

            # use the spot_props_array to find fiducials, create a new spot_mask "placed" on the array
            placed_spotmask = image_parser.build_and_place_block_array(spot_props_array, spot_mask, c.params, return_type='region')

            spot_props = image_parser.generate_props(placed_spotmask, intensity_image_=im_crop)
            bg_props = image_parser.generate_props(placed_spotmask, intensity_image_=background)

            spot_labels = [p.label for p in spot_props]
            bg_props = image_parser.select_props(bg_props, attribute="label", condition="is_in", condition_value=spot_labels)

            props_placed_by_loc = image_parser.generate_props_dict(spot_props,
                                                      c.params['rows'],
                                                      c.params['columns'],
                                                      min_area=100)
            bgprops_by_loc = image_parser.generate_props_dict(bg_props,
                                                 c.params['rows'],
                                                 c.params['columns'],
                                                 min_area=100)
        elif method == 'interp':
            bg_props = image_parser.generate_props(spot_mask, intensity_image_=background)
            eccentricities = np.array([prop.eccentricity for prop in spot_props])
            eccent_ub = eccentricities.mean() + 2 * eccentricities.std()
            # spot_props = select_props(spot_props, attribute="area", condition="greater_than", condition_value=300)
            spot_props = image_parser.select_props(spot_props, attribute="eccentricity", condition="less_than",
                                      condition_value=eccent_ub)
            spot_labels = [p.label for p in spot_props]
            bg_props = image_parser.select_props(bg_props, attribute="label", condition="is_in", condition_value=spot_labels)

            props_placed_by_loc = image_parser.grid_from_centroids(spot_props,
                                               im_crop,
                                               c.params['rows'],
                                               c.params['columns'],
                                               )
            # This call to generate_props_dict is excessive.
            # Both spot_props and bgprops can be assigned locations in previous call.

            bgprops_by_loc = image_parser.grid_from_centroids(bg_props,
                                                 background,
                                                 c.params['rows'],
                                                 c.params['columns'],
                                                 )
        props_array_placed = image_parser.assign_props_to_array(spot_props_array, props_placed_by_loc)
        bgprops_array = image_parser.assign_props_to_array(bgprops_array, bgprops_by_loc)

        # todo: further calculations using bgprops, spot_props here
        # TODO: compute spot and background intensities,
        #  and then show them on a plate like graphic (visualize_elisa_spots).
        od_well, int_well, bg_well = image_parser.compute_od(props_array_placed, bgprops_array)

        pd_OD = pd.DataFrame(od_well)
        pd_OD.to_excel(xl_writer_od, sheet_name=image_name[:-4])

        stop = time.time()
        print(f"\ttime to process={stop-start}")

        # SAVE FOR DEBUGGING
        if debug:
            well_path = os.path.join(c.RUN_PATH)
            os.makedirs(c.RUN_PATH, exist_ok=True)
            output_name = os.path.join(well_path, image_name[:-4])

            # Save spot and background intensities.
            pd_int = pd.DataFrame(int_well)
            pd_int.to_excel(xlwriter_int, sheet_name=image_name[:-4])
            pd_bg = pd.DataFrame(bg_well)
            pd_bg.to_excel(xlwriter_bg, sheet_name=image_name[:-4])

            # Save mask of the well, cropped grayscale image, cropped spot segmentation.
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

            # # This plot shows which spots have been assigned what index.
            debug_plots.plot_spot_assignment(od_well, int_well, bg_well,
                                  im_crop, props_placed_by_loc, bgprops_by_loc,
                                  image_name, output_name, c.params)

            #   save spots
            # save_all_wells(spot_props_array, spot_ids, well_path, image_name[:-4])

            #   save a composite of all spots, where spots are from source or from region prop
            debug_plots.save_composite_spots(im_crop, props_array_placed, well_path, image_name[:-4], from_source=True)

            stop2 = time.time()
            print(f"\ttime to save debug={stop2-stop}")
    if debug:
        xlwriter_int.close()
        xlwriter_bg.close()
    xl_writer_od.close()
