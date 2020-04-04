# bchhun, {2020-04-02}

import sys, getopt, os

from array_analyzer.extract.image_parser import *
from array_analyzer.extract.txt_parser import *
from array_analyzer.load.xlsx_report import *
from array_analyzer.extract.img_processing import *
from array_analyzer.load.debug_images import *
from array_analyzer.transform.property_filters import *

import time
from datetime import datetime
import skimage.io as io
import matplotlib.pyplot as plt
import pandas as pd


def interp(input_folder_, output_folder_, method='interp', debug=False):

    xml = [f for f in os.listdir(input_folder_) if '.xml' in f]
    if len(xml) > 1:
        raise IOError("more than one .xml file found, aborting")
    xml_path = input_folder_+os.sep+xml[0]

    # parsing .xml
    fiduc, spots, repl, params = create_xml_dict(xml_path)

    # creating our arrays
    spot_ids = create_array(params['rows'], params['columns'])
    antigen_array = create_array(params['rows'], params['columns'])

    # adding .xml info to these arrays
    spot_ids = populate_array_id(spot_ids, spots)
    # spot_ids = populate_array_fiduc(spot_ids, fiduc)

    antigen_array = populate_array_antigen(antigen_array, spot_ids, repl)

    # save a sub path for this processing run
    run_path = output_folder_ + os.sep + f'{datetime.now().month}_{datetime.now().day}_{datetime.now().hour}_{datetime.now().minute}_{datetime.now().second}'

    # Write an excel file that can be read into jupyter notebook with minimal parsing.
    xlwriterOD = pd.ExcelWriter(os.path.join(run_path, 'ODs.xlsx'))
    pdantigen = pd.DataFrame(antigen_array)
    pdantigen.to_excel(xlwriterOD, sheet_name='antigens')

    if not os.path.isdir(run_path):
        os.mkdir(run_path)

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

    if debug:
        well_path = os.path.join(run_path)
        os.makedirs(well_path, exist_ok=True)

    for well in wellimages:
        start = time.time()
        image, image_name = read_to_grey(input_folder_, well)
        print(image_name)

        if debug:
            output_name = os.path.join(well_path, image_name[:-4])

        spot_props_array = create_array(params['rows'], params['columns'], dtype=object)
        bgprops_array = create_array(params['rows'], params['columns'], dtype=object)

        # finding center of well and cropping
        cx, cy, r, well_mask = find_well_border(image, detmethod='region', segmethod='otsu')
        im_crop = crop_image(image, cx, cy, r, border_=0)

        # find center of spots from crop
        spot_mask = thresh_and_binarize(im_crop, method='bright_spots')
        if debug:
            io.imsave(output_name + "_well_mask.png",
                      (255 * well_mask).astype('uint8'))
            io.imsave(output_name + "_crop.png",
                      (255 * im_crop).astype('uint8'))
            io.imsave(output_name + "_crop_binary.png",
                  (255 * spot_mask).astype('uint8'))
        background = get_background(im_crop, fit_order=2)

        if debug:
            im_bg_overlay = np.stack([background, im_crop, background], axis=2)
            io.imsave(output_name + "_crop_bg_overlay.png",
                      (255 * im_bg_overlay).astype('uint8'))

        spot_props = generate_props(spot_mask, intensity_image_=im_crop)

        if method == 'fit':
            spot_props = select_props(spot_props, attribute="area", condition="greater_than", condition_value=300)
            fiducial_locations = [(0, 0), (0, 1), (0, 5), (7, 0), (7, 5)]
            pix_size = 0.0049 # in mm
            props_by_loc = find_fiducials_markers(spot_props,
                                                  fiducial_locations,
                                                  params['rows'],
                                                  params['columns'],
                                                  params['v_pitch'],
                                                  params['h_pitch'],
                                                  im_crop.shape,
                                                  pix_size)

            spot_props_array = assign_props_to_array_2(spot_props_array, props_by_loc)

            # use the spot_props_array to find fiducials, create a new spot_mask "placed" on the array
            placed_spotmask = build_and_place_block_array(spot_props_array, spot_mask, params, return_type='region')

            spot_props = generate_props(placed_spotmask, intensity_image_=im_crop)
            bg_props = generate_props(placed_spotmask, intensity_image_=background)

            spot_labels = [p.label for p in spot_props]
            bg_props = select_props(bg_props, attribute="label", condition="is_in", condition_value=spot_labels)

            props_placed_by_loc = generate_props_dict(spot_props,
                                                      params['rows'],
                                                      params['columns'],
                                                      min_area=100)
            bgprops_by_loc = generate_props_dict(bg_props,
                                                 params['rows'],
                                                 params['columns'],
                                                 min_area=100)
        elif method == 'interp':
            bg_props = generate_props(spot_mask, intensity_image_=background)
            eccentricities = np.array([prop.eccentricity for prop in spot_props])
            eccent_ub = eccentricities.mean() + 2 * eccentricities.std()
            # spot_props = select_props(spot_props, attribute="area", condition="greater_than", condition_value=300)
            spot_props = select_props(spot_props, attribute="eccentricity", condition="less_than",
                                      condition_value=eccent_ub)
            spot_labels = [p.label for p in spot_props]
            bg_props = select_props(bg_props, attribute="label", condition="is_in", condition_value=spot_labels)

            props_placed_by_loc = grid_from_centroids(spot_props,
                                               im_crop,
                                               params['rows'],
                                               params['columns'],
                                               )
            # This call to generate_props_dict is excessive.
            # Both spot_props and bgprops can be assigned locations in previous call.

            bgprops_by_loc = grid_from_centroids(bg_props,
                                                 background,
                                                 params['rows'],
                                                 params['columns'],
                                                 )
        props_array_placed = assign_props_to_array(spot_props_array, props_placed_by_loc)
        bgprops_array = assign_props_to_array(bgprops_array, bgprops_by_loc)

        # todo: further calculations using bgprops, spot_props here
        # TODO: compute spot and background intensities,
        #  and then show them on a plate like graphic (visualize_elisa_spots).
        od_well, i_well, bg_well = compute_od(props_array_placed, bgprops_array)

        pd_OD = pd.DataFrame(od_well)
        pd_OD.to_excel(xlwriterOD, sheet_name=image_name[:-4])

        stop = time.time()
        print(f"\ttime to process={stop-start}")

        # SAVE FOR DEBUGGING
        if debug:

            # This plot shows which spots have been assigned what index.
            plot_spot_assignment(od_well, i_well, bg_well,
                                 im_crop, props_placed_by_loc, bgprops_by_loc,
                                 image_name, output_name, params)

            #   save spots
            # save_all_wells(spot_props_array, spot_ids, well_path, image_name[:-4])

            #   save a composite of all spots, where spots are from source or from region prop
            save_composite_spots(im_crop, props_array_placed, well_path, image_name[:-4], from_source=True)
            save_composite_spots(im_crop, props_array_placed, well_path, image_name[:-4], from_source=False)

            stop2 = time.time()
            print(f"\ttime to save debug={stop2-stop}")

    xlwriterOD.close()
