# bchhun, {2020-04-02}

from array_analyzer.extract.txt_parser import *
from array_analyzer.extract.img_processing import *
from array_analyzer.load.debug_plots import *
from array_analyzer.transform.property_filters import *
import array_analyzer.transform.array_generation as array_gen
import re
import time
from datetime import datetime
import skimage.io as io
import pandas as pd

# SCENION_SPOT_DIST = 82

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
    xlwriterOD = pd.ExcelWriter(os.path.join(run_path, 'python_median_ODs.xlsx'))
    if debug:
        xlwriter_int = pd.ExcelWriter(os.path.join(run_path, 'python_median_intensities.xlsx'))
        xlwriter_bg = pd.ExcelWriter(os.path.join(run_path, 'python_median_backgrounds.xlsx'))

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
    # wellimages = ['B8.png', 'D5.png']
    # wellimages = ['A12.png', 'A11.png', 'A8.png', 'A1.png']
    # wellimages = ['A4.png']
    # wellimages = ['E11.png']
    well_path = None
    output_name = None

    for well in wellimages:
        start = time.time()
        image, image_name = read_to_grey(input_folder_, well)
        print(image_name)

        spot_props_array = create_array(params['rows'], params['columns'], dtype=object)
        bgprops_array = create_array(params['rows'], params['columns'], dtype=object)

        # finding center of well and cropping
        well_center, well_radi, well_mask = find_well_border(image, detmethod='region', segmethod='otsu')
        im_crop, _ = \
            crop_image_at_center(image, well_center, 2 * well_radi, 2 * well_radi)

        # find center of spots from crop
        spot_mask = thresh_and_binarize(im_crop, method='bright_spots')
        background = get_background(im_crop, fit_order=2)


        spot_props = generate_props(spot_mask, intensity_image_=im_crop)

        if debug:
            well_path = os.path.join(run_path)
            os.makedirs(run_path, exist_ok=True)
            output_name = os.path.join(well_path, image_name[:-4])

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

        bg_props = generate_props(spot_mask, intensity_image_=background)
        # eccentricities = np.array([prop.eccentricity for prop in spot_props])
        # eccent_ub = eccentricities.mean() + 2.5 * eccentricities.std()
        # spot_props = select_props(spot_props, attribute="area", condition="greater_than", condition_value=300)
        # spot_props = select_props(spot_props, attribute="eccentricity", condition="less_than",
        #                           condition_value=eccent_ub)
        spot_labels = [p.label for p in spot_props]
        bg_props = select_props(bg_props, attribute="label", condition="is_in", condition_value=spot_labels)

        crop_coords = \
            grid_from_centroids(spot_props,
                               im_crop,
                               background,
                               params['rows'],
                               params['columns'],
                               )
        props_by_loc, bgprops_by_loc = \
            array_gen.get_spot_intensity(
                coords=crop_coords,
                im_int=im_crop,
                background=background,
                params=params
            )
        props_array_placed = assign_props_to_array(spot_props_array, props_by_loc)
        bgprops_array = assign_props_to_array(bgprops_array, bgprops_by_loc)

        od_well, int_well, bg_well = compute_od(props_array_placed, bgprops_array)

        pd_OD = pd.DataFrame(od_well)
        pd_OD.to_excel(xlwriterOD, sheet_name=image_name[:-4])

        stop = time.time()
        print(f"\ttime to process={stop-start}")

        # SAVE FOR DEBUGGING
        if debug:
            # Save spot and background intensities.
            pd_int = pd.DataFrame(int_well)
            pd_int.to_excel(xlwriter_int, sheet_name=image_name[:-4])
            pd_bg = pd.DataFrame(bg_well)
            pd_bg.to_excel(xlwriter_bg, sheet_name=image_name[:-4])

            # # This plot shows which spots have been assigned what index.
            plot_centroid_overlay(
                im_crop,
                params,
                props_by_loc,
                bgprops_by_loc,
                output_name,
            )
            plot_od(
                od_well,
                int_well,
                bg_well,
                output_name,
            )
            # #   save a composite of all spots, where spots are from source or from region prop
            save_composite_spots(im_crop, props_array_placed, output_name, from_source=True)
            save_composite_spots(im_crop, props_array_placed, output_name, from_source=False)

            stop2 = time.time()
            print(f"\ttime to save debug={stop2-stop}")
    if debug:
        xlwriter_int.close()
        xlwriter_bg.close()
    xlwriterOD.close()
