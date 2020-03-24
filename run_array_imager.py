# bchhun, {2020-03-23}

"""
Here we will call the methods from the ETL folders

============================
txt_parser workflow:
--------------------
1) xml_to_dict the xml file
2) create ID-array
3) create antigen-array

image_parser workflow:
----------------------
4) read_to_grey(supplied images)

# find center of well
5) thresh and binarize from 4
6) find well border from 5
7) crop image from 4

# find center of spots from crop
8) thresh and binarize from 4
9) clean spot binary from 5 (if using bimodal in 8)
10) generate props from 6
11) generate props dict from 7
12) assign props dict to array from 8

xlsx report generation workflow:
--------------------------------
13) "create base template"
14) "populate main tab" using :
        - workbook from 13
        - "ID-array" from 2
        - "props-array" from 12
        - "well" from "read_to_grey" from 4
15) "populate main replictes" :
        - workbook from 13
        - "props-array" from 12
        - "antigen-array" from 3
        - "well" from "read_to_grey" from 4
16) (populate well-tab) (in progress)
17) (populate well-replicate-tab) (in progress)

18) *repeat 14-17* using next image and well name
19) save .xlsx
==============================

==============================
FULL WORKFLOW

cli
---
- input folder
- output folder

extract
-------
A) search folder for all images, all .xmls to list
B) xlsx_report.create_base_template() step 13

C) txt_parse workflow above to create ID-array, antigen-array
D) image_parser workflow above to loop 4 (read_to_grey)
    within loop:
    E) image_parser steps 5-12

    transform
    ---------
    F) (ANY "transform" methods that will further analyze the data from E)
        (this is set aside as a place to make diagnosis calls and for downstream calculations using spot properties)

    load
    ----
    G) xlsx_report generation workflow steps 14-17

"""
import argparse

from array_imager.extract.image_parser import *
from array_imager.extract.txt_parser import *
from array_imager.load.xlsx_report import *

import time


def main():
    path = "/Users/bryant.chhun/PycharmProjects/array-imager/Plates_given_to_manu/2020-01-15_plate4_AEP_Feb3_6mousesera"
    xml = [f for f in os.listdir(path) if '.xml' in f]
    if len(xml) > 1:
        raise IOError("more than one .xml file found, aborting")

    xml_path = path+os.sep+xml[0]

    # parsing .xml
    fiduc, spots, repl, params = create_xml_dict(xml_path)

    # creating our arrays
    spot_ids = create_array(params['rows'], params['columns'])
    antigen_array = create_array(params['rows'], params['columns'])
    props_array = create_array(params['rows'], params['columns'], dtype=object)

    # adding .xml info to these arrays
    spot_ids = populate_array_id(spot_ids, spots)
    # spot_ids = populate_array_fiduc(spot_ids, fiduc)

    antigen_array = populate_array_antigen(antigen_array, spot_ids, repl)

    xlsx_workbook = create_base_template()

    start = time.time()
    for image, image_name in read_to_grey(path):
        print(image_name)
        # finding center of well and cropping
        binary = thresh_and_binarize(image, method='bimodal')
        cx, cy, r = find_well_border(binary)
        im_crop = crop_image(image, cx, cy, r, border_=200)

        # find center of spots from crop
        binary = thresh_and_binarize(im_crop, method='rosin')
        props = generate_props(binary, intensity_image_=im_crop)
        centroid_map = generate_props_dict(props, params['rows'],
                                           params['columns'],
                                           min_area=100)
        props_array = assign_props_to_array(props_array, centroid_map)

        # xlsx report generation
        xlsx_workbook = populate_main_tab(xlsx_workbook, spot_ids, props_array, image_name[:-4])
        xlsx_workbook = populate_main_replicates(xlsx_workbook, props_array, antigen_array, image_name[:-4])

        stop = time.time()
        print(f"\ttime to process={stop-start}")

    xlsx_workbook.save(path+os.sep+'testrun.xlsx')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_folder", type=str, help='path to folder containing images and .xml')
    parser.add_argument("--output_folder", type=str, help='path to output folder for intermediate images and report')
    args = parser.parse_args()

    main()
