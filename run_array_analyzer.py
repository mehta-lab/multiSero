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
import sys, getopt

from array_analyzer.extract.image_parser import *
from array_analyzer.extract.txt_parser import *
from array_analyzer.load.xlsx_report import *
from array_analyzer.extract.img_processing import *

import time
from datetime import datetime
import skimage.io as io
import matplotlib.pyplot as plt
import pandas as pd

def main(argv):
    inputfolder = ''
    outputfolder = ''
    debug = False
    try:
        options, remainder = getopt.getopt(argv, "hi:o:d", ["help","ifile=", "ofile=", "debug="])
    except getopt.GetoptError:
        print('run_array_analyzer.py -i <inputfolder> -o <outputfolder>')
        sys.exit(2)

    for opt, arg in options:
        if opt == '-h':
            print('run_array_analyzer.py -i <inputfolder> -o <outputfolder>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfolder = arg
        elif opt in ("-o", "--ofile"):
            outputfolder = arg
        elif opt in ("-d", "--debug"):
            print('debug mode on, saving well and spot images')
            debug = True

    if not os.path.isdir(inputfolder):
        raise ValueError("input folder is not a folder or not supplied")

    if not os.path.isdir(outputfolder):
        os.makedirs(outputfolder)

    workflow(inputfolder, outputfolder, debug=debug)


def workflow(input_folder_, output_folder_, debug=False):

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

    xlsx_workbook = create_base_template()


    # save a sub path for this processing run
    run_path = output_folder_ + os.sep + f'run_{datetime.now().hour}_{datetime.now().minute}_{datetime.now().second}'

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
    #TODO: select wells based to analyze based on user input (Bryant)
    #wellimages=['A1.png','A2.png']

    for well in wellimages:
        image, image_name = read_to_grey(input_folder_,well)
        start = time.time()
        print(image_name)
        props_array = create_array(params['rows'], params['columns'], dtype=object)
        bgprops_array = create_array(params['rows'], params['columns'], dtype=object)


        # finding center of well and cropping
        cx, cy, r, well_mask = find_well_border(image, detmethod='region', segmethod='bimodal')
        im_crop = crop_image(image, cx, cy, r, border_=0)

        # find center of spots from crop
        spot_mask = thresh_and_binarize(im_crop, method='rosin')
        # TODO: Fit a grid to identify spots (Bryant, Syuan-Ming)

        background = get_background(im_crop, fit_order=2)
        props = generate_props(spot_mask, intensity_image_=im_crop)
        bg_props = generate_props(spot_mask, intensity_image_=background)

        props = select_props(props, attribute="area", condition="greater_than", condition_value=200)
        props = select_props(props, attribute="eccentricity", condition="less_than", condition_value=0.5)
        spot_labels = [p.label for p in props]
        bg_props = select_props(bg_props, attribute="label", condition="is_in", condition_value=spot_labels)

        props_by_loc = generate_props_dict(props,
                                           params['rows'],
                                           params['columns'],
                                           min_area=100)

        # This call to generate_props_dict is excessive.
        # Both props and bgprops can be assigned locations in previous call.
        bgprops_by_loc  = generate_props_dict(bg_props,
                                           params['rows'],
                                           params['columns'],
                                           min_area=100)


        props_array = assign_props_to_array(props_array, props_by_loc)
        bgprops_array = assign_props_to_array(bgprops_array, bgprops_by_loc)

       # TODO: compute spot and background intensities, and then show them on a plate like graphic (visualize_elisa_spots).
        od_well,i_well,bg_well=compute_od(props_array,bgprops_array)

        pd_OD = pd.DataFrame(od_well)
        pd_OD.to_excel(xlwriterOD, sheet_name=image_name[:-4])

        # Add a sheet to excel file for this well.

        # xlsx report generation
        # xlsx_workbook = populate_main_tab(xlsx_workbook, spot_ids, props_array, image_name[:-4])
        # xlsx_workbook = populate_main_replicates(xlsx_workbook, props_array, antigen_array, image_name[:-4])


        stop = time.time()
        print(f"\ttime to process={stop-start}")

        # SAVE FOR DEBUGGING
        if debug:
            well_path = os.path.join(run_path)
            os.makedirs(well_path, exist_ok=True)
            output_name = os.path.join(well_path, image_name[:-4])
            im_bg_overlay = np.stack([background,im_crop,background], axis=2)

            #   save cropped image and the binary
            io.imsave(output_name + "_crop.png",
                      (255*im_crop).astype('uint8'))
            io.imsave(output_name + "_crop_binary.png",
                      (255 * spot_mask).astype('uint8'))
            io.imsave(output_name + "_well_mask.png",
                      (255 * well_mask).astype('uint8'))
            io.imsave(output_name + "_crop_bg_overlay.png",
                      (255 * im_bg_overlay).astype('uint8'))

            # This plot shows which spots have been assigned what index.
            plt.imshow(im_crop)
            plt.colorbar()
            for r in np.arange(params['rows']):
                for c in np.arange(params['columns']):
                    try:
                        ceny, cenx = props_by_loc[(r,c)].centroid
                    except:
                        spot_text = '(' + str(r) + ',' + str(c) + ')'
                        print(spot_text+ 'not found')
                    else:
                        cenybg, cenxbg = bgprops_by_loc[(r,c)].centroid
                        plt.plot(cenx,ceny,'w+')
                        plt.plot(cenxbg,cenybg,'wx')
                        spot_text='(' + str(r) + ',' + str(c) + ')'
                        plt.text(cenx,ceny-5,spot_text, va='bottom', ha='center', color='w')
                        plt.text(0,0,image_name[:-4]+',spot count='+str(len(props_by_loc)))
            figcentroid=plt.gcf()
            centroids_debug=output_name+'_overlayCentroids.png'
            figcentroid.savefig(centroids_debug)
            plt.show()

            plt.figure(figsize=(6,1.5))
            plt.subplot(131)
            plt.imshow(i_well,cmap='gray')
            plt.colorbar()
            plt.title('intensity')

            plt.subplot(132)
            plt.imshow(bg_well, cmap='gray')
            plt.colorbar()
            plt.title('background')

            plt.subplot(133)
            plt.imshow(od_well, cmap='gray')
            plt.colorbar()
            plt.title('OD')

            figOD = plt.gcf()
            od_debug = output_name + '_od.png'
            figOD.savefig(od_debug)
            plt.show()




            #   save spots
            # for row in range(props_array.shape[0]):
            #     for col in range(props_array.shape[1]):
            #
            #         cell = spot_ids[row][col]
            #         if cell == '':
            #             continue
            #
            #         prop = props_array[row][col]
            #         if prop is not None:
            #             io.imsave(well_path + os.sep + image_name[:-4] + f"_spot_{cell}.png",
            #                       (255*prop.intensity_image).astype('uint8'))
            #         else:
            #             io.imsave(well_path + os.sep + image_name[:-4] + f"_spot_{cell}.png",
            #                       (255*np.ones((32, 32)).astype('uint8')))

    # SAVE COMPLETED WORKBOOK
    # xlsx_workbook.save(run_path + os.sep +
    #                   f'testrun_{datetime.now().year}_'
                       f'{datetime.now().month}{datetime.now().day}_'
                       f'{datetime.now().hour}{datetime.now().minute}.xlsx')

    xlwriterOD.close()

if __name__ == "__main__":
    input_path = '/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/Plates_given_to_manu/2020-01-15_plate4_AEP_Feb3_6mousesera'
    output_path = '/Users/shalin.mehta/Documents/images_local/2020-01-15_plate4_AEP_Feb3_6mousesera/'

    #path = '/Users/bryant.chhun/PycharmProjects/array-imager/Plates_given_to_manu/2020-01-15_plate4_AEP_Feb3_6mousesera'
    # path = '/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/Plates_given_to_manu/2020-01-15_plate4_AEP_Feb3_6mousesera'

    input = ['-i', input_path, '-o', output_path, '-d']

    # main(sys.argv[1:])
    main(input)
