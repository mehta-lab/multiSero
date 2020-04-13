import sys, getopt, os

from array_analyzer.extract.image_parser import *

import time
from datetime import datetime
import skimage.io as io
import pandas as pd
import string


def well_analysis(input_folder_, output_folder_, debug=False):

    start = time.time()
    # save a sub path for this processing run
    run_path = output_folder_ + os.sep + f'{datetime.now().month}_{datetime.now().day}_{datetime.now().hour}_{datetime.now().minute}_{datetime.now().second}'

    # Read plate info
    plate_info = pd.read_excel(os.path.join(input_folder_, 'Plate_Info.xlsx'), usecols='A:M', sheet_name=None, index_col=0)

    # Write an excel file that can be read into jupyter notebook with minimal parsing.
    xlwriter_int = pd.ExcelWriter(os.path.join(run_path, 'intensities.xlsx'))

    if not os.path.isdir(run_path):
        os.mkdir(run_path)

    # ================
    # loop over images => good place for multiproc?  careful with columns in report
    # ================
    well_dirs = [d
                 for d in os.listdir(input_folder_)
                 if re.match(r'[A-P][0-9]{1,2}-Site_0', d)]
    # sort by letter, then by number (with '10' coming AFTER '9')
    well_dirs.sort(key=lambda x: (x[0], int(x[1:-7])))

    # well_dirs = ['A7-Site_0']

    well_path = None
    output_name = None
    int_well = []

    for well_dir in well_dirs:

        well_image_dir = [file for file in os.listdir(os.path.join(input_folder_, well_dir))
                          if '.png' in file or '.tif' in file or '.jpg' in file][0]
        image, image_name = read_to_grey(os.path.join(input_folder_, well_dir), well_image_dir)
        print(well_dir)

        # get well intensity
        well_mask = get_well_mask(image, segmethod='otsu')
        int_well_ = get_well_intensity(image, well_mask)
        int_well.append(int_well_)

        # SAVE FOR DEBUGGING
        if debug:
            well_path = os.path.join(run_path)
            os.makedirs(run_path, exist_ok=True)
            output_name = os.path.join(well_path, well_dir)

            # Save mask of the well, cropped grayscale image, cropped spot segmentation.
            io.imsave(output_name + "_well_mask.png",
                      (255 * well_mask).astype('uint8'))

            # Save masked image
            img_ = image.copy()
            img_[~well_mask] = 0
            io.imsave(output_name + "_masked_image.png",
                      (img_/256).astype('uint8'))

    # df_int = pd.DataFrame({'Well': pd.Series([well_dir.split(sep='-')[0] for well_dir in well_dirs]),
    #                        'Intensity': pd.Series(int_well)})
    df_int = pd.DataFrame(np.reshape(int_well, (8, 12)), index=list(string.ascii_uppercase[:8]), columns=range(1,13))
    plate_info.update({'intensity': df_int})
    for k, v in plate_info.items():
        v.to_excel(xlwriter_int, sheet_name=k)
    xlwriter_int.close()

    stop = time.time()
    print(f"\ttime to process={stop - start}")
