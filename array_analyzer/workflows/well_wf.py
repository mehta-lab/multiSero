import sys, getopt, os

import array_analyzer.extract.image_parser as image_parser
import array_analyzer.extract.img_processing as processing
import array_analyzer.utils.io_utils as io_utils

import time
from datetime import datetime
import skimage.io as io
import pandas as pd
import string
import re
import numpy as np


def well_analysis(input_folder_, output_folder_, method='segmentation', debug=False):
    """
    Workflow that pulls all images scanned on a multi-well plate in a standard ELISA format (one antigen per well)
    It loops over the images in the input_folder (for images acquired using Micro-Manager ONLY).
        Extracts the center of the well, calculates the median intensity of that spot and background, then computes OD.
        Finally, it writes a summary report (.xlsx) containing plate info and the computed values.

    :param input_folder_: str path to experiment directory
    :param output_folder_: str output path to write report and diagnostic images
    :param method: str 'segmentation' or 'crop'.  Methods to estimate the boundaries of the well
    :param debug: bool.  Whether to write debug images.
    :return:
    """

    start = time.time()

    # save a sub path for this processing run
    run_path = output_folder_ + os.sep + f'{datetime.now().month}_{datetime.now().day}_{datetime.now().hour}_{datetime.now().minute}_{datetime.now().second}'

    # Read plate info
    plate_info = pd.read_excel(os.path.join(input_folder_, 'Plate_Info.xlsx'), usecols='A:M', sheet_name=None, index_col=0)

    # Write an excel file that can be read into jupyter notebook with minimal parsing.
    xlwriter_int = pd.ExcelWriter(os.path.join(run_path, 'intensities.xlsx'))

    if not os.path.isdir(run_path):
        os.mkdir(run_path)

    # get well directories
    well_dirs = [d
                 for d in os.listdir(input_folder_)
                 if re.match(r'[A-P][0-9]{1,2}-Site_0', d)]
    # sort by letter, then by number (with '10' coming AFTER '9')
    well_dirs.sort(key=lambda x: (x[0], int(x[1:-7])))

    # well_dirs = ['A7-Site_0']

    int_well = []
    for well_dir in well_dirs:

        # read image
        well_image_file = [file for file in os.listdir(os.path.join(input_folder_, well_dir))
                          if '.png' in file or '.tif' in file or '.jpg' in file][0]
        image, image_name = io_utils.read_to_grey(
            os.path.join(input_folder_, well_dir),
            well_image_file,
        )
        print(well_dir)

        # measure intensity
        if method == 'segmentation':
            # segment well using otsu thresholding
            well_mask = image_parser.get_well_mask(image, segmethod='otsu')
            int_well_ = image_parser.get_well_intensity(image, well_mask)

        elif method == 'crop':
            # get intensity at square crop in the middle of the image
            img_size = image.shape
            radius = np.floor(0.1 * np.min(img_size)).astype('int')
            cx = np.floor(img_size[1]/2).astype('int')
            cy = np.floor(img_size[0]/2).astype('int')
            im_crop = processing.crop_image(image, cx, cy, radius, border_=0)

            well_mask = np.ones_like(im_crop, dtype='bool')
            int_well_ = image_parser.get_well_intensity(im_crop, well_mask)

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
            if method == 'segmentation':
                img_ = image.copy()
                img_[~well_mask] = 0
            elif method == 'crop':
                img_ = im_crop.copy()
                img_[~well_mask] = 0
            else:
                raise NotImplementedError(f'method of type {method} not supported')
            io.imsave(output_name + "_masked_image.png",
                      (img_/256).astype('uint8'))

    df_int = pd.DataFrame(np.reshape(int_well, (8, 12)), index=list(string.ascii_uppercase[:8]), columns=range(1,13))
    plate_info.update({'intensity': df_int})

    # compute optical density
    sample_info = plate_info['sample']
    blanks = np.any(np.dstack((sample_info == 'Blank', sample_info == 'blank', sample_info == 'BLANK')), axis=2)
    if blanks.any():
        # blank intensity is averaged over all blank wells
        int_blank = np.mean(df_int.to_numpy()[blanks])
        df_od = np.log10(int_blank / df_int)
        plate_info.update({'od': df_od})

    # save analysis results
    for k, v in plate_info.items():
        v.to_excel(xlwriter_int, sheet_name=k)
    xlwriter_int.close()

    stop = time.time()
    print(f"\ttime to process={stop - start}")
