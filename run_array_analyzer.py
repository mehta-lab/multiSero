# bchhun, {2020-03-23}

import getopt
import sys

from array_analyzer.transform.property_filters import *
from array_analyzer.workflows import icp_wf, segmentation_wf


def main(argv):
    inputfolder = ''
    outputfolder = ''
    debug = False
    method = 'fit'
    try:
        options, remainder = getopt.getopt(argv, "hi:o:dm:",
                                           ["help","ifile=", "ofile=", "debug=", "method="])
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
        elif opt in ("-m", "--method"):
            method = arg
            assert method in ['fit', 'interp'], \
                ValueError('"method" has to be "fit" or "interp"')
        elif opt in ("-d", "--debug"):
            print('debug mode on, saving well and spot images')
            debug = True

    if not os.path.isdir(inputfolder):
        raise ValueError("input folder is not a folder or not supplied")

    if not os.path.isdir(outputfolder):
        os.makedirs(outputfolder)

    if method == 'fit':
        icp_wf.icp(inputfolder, outputfolder, debug)
    elif method == 'interp':
        segmentation_wf.seg(inputfolder, outputfolder, method='interp', debug=debug)
    else:
        raise KeyError(f"method {method} is not implemented")


if __name__ == "__main__":
    # Fluplate - old imaged on scienion
    # input_path = '/Volumes/GoogleDrive/My Drive/ELISAarrayReader/' \
    #              'images_scienion/Plates_given_to_manu/2020-01-15_plate4_AEP_Feb3_6mousesera'

    # Fluplate - old imaged on octopi
    # input_path = "/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_octopi/20200325AdamsPlate/Averaged/500us"

    # Fluplate - new imaged on scienion 3/30/2020
    input_path = '/Volumes/GoogleDrive/My Drive/ELISAarrayReader/' \
                 'images_scienion/2020-03-30-15-08-05-COVD_March25_fluplatetest_AdamsPlate'

    # local drives to write results
    # output_path = '/Users/shalin.mehta/Documents/images_local/2020-01-15_plate4_AEP_Feb3_6mousesera/'
    # output_path = '/Users/ivan.ivanov/Documents/images_local/' \
    #               'Plates_given_to_manu/2020-01-15_plate4_AEP_Feb3_6mousesera'
    output_path = '/Users/bryant.chhun/Desktop/Data/array-imager/Plates_given_to_manu/expt_merge_seg_icp/icp'

    method = 'fit'  # 'fit' or 'interp'
    flags = ['-i', input_path, '-o', output_path, '-d', '-m', method]

    main(flags)
    # main(sys.argv[1:])
