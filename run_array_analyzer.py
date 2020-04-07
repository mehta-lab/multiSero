# bchhun, {2020-03-23}

import getopt
import sys

from array_analyzer.transform.property_filters import *
from array_analyzer.workflows import icp_wf, interpolation_wf

FIDUCIALS = [(0, 0), (0, 1), (0, 5), (7, 0), (7, 5)]
SCENION_SPOT_DIST = 82


def main(argv):
    inputfolder = ''
    outputfolder = ''
    debug = False
    method = 'interp'
    try:
        options, remainder = getopt.getopt(argv, "hi:o:dm:",
                                           ["help", "ifile=", "ofile=", "debug=", "method="])
    except getopt.GetoptError:
        print('run_array_analyzer.py -i <inputfolder> -o <outputfolder> -m <method>')
        sys.exit(2)

    for opt, arg in options:
        if opt == '-h':
            print('run_array_analyzer.py -i <inputfolder> -o <outputfolder> -m <method>')
            print('\t inputfolder: path to folder containing .xml and images')
            print('\t outputfolder: path to folder for writing report and troubleshooting images')
            print('\t method: one of "fit" or "interp"')
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
        icp_wf.point_registration(inputfolder, outputfolder, debug)
    elif method == 'interp':
        interpolation_wf.interp(inputfolder, outputfolder, method='interp', debug=debug)
    else:
        raise KeyError(f"method {method} is not implemented")


if __name__ == "__main__":
    # Fluplate - old imaged on scienion
    # input_path = '/Volumes/GoogleDrive/My Drive/ELISAarrayReader/' \
    #              'images_scienion/Plates_given_to_manu/2020-01-15_plate4_AEP_Feb3_6mousesera'

    # Fluplate - old imaged on octopi
    # input_path = "/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_octopi/20200325AdamsPlate/Averaged/500us"

    # Fluplate - new imaged on scienion 3/30/2020
    # input_path = f'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-03-30-15-08-05-COVD_March25_fluplatetest_AdamsPlate'

    # Fluplate - blinded test imaged on scienion 4/4/2020
    input_path = f'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-04-04-14-18-32-COVID_April4_flusecondplate'

    datadir=os.path.basename(os.path.normpath(input_path))
    # Write results to local folder, copy them to drive if useful.
    output_path = os.path.join(os.path.expanduser('~/Documents/images_local/'),datadir)

    method = 'interp'  # 'fit' or 'interp'
    flags = ['-i', input_path, '-o', output_path, '-d', '-m', method]

    # main(flags)
    main(sys.argv[1:])
