import argparse
import os

import array_analyzer.workflows.registration_workflow as registration_wf
import array_analyzer.workflows.interpolation_wf as interpolation_wf
import array_analyzer.workflows.well_wf as well_wf
import array_analyzer.extract.constants as constants


def parse_args():
    """
    Parse command line arguments for CLI.

    :return: namespace containing the arguments passed.
    """
    parser = argparse.ArgumentParser()

    # make sure that only extract_od or analyze_od stages are passed.
    stage = parser.add_mutually_exclusive_group(required=True)
    stage.add_argument(
        '-e', '--extract_od',
        action= 'store_const',
        const=True,
    )
    stage.add_argument(
        '-a', '--analyze_od',
        action='store_const',
        const=True,
    )
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help="Input directory path",
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help="Output directory path",
    )
    parser.add_argument(
        '-wf', '--workflow',
        type=str,
        choices=['well_segmentation', 'well_crop', 'array_interp', 'array_fit'],
        default='array_interp',
        help="Workflow to automatically identify and extract intensities from experiment.  "
             "'Well' experiments are for standard ELISA.  "
             "'Array' experiments are for ELISA assays using antigen arrays printed with Scienion Array Printer ",
    )
    parser.add_argument(
        '-d', '--debug',
        dest='debug',
        action='store_true',
        help="Write debug plots of well and spots",
    )
    parser.set_defaults(debug=False)
    parser.add_argument(
        '-m', '--metadata',
        type=str,
        choices=['xml', 'xlsx', 'well'],
        help="specify the file extension for the experiment metadata"
    )
    parser.set_defaults(metadata='xlsx')

    return parser.parse_args()


def extract_od(input_dir, output_dir, workflow):
    """
    For each image in input directory, run either interpolation (default)
    or registration of fiducials workflow.
    An xlsx file (and potentially debug plots) will be written to output directory.

    :param str input_dir: Input directory path
    :param str output_dir: Output directory path
    :param str workflow: str, one of 'array_interp', 'array_fit', 'well_segmentation', 'well_crop'
        <plate>_<method> format:
            <plate> describes the printing style of the antigen (array or ELISA)
            <method> describes the spot segmentation and extraction approach
    """

    if not os.path.isdir(input_dir):
        raise ValueError("input directory is not a directory or doesn't exist")

    os.makedirs(output_dir, exist_ok=True)

    if workflow == 'array_interp':
        interpolation_wf.interp(
            input_dir,
            output_dir,
        )
    elif workflow == 'array_fit':
        registration_wf.point_registration(
            input_dir,
            output_dir,
        )
    elif workflow == 'well_segmentation':
        well_wf.well_analysis(
            input_dir,
            output_dir,
            method='segmentation',
        )
    elif workflow == 'well_crop':
        well_wf.well_analysis(
            input_dir,
            output_dir,
            method='crop',
        )


if __name__ == '__main__':
    args = parse_args()

    constants.METADATA_EXTENSION = args.metadata
    constants.DEBUG = args.debug

    if args.extract_od:
        extract_od(
            input_dir=args.input,
            output_dir=args.output,
            workflow=args.workflow,
        )
    elif args.analyze_od:
        raise NotImplementedError(
            'Automated interpretation is coming. '
            'See the interpretation folder for examples of interpretation scripts.',
        )
