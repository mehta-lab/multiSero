import argparse
import os

import array_analyzer.workflows.registration_workflow as registration_wf
import array_analyzer.workflows.interpolation_wf as interpolation_wf
import array_analyzer.workflows.well_wf as well_wf
import array_analyzer.extract.constants as c


def parse_args():
    """
    Parse command line arguments for CLI.

    :return: namespace containing the arguments passed.
    """
    parser = argparse.ArgumentParser()
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
        help="Analysis workflow",
    )
    parser.add_argument(
        '-d', '--debug',
        dest='debug',
        action='store_true',
        help="Write debug plots of well and spots",
    )
    parser.add_argument(
        '-m', '--metadata',
        type=str,
        required=True,
        choices=['xml', 'csv', 'xlsx'],
        help="specify the file extension for the experiment metadata"
    )
    parser.set_defaults(debug=False)

    return parser.parse_args()


def run_workflow(input_dir, output_dir, workflow, debug=False):
    """
    For each image in input directory, run either interpolation (default)
    or registration of fiducials workflow.
    An xlsx file (and potentially debug plots) will be written to output directory.

    :param str input_dir: Input directory path
    :param str output_dir: Output directory path
    :param str workflow: 'well_segmentation', 'well_crop', 'array_interp', 'array_fit'
    :param bool debug: Write debug plots to output directory
    """

    if not os.path.isdir(input_dir):
        raise ValueError("input directory is not a directory or doesn't exist")

    os.makedirs(output_dir, exist_ok=True)

    if workflow == 'array_interp':
        interpolation_wf.interp(
            input_dir,
            output_dir,
            method='interp',
            debug=debug,
        )
    elif workflow == 'array_fit':
        registration_wf.point_registration(
            input_dir,
            output_dir,
            debug=debug,
        )
    elif workflow == 'well_segmentation':
        well_wf.well_analysis(
            input_dir,
            output_dir,
            method='segmentation',
            debug=debug,
        )
    elif workflow == 'well_crop':
        well_wf.well_analysis(
            input_dir,
            output_dir,
            method='crop',
            debug=debug,
        )


if __name__ == '__main__':
    args = parse_args()
    c.METADATA_EXTENSION = args.metadata

    run_workflow(
        input_dir=args.input,
        output_dir=args.output,
        workflow=args.workflow,
        debug=args.debug,
    )
