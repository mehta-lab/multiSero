import argparse
import os

import array_analyzer.workflows.registration_wf as registration_wf
import array_analyzer.workflows.interpolation_wf as interpolation_wf
import array_analyzer.workflows.well_wf as well_wf


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
        '-f', '--format',
        type=str,
        choices=['well', 'array'],
        default='array',
        help="Experiment format: well plate or printed arrays",
    )
    parser.add_argument(
        '-ad', '--array_detection',
        type=str,
        choices=['interp', 'fit'],
        default='interp',
        help="Method for detected spots: interpolation or registration",
    )
    parser.add_argument(
        '-d', '--debug',
        dest='debug',
        action='store_true',
        help="Write debug plots of well and spots",
    )
    parser.set_defaults(debug=False)

    return parser.parse_args()


def run_workflow(input_dir, output_dir, format, array_detection, debug=False):
    """
    For each image in input directory, run either interpolation (default)
    or registration of fiducials workflow.
    An xlsx file (and potentially debug plots) will be written to output directory.

    :param str input_dir: Input directory path
    :param str output_dir: Output directory path
    :param str format: array or well
    :param str array_detection: interp (interpolation) or fit (registration of fiducials)
    :param bool debug: Write debug plots to output directory
    """

    if not os.path.isdir(input_dir):
        raise ValueError("input directory is not a directory or doesn't exist")

    os.makedirs(output_dir, exist_ok=True)

    if format == 'array':
        if array_detection == 'fit':
            registration_wf.point_registration(
                input_dir,
                output_dir,
                debug=debug,
            )
        elif array_detection == 'interp':
            interpolation_wf.interp(
                input_dir,
                output_dir,
                method='interp',
                debug=debug,
            )
    elif format == 'well':
        well_wf.well_analysis(
            input_dir,
            output_dir,
            debug=debug,
        )


if __name__ == '__main__':
    args = parse_args()
    run_workflow(
        input_dir=args.input,
        output_dir=args.output,
        format=args.format,
        array_detection=args.array_detection,
        debug=args.debug,
    )
