import os
import pytest
from unittest.mock import patch

import pysero as pysero
import array_analyzer.extract.constants as constants


def test_parse_args():
    with patch('argparse._sys.argv',
               ['python',
                '-e',
                '--input', 'input_dir_name',
                '--output', 'output_dir_name',
                '-wf', 'array_fit',
                '-m', 'xml',
                '--debug']):
        parsed_args = pysero.parse_args()
        assert parsed_args.input == 'input_dir_name'
        assert parsed_args.output == 'output_dir_name'
        assert parsed_args.debug is True
        assert parsed_args.workflow == 'array_fit'


def test_parse_args_invalid_method():
    with patch('argparse._sys.argv',
               ['python',
                '--extract_od',
                '--input', 'input_dir_name',
                '--output', 'output_dir_name',
                '--workflow', 'magic']):
        with pytest.raises(BaseException):
            pysero.parse_args()


def test_parse_args_mutially_exclusive():
    with patch('argparse._sys.argv',
               ['python',
                '--extract_od',
                '--analyze_od',
                '--input', 'input_dir_name',
                '--output', 'output_dir_name']):
        with pytest.raises(BaseException):
            pysero.parse_args()


@pytest.mark.slow
def test_extract_od_scienion(scienion_dir):

    constants.METADATA_EXTENSION = 'xlsx'
    constants.DEBUG = True

    # Run array fit workflow in debug mode
    pysero.extract_od(
        input_dir=scienion_dir[0],
        output_dir=scienion_dir[1],
        workflow='array_fit',
    )

    # Get run dir
    run_dir = os.path.join(scienion_dir[1], os.listdir(scienion_dir[1])[0])
    output_list = os.listdir(run_dir)
    # Check that debug images are present
    well_names = ['A1', 'A2', 'B11', 'B12']
    debug_names = ['_od.png',
                   '_registration.png' ,
                   '_composite_spots_img.png',
                   '_crop_bg_overlay.png']
    for well_name in well_names:
        for debug_name in debug_names:
            file_name = well_name + debug_name
            assert file_name in output_list
    # Check that xlsx files are present
    assert 'median_backgrounds.xlsx' in output_list
    assert 'median_intensities.xlsx' in output_list
    assert 'median_ODs.xlsx' in output_list
    # TODO: Open xlsx and check contents
