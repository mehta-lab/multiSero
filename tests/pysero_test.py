import pytest
from unittest.mock import patch

import pysero as pysero


def test_parse_args():
    with patch('argparse._sys.argv',
               ['python',
                '-e',
                '--input', 'input_dir_name',
                '--output', 'output_dir_name',
                '-wf', 'array_fit',
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
