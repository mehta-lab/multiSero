import os
from testfixtures import TempDirectory
import unittest
from unittest.mock import patch

import pysero as pysero


class TestPysero(unittest.TestCase):
    # TODO: Add tests for entire workflow

    def setUp(self):
        # TODO: Setup test images too
        self.tempdir = TempDirectory()
        self.temp_dir = self.tempdir.path
        self.input_dir = os.path.join(self.temp_dir, 'inputs')
        self.output_dir = os.path.join(self.temp_dir, 'outputs')
        self.tempdir.makedir(self.input_dir)
        self.tempdir.makedir(self.output_dir)

    def tearDown(self):
        """
        Tear down temporary folder and file structure
        """
        TempDirectory.cleanup_all()

    def test_parse_args(self):
        with patch('argparse._sys.argv',
                   ['python',
                    '-e',
                    '--input', self.input_dir,
                    '--output', self.output_dir,
                    '-wf', 'array_fit',
                    '--debug']):
            parsed_args = pysero.parse_args()
            self.assertEqual(parsed_args.input, self.input_dir)
            self.assertEqual(parsed_args.output, self.output_dir)
            self.assertTrue(parsed_args.debug)
            self.assertEqual(parsed_args.workflow, 'array_fit')

    def test_parse_args_invalid_method(self):
        with patch('argparse._sys.argv',
                   ['python',
                    '--extract_od',
                    '--input', self.input_dir,
                    '--output', self.output_dir,
                    '--workflow', 'magic']):
            with self.assertRaises(BaseException) as c:
                pysero.parse_args()

    def test_parse_args_mutially_exclusive(self):
        with patch('argparse._sys.argv',
                   ['python',
                    '--extract_od',
                    '--analyze_od',
                    '--input', self.input_dir,
                    '--output', self.output_dir]):
            with self.assertRaises(BaseException) as c:
                pysero.parse_args()
