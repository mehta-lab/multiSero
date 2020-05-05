import os
import numpy as np
from datetime import datetime
import pandas as pd
import shutil

import array_analyzer.extract.txt_parser as txt_parser
import array_analyzer.extract.constants as constants


class MetaData:

    def __init__(self, input_folder_, output_folder_):
        """
        Parses metadata spreadsheets then populates all necessary ARRAY data structures
        Extracts all necessary constants and assigns them in the constants.py namespace

        :param input_folder_: str full path to metadata spreadsheet
        :param output_folder_: str full path to output folder for reports and diagnostics
        """

        self.fiduc, self.spots, self.repl, self.params = None, None, None, None
        self.xlsx_path = None
        self.xml_path = None
        constants.INPUT_FOLDER = input_folder_
        constants.OUTPUT_FOLDER = output_folder_

        # parse fiducials, spot types, antigens, and hardware parameters from metadata
        if constants.METADATA_EXTENSION == 'xml':
            # check that exactly one .xml exists
            xml = [f for f in os.listdir(input_folder_) if '.xml' in f]
            if len(xml) > 1:
                raise IOError("more than one .xml file found, aborting")
            self.xml_path = os.path.join(input_folder_, xml[0])

            # parsing .xml
            self.fiduc, self.spots, self.repl, self.params = txt_parser.create_xml_dict(self.xml_path)

        elif constants.METADATA_EXTENSION == 'well':
            self._set_run_path(output_folder_)
            return

        elif constants.METADATA_EXTENSION == 'csv':
            # check that three .csvs exist
            three_csvs = ['array_format_antigen', 'array_format_type', 'array_parameters']
            csvs = [f for f in os.listdir(input_folder_) if '.csv' in f]
            if len(csvs) != 3:
                raise IOError("incorrect number of .csv files found, aborting")
            for target in three_csvs:
                if True not in [target in file for file in csvs]:
                    raise IOError(f".csv file with substring {target} is missing")

            csv_paths = [os.path.join(input_folder_, one_csv) for one_csv in csvs]

            # parsing .csv
            self.fiduc, _, self.repl, self.params = txt_parser.create_csv_dict(csv_paths)

        elif constants.METADATA_EXTENSION == 'xlsx':
            # check that exactly one .xlsx exists
            # xlsxs = [f for f in os.listdir(input_folder_) if '.xlsx' in f]
            # if len(xlsxs) > 1:
            #     raise IOError("more than one .xlsx file found, aborting")
            # xlsx_path = os.path.join(input_folder_, xlsxs[0])

            # check that properly named .xlsx exists
            if not os.path.isfile(os.path.join(input_folder_, 'pysero_output_data_metadata.xlsx')):
                raise IOError("required metadata file named 'pysero_output_data_metadata.xlsx' does not exist")
            self.xlsx_path = os.path.join(input_folder_, 'pysero_output_data_metadata.xlsx')

            # check that the xlsx file contains necessary worksheets
            sheets = pd.read_excel(self.xlsx_path, sheet_name=None)
            if 'imaging_and_array_parameters' not in sheets.keys():
                raise IOError("sheet by name 'imaging_and_array_parameters' not present in excel file, aborting")
            if 'antigen_array' not in sheets.keys():
                raise IOError("sheet by name 'array_antigens' not present in excel file, aborting")

            # parsing .xlsx
            self.fiduc, _, self.repl, self.params = txt_parser.create_xlsx_dict(self.xlsx_path)

            # parsing .xlsx using pandas !! not tested or finished yet
            # self.fiduc, _, self.repl, self.params = txt_parser.create_xlsx_array(xlsx_path)
            # c.FIDUCIAL_ARRAY = self.fiduc
            # c.ANTIGEN_ARRAY = self.repl

        else:
            raise NotImplementedError(f"metadata with extension {constants.METADATA_EXTENSION} is not supported")

        # set hardware and array parameters
        self._assign_params(self.params)

        # setting constant arrays
        if constants.METADATA_EXTENSION == 'xml':
            self._create_spot_id_array()
            self._create_spot_type_array()
        self._create_fiducials_array()
        self._create_antigen_array()

        # setting location of fiducials and other useful parameters
        self._calculate_fiduc_coords()
        self._calculate_fiduc_idx()
        self._calc_spot_dist()
        self._set_run_path(output_folder_)
        self._copy_metadata_to_output()

        # setting 96-well constants
        self._calc_image_to_well()
        self._calc_empty_plate_const()

    @staticmethod
    def _assign_params(params_):
        constants.params['rows'] = int(params_['rows'])
        constants.params['columns'] = int(params_['columns'])
        constants.params['v_pitch'] = float(params_['v_pitch'])
        constants.params['h_pitch'] = float(params_['h_pitch'])
        constants.params['spot_width'] = float(params_['spot_width'])
        if constants.METADATA_EXTENSION == 'xml':
            constants.params['pixel_size'] = constants.params['pixel_size_scienion']
        else:
            constants.params['pixel_size'] = float(params_['pixel_size'])

    def _create_spot_id_array(self):
        """
        Creates an empty ndarray of strings, whose rows and columns match the printed array's rows/cols
            Sets the ARRAY constant corresponding to "Spot-ID" based on metadata
        *** note: this is used ONLY for .xml metadata files generated by the sciReader ***
        :return:
        """
        self.spot_ids = np.empty(
            shape=(constants.params['rows'], constants.params['columns']),
            dtype='U100',
        )
        constants.SPOT_ID_ARRAY = txt_parser.populate_array_id(
            self.spot_ids,
            self.spots,
        )

    def _create_spot_type_array(self):
        """
        Creates an empty ndarray of strings, whose rows and columns match the printed array's rows/cols
            Sets the ARRAY constant corresponding to "Spot Type" based on metadata
            "Spot Type" are values like "Positive Control", "Diagnostic"
        :return:
        """
        self.spot_type = np.empty(
            shape=(constants.params['rows'], constants.params['columns']),
            dtype='U100',
        )
        constants.SPOT_TYPE_ARRAY = txt_parser.populate_array_spots_type(
            self.spot_type,
            self.spots,
            self.fiduc,
        )

    def _create_fiducials_array(self):
        """
        Creates an empty ndarray of strings, whose rows and columns match the printed array's rows/cols
            Sets the ARRAY constant corresponding to "Fiducial" based on metadata
            "Fiducial" are values like "Fiducial" or "Reference, Diagnostic" and are used to align array spots
        :return:
        """
        self.fiducials_array = np.empty(
            shape=(constants.params['rows'], constants.params['columns']),
            dtype='U100',
        )
        constants.FIDUCIAL_ARRAY = txt_parser.populate_array_fiduc(
            self.fiducials_array,
            self.fiduc,
        )

    def _create_antigen_array(self):
        """
        Creates an empty ndarray of strings, whose rows and columns match the printed array's rows/cols
            Assigns the corresponding "Antigen" based on metadata
            "Antigens" are descriptive values of the antigen at each spot location.

        This is the only array creator that requires one of the above arrays (SPOT_ID_ARRAY, .xml meta ONLY)
            multiple "spot-ids" can contain the same "antigen".  Thus it's required to pass the "SPOT_ID_ARRAY"
        :return:
        """
        self.antigen_array = np.empty(
            shape=(constants.params['rows'], constants.params['columns']),
            dtype='U100',
        )
        if constants.METADATA_EXTENSION == 'xml':
            if constants.SPOT_ID_ARRAY.size == 0:
                raise AttributeError("attempting to create antigen array "
                                     "before SPOT_ID_ARRAY is assigned")
            constants.ANTIGEN_ARRAY = txt_parser.populate_array_antigen_xml(
                self.antigen_array,
                constants.SPOT_ID_ARRAY,
                self.repl,
            )
        elif constants.METADATA_EXTENSION == 'csv' or constants.METADATA_EXTENSION == 'xlsx':
            constants.ANTIGEN_ARRAY = txt_parser.populate_array_antigen(
                self.antigen_array,
                self.repl,
            )

    @staticmethod
    def _calculate_fiduc_coords():
        """
        Calculate and set the fiducial coordinates like:
            FIDUCIALS = [(0, 0), (0, 1), (0, 5), (7, 0), (7, 5)]
            fiducial coordinates are labeled as "Reference, Diagnostic"
        :return:
        """
        x, y = np.where(constants.FIDUCIAL_ARRAY == 'Reference, Diagnostic')
        if x.size == 0 or y.size == 0:
            x, y = np.where(constants.FIDUCIAL_ARRAY == 'Fiducial')
        constants.FIDUCIALS = list(zip(x, y))

    @staticmethod
    def _calculate_fiduc_idx():
        """
        Calculate fiducial index like
            FIDUCIALS_IDX = [0, 5, 6, 30, 35]\


            FIDUCIALS_IDX = [0, 7, 8, 40, 47] for 8 columns
        :return:
        """
        constants.FIDUCIALS_IDX = list(np.where(
            constants.FIDUCIAL_ARRAY.flatten() == 'Reference, Diagnostic')[0])
        if len(constants.FIDUCIALS_IDX) == 0:
            constants.FIDUCIALS_IDX = list(np.where(
                constants.FIDUCIAL_ARRAY.flatten() == 'Fiducial')[0])

    @staticmethod
    def _calc_spot_dist():
        """
        Calculate distance between spots in both pixels and microns
        :return:
        """
        v_pitch_mm = constants.params['v_pitch']
        h_pitch_mm = constants.params['h_pitch']
        pix_size = constants.params['pixel_size']

        # assuming similar v_pitch and h_pitch, average to get the SPOT_DIST
        v_pitch_pix = v_pitch_mm/pix_size
        h_pitch_pix = h_pitch_mm/pix_size
        constants.SPOT_DIST_PIX = np.mean([v_pitch_pix, h_pitch_pix]).astype('uint8')

        # convert the SPOT_DIST to microns, 0 - 255
        constants.SPOT_DIST_UM = np.mean([v_pitch_mm * 1000, h_pitch_mm * 1000]).astype('uint8')

    # set filesaving run_path
    def _set_run_path(self, output_folder):
        """
        Create the output folder for this analysis run
        folder is unique to the second, and can contain both reports and diagnostics
        :param output_folder: str path to output folder specified at CLI
        :return:
        """
        constants.RUN_PATH = os.path.join(
            output_folder,
            ''.join(['pysero_',
                     os.path.basename(os.path.normpath(constants.INPUT_FOLDER)),
                     '_',
                     f"{datetime.now().year:04d}",
                     f"{datetime.now().month:02d}",
                     f"{datetime.now().day:02d}",
                     '_',
                     f"{datetime.now().hour:02d}",
                     f"{datetime.now().minute:02d}"]
                    )
        )
        os.makedirs(constants.RUN_PATH, exist_ok=True)

    def _copy_metadata_to_output(self):
        if constants.METADATA_EXTENSION == 'xlsx':
            shutil.copy2(self.xlsx_path, constants.RUN_PATH)
        elif constants.METADATA_EXTENSION == 'xml':
            shutil.copy2(self.xml_path, constants.RUN_PATH)

    # create image-to-well mapping dictionary
    @staticmethod
    def _calc_image_to_well():
        """
        Calculate the mapping from ImageName: (row, col) position in the plate.
        :return:
        """
        # assuming file names are "rowcol" or "A1" - "H12"
        rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        cols = list(range(1, 13))
        for r_idx, row in enumerate(rows):
            for col in cols:
                constants.IMAGE_TO_WELL[row+str(col)] = (r_idx+1, col)

    @staticmethod
    def _calc_empty_plate_const():
        """
        initialize report arrays assuming a 96-well plate format
        each element of these arrays contains a sub-array of antigens
        :return:
        """
        constants.WELL_BG_ARRAY = np.empty((8, 12), dtype=object)
        constants.WELL_INT_ARRAY = np.empty((8, 12), dtype=object)
        constants.WELL_OD_ARRAY = np.empty((8, 12), dtype=object)
