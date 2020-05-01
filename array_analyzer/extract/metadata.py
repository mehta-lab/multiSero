import os
import numpy as np
from datetime import datetime
import pandas as pd

import array_analyzer.extract.txt_parser as txt_parser
import array_analyzer.extract.constants as c


class MetaData:

    def __init__(self, input_folder_, output_folder_):
        """
        Parses metadata spreadsheets then populates all necessary ARRAY data structures
        Extracts all necessary constants and assigns them in the constants.py namespace

        :param input_folder_: str full path to metadata spreadsheet
        :param output_folder_: str full path to output folder for reports and diagnostics
        """

        self.fiduc, self.spots, self.repl, self.params = None, None, None, None

        # parse fiducials, spot types, antigens, and hardware parameters from metadata
        if c.METADATA_EXTENSION == 'xml':
            # check that exactly one .xml exists
            xml = [f for f in os.listdir(input_folder_) if '.xml' in f]
            if len(xml) > 1:
                raise IOError("more than one .xml file found, aborting")
            xml_path = os.path.join(input_folder_, xml[0])

            # parsing .xml
            self.fiduc, self.spots, self.repl, self.params = txt_parser.create_xml_dict(xml_path)

        elif c.METADATA_EXTENSION == 'well':
            self._set_run_path(output_folder_)
            return

        elif c.METADATA_EXTENSION == 'csv':
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

        elif c.METADATA_EXTENSION == 'xlsx':
            # check that exactly one .xlsx exists
            # xlsxs = [f for f in os.listdir(input_folder_) if '.xlsx' in f]
            # if len(xlsxs) > 1:
            #     raise IOError("more than one .xlsx file found, aborting")
            # xlsx_path = os.path.join(input_folder_, xlsxs[0])

            # check that properly named .xlsx exists
            if not os.path.isfile(os.path.join(input_folder_, 'Metadata_and_Plate_configuration.xlsx')):
            # if not os.path.isfile(os.path.join(input_folder_, 'pysero_output_data_metadata.xlsx')):
                raise IOError("required metadata file named 'Metadata_and_Plate_configuration.xlsx' does not exist")
            xlsx_path = os.path.join(input_folder_, 'Metadata_and_Plate_configuration.xlsx')

            # check that the xlsx file contains necessary worksheets
            sheets = pd.read_excel(xlsx_path, sheet_name=None)
            if 'imaging_and_array_parameters' not in sheets.keys():
                raise IOError("sheet by name 'imaging_and_array_parameters' not present in excel file, aborting")
            if 'array_antigens' not in sheets.keys():
                raise IOError("sheet by name 'array_antigens' not present in excel file, aborting")

            # parsing .xlsx
            self.fiduc, _, self.repl, self.params = txt_parser.create_xlsx_dict(xlsx_path)

            # parsing .xlsx using pandas
            # self.fiduc, _, self.repl, self.params = txt_parser.create_xlsx_array(xlsx_path)
            # c.FIDUCIAL_ARRAY = self.fiduc
            # c.ANTIGEN_ARRAY = self.repl

        else:
            raise NotImplementedError(f"metadata with extension {c.METADATA_EXTENSION} is not supported")

        # setting constants
        # if fiduc_:
        #     c.fiducials = fiduc_
        # if spots_:
        #     c.spots = spots_
        # if repl_:
        #     c.replicates = repl_
        # if params_:
        c.params['rows'] = int(self.params['rows'])
        c.params['columns'] = int(self.params['columns'])
        c.params['v_pitch'] = float(self.params['v_pitch'])
        c.params['h_pitch'] = float(self.params['h_pitch'])
        c.params['spot_width'] = float(self.params['spot_width'])
        if c.METADATA_EXTENSION == 'xml':
            c.params['pixel_size'] = c.params['pixel_size_scienion']
            # these params are present in .xml but not used
            # c.params['bg_offset'] = int(params_['bg_offset'])
            # c.params['bg_thickness'] = int(params_['bg_thickness'])
            # c.params['max_diam'] = int(params_['max_diam'])
            # c.params['min_diam'] = int(params_['min_diam'])
        else:
            c.params['pixel_size'] = float(self.params['pixel_size'])

        # setting constant arrays
        if c.METADATA_EXTENSION == 'xml':
            self._create_spot_id_array()
            self._create_spot_type_array()
        self._create_fiducials_array()
        self._create_antigen_array()

        # setting location of fiducials and other useful parameters
        self._calculate_fiduc_coords()
        self._calculate_fiduc_idx()
        self._calc_spot_dist()
        self._set_run_path(output_folder_)

        # setting 96-well constants
        self._calc_image_to_well()
        self._calc_empty_plate_const()

    def _create_spot_id_array(self):
        """
        Creates an empty ndarray of strings, whose rows and columns match the printed array's rows/cols
            Sets the ARRAY constant corresponding to "Spot-ID" based on metadata
        *** note: this is used ONLY for .xml metadata files generated by the sciReader ***
        :return:
        """
        self.spot_ids = np.empty(shape=(c.params['rows'], c.params['columns']), dtype='U100')
        c.SPOT_ID_ARRAY = txt_parser.populate_array_id(self.spot_ids, self.spots)

    def _assign_params(self, params_):
        c.params['rows'] = int(params_['rows'])
        c.params['columns'] = int(params_['columns'])
        c.params['v_pitch'] = float(params_['v_pitch'])
        c.params['h_pitch'] = float(params_['h_pitch'])
        c.params['spot_width'] = float(params_['spot_width'])
        if c.METADATA_EXTENSION == 'xml':
            c.params['pixel_size'] = c.params['pixel_size_scienion']
            # these params are present in .xml but not used
            # c.params['bg_offset'] = int(params_['bg_offset'])
            # c.params['bg_thickness'] = int(params_['bg_thickness'])
            # c.params['max_diam'] = int(params_['max_diam'])
            # c.params['min_diam'] = int(params_['min_diam'])
        else:
            c.params['pixel_size'] = float(params_['pixel_size'])

    def _create_spot_type_array(self):
        """
        Creates an empty ndarray of strings, whose rows and columns match the printed array's rows/cols
            Sets the ARRAY constant corresponding to "Spot Type" based on metadata
            "Spot Type" are values like "Positive Control", "Diagnostic"
        :return:
        """
        self.spot_type = np.empty(shape=(c.params['rows'], c.params['columns']), dtype='U100')
        c.SPOT_TYPE_ARRAY = txt_parser.populate_array_spots_type(self.spot_type, self.spots, self.fiduc)

    def _create_fiducials_array(self):
        """
        Creates an empty ndarray of strings, whose rows and columns match the printed array's rows/cols
            Sets the ARRAY constant corresponding to "Fiducial" based on metadata
            "Fiducial" are values like "Fiducial" or "Reference, Diagnostic" and are used to align array spots
        :return:
        """
        self.fiducials_array = np.empty(shape=(c.params['rows'], c.params['columns']), dtype='U100')
        c.FIDUCIAL_ARRAY = txt_parser.populate_array_fiduc(self.fiducials_array, self.fiduc)

    def _create_antigen_array(self):
        """
        Creates an empty ndarray of strings, whose rows and columns match the printed array's rows/cols
            Assigns the corresponding "Antigen" based on metadata
            "Antigens" are descriptive values of the antigen at each spot location.

        This is the only array creator that requires one of the above arrays (SPOT_ID_ARRAY, .xml meta ONLY)
            multiple "spot-ids" can contain the same "antigen".  Thus it's required to pass the "SPOT_ID_ARRAY"
        :return:
        """
        self.antigen_array = np.empty(shape=(c.params['rows'], c.params['columns']), dtype='U100')
        if c.METADATA_EXTENSION == 'xml':
            if c.SPOT_ID_ARRAY.size == 0:
                raise AttributeError("attempting to create antigen array before SPOT_ID_ARRAY is assigned")
            c.ANTIGEN_ARRAY = txt_parser.populate_array_antigen_xml(self.antigen_array, c.SPOT_ID_ARRAY, self.repl)
        elif c.METADATA_EXTENSION == 'csv' or c.METADATA_EXTENSION == 'xlsx':
            c.ANTIGEN_ARRAY = txt_parser.populate_array_antigen(self.antigen_array, self.repl)

    def _calculate_fiduc_coords(self):
        """
        Calculate and set the fiducial coordinates like:
            FIDUCIALS = [(0, 0), (0, 1), (0, 5), (7, 0), (7, 5)]
            fiducial coordinates are labeled as "Reference, Diagnostic"
        :return:
        """
        x, y = np.where(c.FIDUCIAL_ARRAY == 'Reference, Diagnostic')
        if x.size == 0 or y.size == 0:
            x, y = np.where(c.FIDUCIAL_ARRAY == 'Fiducial')
        c.FIDUCIALS = list(zip(x, y))

    def _calculate_fiduc_idx(self):
        """
        Calculate fiducial index like
            FIDUCIALS_IDX = [0, 5, 6, 30, 35]\


            FIDUCIALS_IDX = [0, 7, 8, 40, 47] for 8 columns
        :return:
        """
        c.FIDUCIALS_IDX = list(np.where(c.FIDUCIAL_ARRAY.flatten() == 'Reference, Diagnostic')[0])
        if len(c.FIDUCIALS_IDX) == 0:
            c.FIDUCIALS_IDX = list(np.where(c.FIDUCIAL_ARRAY.flatten() == 'Fiducial')[0])

    def _calc_spot_dist(self):
        """
        Calculate distance between spots in both pixels and microns
        :return:
        """
        v_pitch_mm = c.params['v_pitch']
        h_pitch_mm = c.params['h_pitch']
        pix_size = c.params['pixel_size']

        # assuming similar v_pitch and h_pitch, average to get the SPOT_DIST
        v_pitch_pix = v_pitch_mm/pix_size
        h_pitch_pix = h_pitch_mm/pix_size
        c.SPOT_DIST_PIX = np.mean([v_pitch_pix, h_pitch_pix]).astype('uint8')

        # convert the SPOT_DIST to microns, 0 - 255
        c.SPOT_DIST_UM = np.mean([v_pitch_mm * 1000, h_pitch_mm * 1000]).astype('uint8')

    # set filesaving run_path
    def _set_run_path(self, output_folder):
        """
        Create the output folder for this analysis run
        folder is unique to the second, and can contain both reports and diagnostics
        :param output_folder: str path to output folder specified at CLI
        :return:
        """
        c.RUN_PATH = os.path.join(
            output_folder,
            '_'.join([str(datetime.now().month),
                      str(datetime.now().day),
                      str(datetime.now().hour),
                      str(datetime.now().minute),
                      str(datetime.now().second)]),
        )
        if not os.path.isdir(c.RUN_PATH):
            os.mkdir(c.RUN_PATH)

    # create image-to-well mapping dictionary
    def _calc_image_to_well(self):
        """
        Calculate the mapping from ImageName: (row, col) position in the plate.

        :return:
        """
        # assuming file names are "rowcol" or "A1" - "H12"
        rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        cols = list(range(1, 13))
        for r_idx, row in enumerate(rows):
            for col in cols:
                c.IMAGE_TO_WELL[row+str(col)] = (r_idx+1, col)

    def _calc_empty_plate_const(self):
        """
        initialize report arrays assuming a 96-well plate format
        :return:
        """
        c.WELL_BG_ARRAY = np.empty((8, 12), dtype=object)
        c.WELL_INT_ARRAY = np.empty((8, 12), dtype=object)
        c.WELL_OD_ARRAY = np.empty((8, 12), dtype=object)
