import os
import numpy as np
from datetime import datetime

import array_analyzer.extract.txt_parser as txt_parser
import array_analyzer.extract.constants as c


class MetaData:

    def __init__(self, input_folder_, output_folder_):
        """

        :param input_folder_:
        :param output_folder_:
        """

        fiduc_, spots_, repl_, params_ = None, None, None, None
        if c.METADATA_EXTENSION == 'xml':
            # check that exactly one .xml is in the data folder
            xml = [f for f in os.listdir(input_folder_) if '.xml' in f]
            if len(xml) > 1:
                raise IOError("more than one .xml file found, aborting")
            xml_path = os.path.join(input_folder_, xml[0])

            # parsing .xml
            fiduc_, spots_, repl_, params_ = txt_parser.create_xml_dict(xml_path)

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
            fiduc_, _, repl_, params_ = txt_parser.create_csv_dict(csv_paths)

        elif c.METADATA_EXTENSION == 'xlsx':
            pass

        # setting constants
        c.fiducials = fiduc_
        c.spots = spots_
        c.replicates = repl_

        c.params['rows'] = params_['rows']
        c.params['columns'] = params_['columns']
        c.params['v_pitch'] = params_['v_pitch']
        c.params['h_pitch'] = params_['h_pitch']
        c.params['spot_width'] = params_['spot_width']
        c.params['bg_offset'] = params_['bg_offset']
        c.params['bg_thickness'] = params_['bg_thickness']
        c.params['max_diam'] = params_['max_diam']
        c.params['min_diam'] = params_['min_diam']

        # setting constant arrays
        self._create_spot_id_array()
        self._create_spot_type_array()
        self._create_fiducials_array()
        self._create_antigen_array()

        # setting location of fiducials and other useful parameters
        self._calculate_fiduc_coords()
        self._calculate_fiduc_idx()
        self._calc_scienion_spot_dist()

        self._set_run_path(output_folder_)

    def _create_spot_id_array(self):
        self.spot_ids = np.empty(shape=(c.params['rows'], c.params['columns']), dtype='U100')
        c.SPOT_ID_ARRAY = txt_parser.populate_array_id(self.spot_ids, c.spots)

    def _create_spot_type_array(self):
        self.spot_type = np.empty(shape=(c.params['rows'], c.params['columns']), dtype='U100')
        c.SPOT_TYPE_ARRAY = txt_parser.populate_array_spots_type(self.spot_type, c.spots, c.fiducials)

    def _create_fiducials_array(self):
        self.fiducials_array = np.empty(shape=(c.params['rows'], c.params['columns']), dtype='U100')
        c.FIDUCIAL_ARRAY = txt_parser.populate_array_fiduc(self.fiducials_array, c.fiducials)

    def _create_antigen_array(self):
        """
        This is the only array creator that requires one of the above arrays
            multiple "spot-ids" can contain the same "antigen".  Thus it's required to pass the "SPOT_ID_ARRAY"
        :return:
        """
        self.antigen_array = np.empty(shape=(c.params['rows'], c.params['columns']), dtype='U100')
        c.ANTIGEN_ARRAY = txt_parser.populate_array_antigen(self.antigen_array, c.SPOT_ID_ARRAY, c.replicates)

    def _calculate_fiduc_coords(self):
        """
        calculate and set the fiducial coordinates like:
            FIDUCIALS = [(0, 0), (0, 1), (0, 5), (7, 0), (7, 5)]
            fiducial coordinates are labeled as "Reference, Diagnostic"
        :return:
        """
        x, y = np.where(c.FIDUCIAL_ARRAY == 'Reference, Diagnostic')
        c.FIDUCIALS = list(zip(x, y))

    def _calculate_fiduc_idx(self):
        """
        calculate fiducial index like
            FIDUCIALS_IDX = [0, 5, 6, 30, 35]
            FIDUCIALS_IDX_8COLS = [0, 7, 8, 40, 47]
        :return:
        """
        c.FIDUCIALS_IDX = list(np.where(c.FIDUCIAL_ARRAY.flatten() == 'Reference, Diagnostic')[0])

    def _calc_scienion_spot_dist(self):
        """
        calculate distance between spots in pixels for the scienion camera
        :return:
        """
        v_pitch_mm = c.params['v_pitch']
        h_pitch_mm = c.params['h_pitch']
        pix_size = c.params['pixel_size_scienion']

        v_pitch_pix = v_pitch_mm/pix_size
        h_pitch_pix = h_pitch_mm/pix_size
        c.SCENION_SPOT_DIST = np.mean([v_pitch_pix, h_pitch_pix]).astype('uint8')

    def _calc_octopi_spot_dist(self):
        """
        calculate distance between spots in pixels for the octopi camera
        :return:
        """
        v_pitch_mm = c.params['v_pitch']
        h_pitch_mm = c.params['h_pitch']
        pix_size = c.params['pixel_size_octopi']

        v_pitch_pix = v_pitch_mm / pix_size
        h_pitch_pix = h_pitch_mm / pix_size
        c.SCENION_SPOT_DIST = np.mean([v_pitch_pix, h_pitch_pix]).astype('uint8')

    # set filesaving run_path
    def _set_run_path(self, output_folder):
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
