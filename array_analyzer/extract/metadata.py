import os
import numpy as np

import array_analyzer.extract.txt_parser as txt_parser
import array_analyzer.extract.constants as c


class MetaData:

    def __init__(self, input_folder_, output_folder_):
        xml = [f for f in os.listdir(input_folder_) if '.xml' in f]
        if len(xml) > 1:
            raise IOError("more than one .xml file found, aborting")
        xml_path = input_folder_ + os.sep + xml[0]

        # parsing .xml
        fiduc_, spots_, repl_, params_ = txt_parser.create_xml_dict(xml_path)

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

        self._create_spot_id_array()
        self._create_spot_type_array()
        self._create_fiducials_array()
        self._create_antigen_array()

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

    def _calc_fiducial_coords(self):
        """
        fiducials are defined as ""
        :return:
        """

    # calculate and set the fiducial coordinates like:
    # FIDUCIALS = [(0, 0), (0, 1), (0, 5), (7, 0), (7, 5)]
    def calc_fiduc_coords(self):
        pass

    # spot distance for ICP (in pixels?)
    def calc_scienion_spot_dist(self):
        # use vpitch, hpitch, pixel size in params to calculate and set spot dist
        pass

