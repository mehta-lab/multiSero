import os
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

