import array_analyzer.extract.constants as c
from copy import deepcopy
import numpy as np
import pandas as pd


def write_od_to_plate(data, well_name, array_type):
    """
    Given data from image_parser.compute_od, write the output to an array
        representing its WELL position
    :param data: np.ndarray output from image_parser.compute_od
    :param well_name: str well name
    :param array_type: str data type
    :return:
    """
    if array_type not in ['od', 'int', 'bg']:
        raise AttributeError(f"array type {array_type} not implemented!")

    if well_name in c.IMAGE_TO_WELL:
        (row, col) = c.IMAGE_TO_WELL[well_name]
    else:
        raise AttributeError(f"well name {well_name} is not recognized")

    if array_type == 'od':
        c.WELL_OD_ARRAY[row-1, col-1] = data
    if array_type == 'int':
        c.WELL_INT_ARRAY[row-1, col-1] = data
    if array_type == 'bg':
        c.WELL_BG_ARRAY[row-1, col-1] = data


def write_antigen_report(writer, array_type):
    well_to_image = {v: k for k, v in c.IMAGE_TO_WELL.items()}
    for antigen_position, antigen in np.ndenumerate(c.ANTIGEN_ARRAY):
        if antigen == '' or antigen is None:
            continue
        print(f"writing antigen {antigen} to excel sheets")

        sheet = deepcopy(c.WELL_OUTPUT_TEMPLATE)

        # loop all wells and write OD, INT, BG of this antigen
        if array_type == 'od':
            for od_position, od_well in np.ndenumerate(c.WELL_OD_ARRAY):
                od_val = od_well[antigen_position[0], antigen_position[1]]
                well_name = well_to_image[(od_position[0] + 1, od_position[1] + 1)]
                sheet[well_name[0]][int(well_name[1:])] = od_val
        elif array_type == 'int':
            for int_position, int_well in np.ndenumerate(c.WELL_INT_ARRAY):
                int_val = int_well[antigen_position[0], antigen_position[1]]
                well_name = well_to_image[(int_position[0] + 1, int_position[1] + 1)]
                sheet[well_name[0]][int(well_name[1:])] = int_val
        elif array_type == 'bg':
            for bg_position, bg_well in np.ndenumerate(c.WELL_BG_ARRAY):
                bg_val = bg_well[antigen_position[0], antigen_position[1]]
                well_name = well_to_image[(bg_position[0] + 1, bg_position[1] + 1)]
                sheet[well_name[0]][int(well_name[1:])] = bg_val
        else:
            raise AttributeError(f"report array type {array_type} not supported")

        # write the outputs from the above three to a worksheet
        od_sheet_df = pd.DataFrame(sheet).T

        od_sheet_df.to_excel(writer, sheet_name=f'{array_type}_{antigen_position}_{antigen}')
