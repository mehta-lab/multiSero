import array_analyzer.extract.constants as constants
from copy import deepcopy
import numpy as np
import pandas as pd
import warnings


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

    if well_name in constants.IMAGE_TO_WELL:
        (row, col) = constants.IMAGE_TO_WELL[well_name]
    else:
        raise AttributeError(f"well name {well_name} is not recognized")

    if data is None:
        print("\tDATA IS NONE FOR THIS WELL")
        data = np.ones(shape=(constants.params['rows'], constants.params['columns']))

    if array_type == 'od':
        constants.WELL_OD_ARRAY[row-1, col-1] = data
    if array_type == 'int':
        constants.WELL_INT_ARRAY[row-1, col-1] = data
    if array_type == 'bg':
        constants.WELL_BG_ARRAY[row-1, col-1] = data


def write_antigen_report(writer, array_type):
    """
    Creates and writes a single Excel worksheet containing:
        - one of OD, INT, or BG values for a SINGLE antigen, in row-col format
        - row-col format is based on 96-well plate: c.WELL_OUTPUT_TEMPLATE
    :param writer: pd.ExcelWriter object
    :param array_type: str one of 'od', 'int', 'bg'
    :return:
    """
    #todo: loops over all antigens for every well.  So if 6x6 arrays of antigens, 6x6x96 calls
    # is there a more efficent way to do this?
    well_to_image = {v: k for k, v in constants.IMAGE_TO_WELL.items()}
    for antigen_position, antigen in np.ndenumerate(constants.ANTIGEN_ARRAY):
        if antigen == '' or antigen is None:
            continue
        if constants.DEBUG:
            print(f"writing antigen {antigen} to excel sheets")

        sheet = deepcopy(constants.WELL_OUTPUT_TEMPLATE)
        # loop all wells and write OD, INT, BG of this antigen
        if array_type == 'od':
            sheet = write_to_sheet(sheet, constants.WELL_OD_ARRAY, antigen_position, well_to_image)
        elif array_type == 'int':
            sheet = write_to_sheet(sheet, constants.WELL_INT_ARRAY, antigen_position, well_to_image)
        elif array_type == 'bg':
            sheet = write_to_sheet(sheet, constants.WELL_BG_ARRAY, antigen_position, well_to_image)
        else:
            raise AttributeError(f"report array type {array_type} not supported")

        # write the outputs from ONE of the above three to a worksheet
        # (this function is called once for each OD, INT, BG)
        od_sheet_df = pd.DataFrame(sheet).T

        sheet_name = f'{array_type}_{antigen_position[0]}_{antigen_position[1]}_{antigen}'
        if len(sheet_name) >= 31:
            warnings.warn("antigen sheet name is too long, truncating")
            sheet_name = sheet_name[:31]

        od_sheet_df.to_excel(writer,
                             sheet_name=sheet_name)


def write_to_sheet(sheet_, well_array_, antigen_position_, well_to_image_):
    """

    :param sheet_:
    :param well_array_:
    :param antigen_position_:
    :param well_to_image_:
    :return:
    """
    # iterate over all 96 well positions
    for position, well in np.ndenumerate(well_array_):
        # extract intensity at antigen_position within this well
        if well is None:
            val = -999
        else:
            val = well[antigen_position_[0], antigen_position_[1]]

        # get the well name describing this position in form "A#"
        well_name = well_to_image_[(position[0] + 1, position[1] + 1)]

        # write a new sheet named "A#"
        sheet_[well_name[0]][int(well_name[1:])] = val

    return sheet_
