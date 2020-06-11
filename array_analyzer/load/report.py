import array_analyzer.extract.constants as constants
from copy import deepcopy
import numpy as np
import pandas as pd
import warnings


def write_od_to_plate(data, well_name, well_array):
    """
    Given data from image_parser.compute_od, write the output to an array
        representing its WELL position
    :param data: np.ndarray output from image_parser.compute_od
    :param well_name: str well name
    :param well_array: np.ndarray of 8x12 format (96 well plate)
    :return:
    """
    if well_name in constants.IMAGE_TO_WELL:
        (row, col) = constants.IMAGE_TO_WELL[well_name]
    else:
        raise AttributeError(f"well name {well_name} is not recognized")

    well_array[row-1, col-1] = data


def write_antigen_report(writer, well_array, array_type):
    """
    Creates and writes a single Excel worksheet containing:
        - one of OD, INT, or BG values for a SINGLE antigen, in row-col format
        - row-col format is based on 96-well plate: c.WELL_OUTPUT_TEMPLATE
    :param writer: pd.ExcelWriter object
    :param well_array: np.ndarray of 8x12 format (96 well plate), Specific to OD, INT, or BG
    :param array_type: str one of 'od', 'int', 'bg' based on well_array
    :return:
    """
    # todo: loops over all antigens for every well.  So if 6x6 arrays of antigens, 6x6x96 calls
    #  is there a more efficent way to do this?
    well_to_image = {v: k for k, v in constants.IMAGE_TO_WELL.items()}
    for antigen_position, antigen in np.ndenumerate(constants.ANTIGEN_ARRAY):
        if antigen == '' or antigen is None:
            continue
        if constants.DEBUG:
            print(f"writing antigen {antigen} to excel sheets")

        sheet = deepcopy(constants.WELL_OUTPUT_TEMPLATE)

        # loop all wells and write OD, INT or BG of this antigen
        for position, well in np.ndenumerate(well_array):
            # extract intensity at antigen_position within this well
            if well is None:
                val = None
            else:
                val = well[antigen_position[0], antigen_position[1]]

            # get the well name describing this position in form "A#"
            well_name = well_to_image[(position[0] + 1, position[1] + 1)]

            # write a new sheet named "A#"
            sheet[well_name[0]][int(well_name[1:])] = val

        # write the outputs to a worksheet
        sheet_df = pd.DataFrame(sheet).T

        sheet_name = f'{array_type}_{antigen_position[0]}_{antigen_position[1]}_{antigen}'
        if len(sheet_name) >= 31:
            warnings.warn("antigen sheet name is too long, truncating")
            sheet_name = sheet_name[:31]

        sheet_df.to_excel(writer,
                          sheet_name=sheet_name)

