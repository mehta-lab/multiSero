import collections
from copy import deepcopy
import numpy as np
import os
import pandas as pd

import array_analyzer.extract.constants as constants


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


def write_antigen_report(writer, well_array, array_type, logger):
    """
    Creates and writes a single Excel worksheet containing:
        - one of OD, INT, or BG values for a SINGLE antigen, in row-col format
        - row-col format is based on 96-well plate: c.WELL_OUTPUT_TEMPLATE
    :param writer: pd.ExcelWriter object
    :param well_array: np.ndarray of 8x12 format (96 well plate), Specific to OD, INT, or BG
    :param array_type: str one of 'od', 'int', 'bg' based on well_array
    :param logging inst logger: Log total distance in each iteration
    """
    # todo: loops over all antigens for every well.  So if 6x6 arrays of antigens, 6x6x96 calls
    #  is there a more efficent way to do this?
    logger.info("Writing array {}".format(array_type))
    well_to_image = {v: k for k, v in constants.IMAGE_TO_WELL.items()}
    for antigen_position, antigen in np.ndenumerate(constants.ANTIGEN_ARRAY):
        if antigen == '' or antigen is None:
            continue
        logger.debug(f"writing antigen {antigen} to excel sheets")

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
            logger.warning("antigen sheet name is too long, truncating")
            sheet_name = sheet_name[:31]

        sheet_df.to_excel(writer, sheet_name=sheet_name)


class ReportWriter:

    def __init__(self):
        # Dataframe for a whole plate
        cols = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        rows = list(range(1, 13))
        plate_df = pd.DataFrame(None, index=rows, columns=cols)

        report_dict = collections.OrderedDict()
        # Dataframe for antigen positions on grid
        self.antigen_df = pd.DataFrame(columns=['antigen', 'grid_row', 'grid_col'])
        for antigen_position, antigen in np.ndenumerate(constants.ANTIGEN_ARRAY):
            if antigen == '' or antigen is None:
                continue
            # Abbreviate antigen name if too long
            sheet_name = antigen
            if len(sheet_name) >= 31:
                # logger.warning("antigen name is too long for sheet, truncating")
                sheet_name = sheet_name[:31]
            # Add empty plate dataframe to antigen
            report_dict[sheet_name] = plate_df.copy()
            # Add grid row and column for antigen
            idx_row = {'antigen': sheet_name,
                       'grid_row': antigen_position[0],
                       'grid_col': antigen_position[1]}
            self.antigen_df = self.antigen_df.append(idx_row, ignore_index=True)

        self.antigen_names = list(self.antigen_df['antigen'].values)
        self.report_int = deepcopy(report_dict)
        self.report_bg = deepcopy(report_dict)
        self.report_od = deepcopy(report_dict)
        # Report paths
        self.od_path = os.path.join(constants.RUN_PATH, 'median_ODs.xlsx')
        self.int_path = os.path.join(constants.RUN_PATH, 'median_intensities.xlsx')
        self.bg_path = os.path.join(constants.RUN_PATH, 'median_backgrounds.xlsx')

    def get_antigen_df(self):
        """
        Returns dataframe with antigen names and their locations on the grid.

        :return pd.DataFrame antigen_df: Antigen names and grid rows, cols
        """
        return self.antigen_df

    def load_existing_reports(self):
        """
        If doing a rerun, load existing reports and make sure the sheet names
        (keys) match the ones from the current run.
        Rerun wells will be added to existing reports and rewritten.
        """
        assert os.path.isfile(self.od_path), \
            "OD report doesn't exist: {}".format(self.od_path)
        assert os.path.isfile(self.int_path), \
            "Intensity report doesn't exist: {}".format(self.int_path)
        assert os.path.isfile(self.bg_path), \
            "Background report doesn't exist: {}".format(self.bg_path)
        # Read reports and make sure they have the right keys
        ordered_dict = pd.read_excel(self.od_path, sheet_name=None, index_col=0)
        assert list(ordered_dict) == self.antigen_names, \
            "Existing report keys don't match current keys"
        self.report_od = ordered_dict
        ordered_dict = pd.read_excel(self.int_path, sheet_name=None, index_col=0)
        assert list(ordered_dict) == self.antigen_names, \
            "Existing report keys don't match current keys"
        self.report_int = ordered_dict
        ordered_dict = pd.read_excel(self.bg_path, sheet_name=None, index_col=0)
        assert list(ordered_dict) == self.antigen_names, \
            "Existing report keys don't match current keys"
        self.report_bg = ordered_dict

    def assign_well_to_plate(self, well_name, spots_df):
        """
        Takes intensity, background and OD values for a well and
        reassigns them by antigens for the plate.

        :param str well_name: Well name (e.g. 'B12')
        :param pd.DataFrame spots_df: Metrics for all spots in a well
        """
        plate_col = well_name[0]
        plate_row = int(well_name[1:])
        for idx, antigen_row in self.antigen_df.iterrows():
            spots_row = spots_df.loc[
                (spots_df['grid_row'] == antigen_row['grid_row']) &
                (spots_df['grid_col'] == antigen_row['grid_col']),
            ]
            self.report_int[antigen_row['antigen']].at[plate_row, plate_col] = \
                spots_row['intensity_median'].values[0]
            self.report_bg[antigen_row['antigen']].at[plate_row, plate_col] = \
                spots_row['bg_median'].values[0]
            self.report_od[antigen_row['antigen']].at[plate_row, plate_col] = \
                spots_row['od_norm'].values[0]

    def write_reports(self):
        """
        After all wells are run, write plate based reports for OD,
        intensity, and background.
        """
        # Write OD report
        with pd.ExcelWriter(self.od_path) as writer:
            for antigen_name in self.antigen_names:
                sheet_df = self.report_od[antigen_name]
                sheet_df.to_excel(writer, sheet_name=antigen_name)
        # Write intensity report
        with pd.ExcelWriter(self.int_path) as writer:
            for antigen_name in self.antigen_names:
                sheet_df = self.report_int[antigen_name]
                sheet_df.to_excel(writer, sheet_name=antigen_name)
        # Write background report
        with pd.ExcelWriter(self.bg_path) as writer:
            for antigen_name in self.antigen_names:
                sheet_df = self.report_bg[antigen_name]
                sheet_df.to_excel(writer, sheet_name=antigen_name)
