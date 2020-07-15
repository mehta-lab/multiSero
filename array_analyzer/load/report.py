import collections
from copy import deepcopy
import logging
import numpy as np
import os
import pandas as pd

import array_analyzer.extract.constants as constants


class ReportWriter:

    def __init__(self):
        self.logger = logging.getLogger(constants.LOG_NAME)
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
                self.logger.warning("antigen name is too long for sheet, truncating")
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
        self.logger.debug('Loaded existing reports')

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
        self.logger.debug("Wrote plate reports")
