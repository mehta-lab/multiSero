import numpy as np

import array_analyzer.extract.constants as constants
import array_analyzer.load.report as report


def test_report_init():
    antigen_array = np.empty(shape=(2, 3), dtype='U100')
    antigen_array[0, 0] = 'antigen_0_0'
    antigen_array[1, 2] = 'antigen_1_2'
    antigen_array[0, 1] = 'suuuuuuuuper_loooooooong_antigen_name'
    constants.ANTIGEN_ARRAY = antigen_array
    constants.RUN_PATH = 'test_run_dir'
    # Create instance
    reporter = report.ReportWriter()
    # Check paths
    assert reporter.od_path == 'test_run_dir/median_ODs.xlsx'
    assert reporter.int_path == 'test_run_dir/median_intensities.xlsx'
    assert reporter.bg_path == 'test_run_dir/median_backgrounds.xlsx'
    # Check that antigens have correct names and indices
    antigen_df = reporter.antigen_df
    assert antigen_df.shape == (3, 3)
    antigen = antigen_df.loc[(antigen_df['grid_row'] == 0) &
                             (antigen_df['grid_col'] == 0), 'antigen'].values[0]
    assert antigen == 'antigen_0_0'

    antigen = antigen_df.loc[(antigen_df['grid_row'] == 1) &
                             (antigen_df['grid_col'] == 2), 'antigen'].values[0]
    assert antigen == 'antigen_1_2'

    antigen = antigen_df.loc[(antigen_df['grid_row'] == 0) &
                             (antigen_df['grid_col'] == 1), 'antigen'].values[0]
    assert antigen == 'suuuuuuuuper_loooooooong_antige'
    # Check the dicts that will be reports
    assert list(reporter.report_int) == [
        'antigen_0_0',
        'suuuuuuuuper_loooooooong_antige',
        'antigen_1_2',
    ]
    # Check first dataframe
    plate_df = reporter.report_int['antigen_0_0']
    assert plate_df.shape == (12, 8)
    assert list(plate_df) == ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']


def test_antigen_df():
    antigen_array = np.empty(shape=(2, 3), dtype='U100')
    antigen_array[0, 0] = 'antigen_0_0'
    antigen_array[1, 2] = 'antigen_1_2'
    constants.ANTIGEN_ARRAY = antigen_array
    # Create instance
    reporter = report.ReportWriter()
    antigen_df = reporter.get_antigen_df()
    assert antigen_df.shape == (2, 3)
    assert list(antigen_df) == ['antigen', 'grid_row', 'grid_col']
    assert list(antigen_df['antigen'].values) == ['antigen_0_0', 'antigen_1_2']
