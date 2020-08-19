import numpy as np
import pandas as pd


def antigen2D_to_df1D(xlsx_path, sheet, data_col):
    """
    Convert old 2D output format (per antigen) to 1D dataframe
    :param xlsx_path:
    :param sheet:
    :param data_col:
    :return:
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet, index_col=0)
    df = df.unstack().reset_index(name=data_col)  # unpivot (linearize) the table
    df.rename(columns={'level_1': 'antigen_row', 'level_0': 'antigen_col'}, inplace=True)
    df[['antigen_row', 'antigen_col']] = df[['antigen_row', 'antigen_col']].applymap(int)
    df = df[['antigen_row', 'antigen_col', data_col]]
    df.dropna(inplace=True)
    return df


def well2D_to_df1D(xlsx_path, sheet, data_col):
    """
    Convert new 2D output format (per well) to 1D dataframe
    :param xlsx_path:
    :param sheet:
    :param data_col:
    :return:
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet, index_col=0)
    df = df.unstack().reset_index(name=data_col)  # unpivot (linearize) the table
    df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
    df['well_id'] = df.row_id + df.col_id.map(str)
    df = df[['well_id', data_col]]
    return df


def read_plate_info(metadata_xlsx):
    print('Reading the plate info...')

    sheet_names = ['serum ID',
                   'serum dilution',
                   'serum type',
                   'secondary ID',
                   'secondary dilution']
    plate_info_df = pd.DataFrame()
    # get sheet names that are available in metadata
    sheet_names = list(set(metadata_xlsx.sheet_names).intersection(sheet_names))
    for sheet_name in sheet_names:
        sheet_df = pd.read_excel(metadata_xlsx, sheet_name=sheet_name, index_col=0)
        sheet_df = sheet_df.unstack().reset_index(name=sheet_name)  # unpivot (linearize) the table
        sheet_df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
        if plate_info_df.empty:
            plate_info_df = sheet_df
        else:
            plate_info_df = pd.merge(plate_info_df,
                                     sheet_df,
                                     how='left', on=['row_id', 'col_id'])
    plate_info_df['well_id'] = plate_info_df.row_id + plate_info_df.col_id.map(str)
    sheet_names.append('well_id')
    # convert to number and non-numeric to NaN
    plate_info_df['serum dilution'] = \
        plate_info_df['serum dilution'].apply(pd.to_numeric, errors='coerce')
    plate_info_df.dropna(inplace=True)
    if np.all(plate_info_df['serum dilution'] >= 1):
        # convert dilution to concentration
        plate_info_df['serum dilution'] = 1 / plate_info_df['serum dilution']
    plate_info_df.drop(['row_id', 'col_id'], axis=1, inplace=True)
    return plate_info_df


def read_antigen_info(metadata_path):
    print('Reading antigen information...')
    antigen_df = antigen2D_to_df1D(xlsx_path=metadata_path, sheet='antigen_array', data_col='antigen')
    antigen_type_df = antigen2D_to_df1D(xlsx_path=metadata_path, sheet='antigen_type', data_col='antigen type')
    antigen_df = pd.merge(antigen_df, antigen_type_df, how='left', on=['antigen_row', 'antigen_col'])
    return antigen_df


def read_pysero_output(file_path, antigen_df, file_type='od'):
    print('Reading {}...'.format(file_type))
    data_col = {'od': 'OD', 'int': 'intensity', 'bg': 'background'}
    data_df = pd.DataFrame()

    with pd.ExcelFile(file_path) as file:
        for _, row in antigen_df.iterrows():
            sheet_name = '{}_{}_{}_{}'.format(file_type, row['antigen_row'], row['antigen_col'], row['antigen'])
            data_1_antiten_df = well2D_to_df1D(xlsx_path=file, sheet=sheet_name, data_col=data_col[file_type])
            data_1_antiten_df['antigen_row'] = row['antigen_row']
            data_1_antiten_df['antigen_col'] = row['antigen_col']
            data_1_antiten_df['antigen'] = row['antigen']
            data_df = data_df.append(data_1_antiten_df, ignore_index=True)
    return data_df


def read_scn_output(file_path, plate_info_df):
    # Read analysis output from Scienion
    scienion_df = pd.DataFrame()
    with pd.ExcelFile(file_path) as scienion_xlsx:
        for well_id in plate_info_df['well_id']:
            OD_1_antiten_df = pd.read_excel(scienion_xlsx, sheet_name=well_id)
            OD_1_antiten_df['well_id'] = well_id
            scienion_df = scienion_df.append(OD_1_antiten_df, ignore_index=True)
    # parse spot ids
    spot_id_df = scienion_df['ID'].str.extract(r'spot-(\d)-(\d)')
    spot_id_df = spot_id_df.astype(int) - 1  # index starting from 0
    spot_id_df.rename(columns={0: 'antigen_row', 1: 'antigen_col'}, inplace=True)
    scienion_df = pd.concat([spot_id_df, scienion_df], axis=1)
    scienion_df.drop('ID', axis=1, inplace=True)
    # invert the intensity and compute ODs
    df_scn = scienion_df.loc[:, ['antigen_row', 'antigen_col', 'well_id']]
    df_scn['intensity'] = 1 - scienion_df['Median'] / 255
    df_scn['background'] = 1 - scienion_df['Background Median'] / 255
    df_scn['OD'] = np.log10(df_scn['background'] / df_scn['intensity'])
    return df_scn


def slice_df(df, slice_action, column, keys):
    if slice_action is None:
        return df
    elif slice_action == 'keep':
        df = df[df[column].isin(keys)]
    elif slice_action == 'drop':
        df = df[~df[column].isin(keys)]
    else:
        raise ValueError('slice action has to be "keep" or "drop", not "{}"'.format(slice_action))
    return df


def normalize_od_helper(norm_antigen):
    def normalize(df):
        norm_antigen_df = slice_df(df, 'keep', 'antigen', [norm_antigen])
        norm_factor = norm_antigen_df['OD'].mean()
        df['OD'] = df['OD'] / norm_factor
        return df
    return normalize


def normalize_od(df, norm_antigen=None, group='plate'):
    """fit model to x, y data in dataframe.
    Return a dataframe with fit x, y for plotting
    """
    if norm_antigen is None:
        return df
    if group == 'plate':
        groupby_cols = ['plate_id']
    elif group == 'well':
        groupby_cols = ['plate_id', 'well_id']
    else:
        ValueError('normalization group has to be plate or well, not {}'.format(group))
    norm_antigen_df = slice_df(df, 'keep', 'antigen', [norm_antigen])
    df.loc[df['antigen'] == norm_antigen, 'OD'] = norm_antigen_df['OD'] / norm_antigen_df['OD'].mean()
    norm_fn = normalize_od_helper(norm_antigen)
    df = df.groupby(groupby_cols).apply(norm_fn)
    # df = df.reset_index()
    return df


def offset_od_helper(norm_antigen):
    def offset(df):
        norm_antigen_df = slice_df(df, 'keep', 'antigen', [norm_antigen])
        norm_factor = norm_antigen_df['OD'].mean()
        df['OD'] = df['OD'] - norm_factor
        df.loc[df['OD'] < 0, 'OD'] = 0
        return df
    return offset


def offset_od(df, norm_antigen=None, group='plate'):
    """fit model to x, y data in dataframe.
    Return a dataframe with fit x, y for plotting
    """
    if norm_antigen is None:
        return df
    if group == 'plate':
        groupby_cols = ['plate_id']
    elif group == 'well':
        groupby_cols = ['plate_id', 'well_id']
    else:
        ValueError('normalization group has to be plate or well, not {}'.format(group))
    norm_fn = offset_od_helper(norm_antigen)
    df = df.groupby(groupby_cols).apply(norm_fn)
    return df