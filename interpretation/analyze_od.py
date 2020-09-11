import pandas as pd
import numpy as np
import os
import matplotlib
from matplotlib import pyplot as plt
from natsort import natsorted
import seaborn as sns
from interpretation.plotting import fourPL, fit2df, get_roc_df, roc_plot_grid, thr_plot_grid
from interpretation.report_reader import read_plate_info, read_antigen_info, read_pysero_output, read_scn_output, \
    slice_df, normalize_od, offset_od

sns.set_context("talk")
font = {'size': 10, 'weight': 'normal', 'family': 'arial'}
matplotlib.rc('font', **font)
#%% Scienion paths
scn_psr_dirs = [r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200630_1647',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200630_1657',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200630_1708',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200630_1739',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-25-28-COVID_June24_OJassay_plate9_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200701_0933',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-25-28-COVID_June24_OJassay_plate9_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_1002',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-25-28-COVID_June24_OJassay_plate9_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_1028',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-40-02-COVID_June5_OJassay_plate7_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200611_1257',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-40-02-COVID_June5_OJassay_plate7_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200613_1341',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-40-02-COVID_June5_OJassay_plate7_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200709_1024',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-44-32-COVID_June5_OJassay_plate8_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200613_2017',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-44-32-COVID_June5_OJassay_plate8_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200613_2042',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-44-32-COVID_June5_OJassay_plate8_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200613_2051',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-44-32-COVID_June5_OJassay_plate8_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200618_1449',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-21-32-COVID_June24_OJassay_plate4_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200630_1756',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-21-32-COVID_June24_OJassay_plate4_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_0849',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-21-32-COVID_June24_OJassay_plate4_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_0903',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-29-42-COVID_June24_OJassay_plate10_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200701_1047',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-29-42-COVID_June24_OJassay_plate10_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_1110',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-29-42-COVID_June24_OJassay_plate10_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_1141',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-08-14-18-24-59-COVID_August14_OJassay_plate11_images/Stitched data from multiple outputs/pysero_biotin_fiducial_20200818_1526',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-08-14-18-24-59-COVID_August14_OJassay_plate11_images/Stitched data from multiple outputs/pysero_IgG fiducial_20200818_1547'
                ]

scn_scn_dirs = [r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-25-28-COVID_June24_OJassay_plate9_images/',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-40-02-COVID_June5_OJassay_plate7_images/',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-44-32-COVID_June5_OJassay_plate8_images/',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-21-32-COVID_June24_OJassay_plate4_images/',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-29-42-COVID_June24_OJassay_plate10_images/',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-08-14-18-24-59-COVID_August14_OJassay_plate11_images',
                ]

output_dir = r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/OJ_plate3_9_7_8_4_10_11'
fig_path = os.path.join(output_dir, 'pysero_plots')
scn_psr_plate_ids = ['plate_3'] * 4 + ['plate_9'] * 3 + ['plate_7'] * 3 + ['plate_8'] * 4 + ['plate_4'] * 3 + ['plate_10'] * 3 + ['plate_11'] * 2
scn_scn_plate_ids = ['plate_3', 'plate_9', 'plate_7', 'plate_8', 'plate_4', 'plate_10', 'plate_11']
scn_psr_slice_actions = ['keep', 'drop', 'keep', 'drop'] + \
                        [None, 'drop', 'keep'] + \
                        [None, 'drop', 'keep'] + \
                        [None, 'keep', 'drop', 'keep'] + \
                        [None, 'drop', 'keep'] + \
                        [None, 'drop', 'keep'] + \
                        [None,None]
scn_psr_well_ids = [['E3', 'E4'], ['E3', 'E4'], ['F2'], ['F2']] + \
                   [None, ['B7', 'B8', 'B10','B11','D3','D8','D9','D10','D11','F8','F9','F11','H2','H7','H8','H9','H10','H11'],
                ['B7', 'B8', 'B10','B11','D3','D8','D9','D10','D11','F8','F9','F11','H2','H7','H8','H9','H10','H11']] + \
                   [None, ['B4','F2'], ['B4','F2']] + \
                   [None, ['B1','D4','D6','F5'], ['B1', 'B12','D4','D6','F5','H11'], ['B12','H11']] + \
                   [None, ['B3','B5','D1','D2','D5','F5'], ['B3','B5','D1','D2','D5','F5']] + \
                   [None, ['B12','D2','D3','D1','D5','H2','H5'], ['B12','D1','D5','D2','D3','H2']] + \
                    [None,None]
sera_fit_list = ['Pool', 'mab', 'CR3022']
sera_cat_list = ['Pool', 'mab', 'Blank', 'CR3022']
sera_roc_list = sera_cat_list
#%% Nautilus

ntl_dirs = [r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate10_images_655_2020-06-24 18-37-39.863642/0 renamed 16bit/biotin_fiducial/pysero_biotin_fiducial_20200708_1813',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate10_images_655_2020-06-24 18-37-39.863642/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200708_1853',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate10_images_655_2020-06-24 18-37-39.863642/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200709_1431',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate10_images_655_2020-06-24 18-37-39.863642/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200709_1613',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate10_images_655_2020-06-24 18-37-39.863642/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200709_1516',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate3_images_655_2020-06-24 17-58-27.941906/0 renamed 16bit/biotin_fiducial/pysero_biotin_fiducial_20200702_1350',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate3_images_655_2020-06-24 17-58-27.941906/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200727_1954',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate3_images_655_2020-06-24 17-58-27.941906/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200706_0654',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate3_images_655_2020-06-24 17-58-27.941906/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200728_1129',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate9_images_655_2020-06-24 18-26-58.173049/0 renamed/Stitched data/pysero_biotin_fiducial_20200706_0847',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate9_images_655_2020-06-24 18-26-58.173049/0 renamed/Stitched data/pysero_igg_fiducial_20200706_0903',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate9_images_655_2020-06-24 18-26-58.173049/0 renamed/Stitched data/pysero_igg_fiducial_20200706_0937',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate9_images_655_2020-06-24 18-26-58.173049/0 renamed/Stitched data/pysero_igg_fiducial_20200708_1136',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate7_images_2020-06-08 19-54-45.277504/0 renamed/biotin fiducial/pysero_biotin fiducial_20200727_2016',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate7_images_2020-06-08 19-54-45.277504/0 renamed/igg fiducial/pysero_igg fiducial_20200727_2038',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate7_images_2020-06-08 19-54-45.277504/0 renamed/igg fiducial/pysero_igg fiducial_20200728_0822',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate7_images_2020-06-08 19-54-45.277504/0 renamed/igg fiducial/pysero_igg fiducial_20200728_0902',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate8_images_2020-06-08 20-38-40.554230/0 renamed/biotin_fiducial/pysero_biotin_fiducial_20200728_1018',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate8_images_2020-06-08 20-38-40.554230/0 renamed/igg_fiducial/pysero_igg_fiducial_20200728_1044',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate8_images_2020-06-08 20-38-40.554230/0 renamed/igg_fiducial/pysero_igg_fiducial_20200728_1108',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate4_images_655_2020-06-24 18-13-40.894052/0 renamed/biotin_fiducial/pysero_biotin_fiducial_20200728_0921',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate4_images_655_2020-06-24 18-13-40.894052/0 renamed/igg_fiducial/pysero_igg_fiducial_20200709_1411',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate4_images_655_2020-06-24 18-13-40.894052/0 renamed/igg_fiducial/pysero_igg_fiducial_20200728_0956',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-08-14-COVID_Aug14_OJ_Plate11_2020-08-14 19-29-59.049679/pysero_0_renamed_16bit_20200819_1557']

ntl_plate_ids = ['plate_10'] * 5 + ['plate_3'] * 4 + ['plate_9'] * 4 + ['plate_7'] * 4 + ['plate_8'] * 3 + ['plate_4'] * 3 + ['plate_11']

ntl_slice_actions = [None, 'drop','drop','keep','keep'] +\
                  [None,'drop','keep','keep'] +\
                  [None, 'drop','drop','keep'] +\
                  [None, 'keep', 'drop', 'keep'] +\
                  [None,'drop','keep'] +\
                   [None,'drop','keep'] +\
                    [None]

ntl_well_ids = [None,['B2', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11'],['B1', 'B2', 'B4', 'B12', 'D2', 'D4', 'D5', 'D8', 'F3','H4','B3','B5','D1','D2','D3','D12','F1','F2','F4','F5','F12','H1','H2','H3','H12'],['B12', 'H4'],['B2','D8']] + \
           [None, ['D3', 'D4','D7','D8','D9','D10', 'D11', 'F7', 'F8', 'F9','F10', 'F12','H7','H9','H10','H11'],['D3', 'D8','D9','D10','D11','F8','F9','F10','H7','H9','H10','H11'], ['D4','D7','F7','F12']] + \
            [None, ['B1',' B2', 'B5','B7','B8','b9','b10','B11','D3','D5','D7','D8','D9','D10','D11','F7','F8','F9','F10','F11','H7','H8','H9','H10','H11'], ['B5','D5','D6'],['B5','D2','D6']] + \
            [None, ['B2','B12','D3','D12', 'F4','F12','H1','H2','H3','H6'], ['B2','B12','D3','D12', 'F4','F12','H1','H2','H3','H6','B7','B8','B9','B10','B11','D7','D8','D9','D10','D11','F2','F6','F7','F8','F9','F10','F11','H4','H7','H8','H9','H10','H11','H12'], ['B7','B8','B9','B10','B11','D7','D8','D9','D10','D11','F2','F6','F7','F8','F9', 'F10','F11','H4','H7','H8','H9','H10','H11','H12']] + \
            [None, ['B9','B10','B11','B12','D7','D8','D9','D11','F6','F7','F8','F9','F10','F11','H6','H7','H8','H9', 'H11','H12'], ['B9','B10','B11','D7','D8','D9','D11','F6','F7','F8','F9','F10','F11','H6','H7','H8','H9', 'H11','H12']] + \
            [None,['D1','D5','D12','H2'], ['D1','D5','D12','H2']] + \
            [None]
#%%
os.makedirs(fig_path, exist_ok=True)
load_master_report = True
if not load_master_report:
    df_list = []
    scn_df = pd.DataFrame()

    for scn_scn_dir, plate_id, in zip(scn_scn_dirs, scn_scn_plate_ids):
        metadata_path = os.path.join(scn_scn_dir, 'pysero_output_data_metadata.xlsx')
        with pd.ExcelFile(metadata_path) as meta_file:
            antigen_df = read_antigen_info(meta_file)
            plate_info_df = read_plate_info(meta_file)
        plate_info_df['plate_id'] = plate_id
        scn_fname = [f for f in os.listdir(scn_scn_dir) if '_analysis.xlsx' in f]
        scn_path = os.path.join(scn_scn_dir, scn_fname[0])
        scn_df_tmp = read_scn_output(scn_path, plate_info_df)
        # Join Scienion data with plateInfo
        scn_df_tmp = pd.merge(scn_df_tmp,
                          antigen_df,
                          how='left', on=['antigen_row', 'antigen_col'])
        scn_df_tmp = pd.merge(scn_df_tmp,
                          plate_info_df,
                          how='right', on=['well_id'])
        scn_df = scn_df.append(scn_df_tmp, ignore_index=True)
    scn_df['pipeline'] = 'scienion'
    scn_df.dropna(subset=['OD'], inplace=True)

    for data_folder, slice_action, well_id, plate_id in \
            zip(scn_psr_dirs, scn_psr_slice_actions, scn_psr_well_ids, scn_psr_plate_ids):
        print('Load {}...'.format(data_folder))
        metadata_path = os.path.join(data_folder, 'pysero_output_data_metadata.xlsx')
        OD_path = os.path.join(data_folder, 'median_ODs.xlsx')
        int_path = os.path.join(data_folder, 'median_intensities.xlsx')
        bg_path = os.path.join(data_folder, 'median_backgrounds.xlsx')

        with pd.ExcelFile(metadata_path) as meta_file:
            antigen_df = read_antigen_info(meta_file)
            plate_info_df = read_plate_info(meta_file)
        plate_info_df['plate_id'] = plate_id
        OD_df = read_pysero_output(OD_path, antigen_df, file_type='od')
        int_df = read_pysero_output(int_path, antigen_df, file_type='int')
        bg_df = read_pysero_output(bg_path, antigen_df, file_type='bg')
        OD_df = pd.merge(OD_df,
                         antigen_df[['antigen_row', 'antigen_col', 'antigen type']],
                         how='left', on=['antigen_row', 'antigen_col'])
        OD_df = pd.merge(OD_df,
                         plate_info_df,
                         how='right', on=['well_id'])
        pysero_df = pd.merge(OD_df,
                             int_df,
                             how='left', on=['antigen_row', 'antigen_col', 'well_id'])
        pysero_df = pd.merge(pysero_df,
                             bg_df,
                             how='left', on=['antigen_row', 'antigen_col', 'well_id'])
        pysero_df['pipeline'] = 'python'
        # pysero_df = pysero_df.append(scn_df)
        pysero_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        pysero_df.dropna(subset=['OD'], inplace=True)
        pysero_df = slice_df(pysero_df, slice_action, 'well_id', well_id)
        df_list.append(pysero_df)
    #%
    for data_folder, slice_action, well_id, plate_id in \
            zip(ntl_dirs, ntl_slice_actions, ntl_well_ids, ntl_plate_ids):
        print('Load {}...'.format(data_folder))
        metadata_path = os.path.join(data_folder, 'pysero_output_data_metadata.xlsx')
        OD_path = os.path.join(data_folder, 'median_ODs.xlsx')
        int_path = os.path.join(data_folder, 'median_intensities.xlsx')
        bg_path = os.path.join(data_folder, 'median_backgrounds.xlsx')

        with pd.ExcelFile(metadata_path) as meta_file:
            antigen_df = read_antigen_info(meta_file)
            plate_info_df = read_plate_info(meta_file)
        plate_info_df['plate_id'] = plate_id
        OD_df = read_pysero_output(OD_path, antigen_df, file_type='od')
        int_df = read_pysero_output(int_path, antigen_df, file_type='int')
        bg_df = read_pysero_output(bg_path, antigen_df, file_type='bg')
        OD_df = pd.merge(OD_df,
                         antigen_df[['antigen_row', 'antigen_col', 'antigen type']],
                         how='left', on=['antigen_row', 'antigen_col'])
        OD_df = pd.merge(OD_df,
                         plate_info_df,
                         how='right', on=['well_id'])
        pysero_df = pd.merge(OD_df,
                             int_df,
                             how='left', on=['antigen_row', 'antigen_col', 'well_id'])
        pysero_df = pd.merge(pysero_df,
                             bg_df,
                             how='left', on=['antigen_row', 'antigen_col', 'well_id'])
        pysero_df['pipeline'] = 'nautilus'
        pysero_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        pysero_df.dropna(subset=['OD'], inplace=True)
        pysero_df = slice_df(pysero_df, slice_action, 'well_id', well_id)
        df_list.append(pysero_df)
    df_list.append(scn_df)
    #% Concatenate dataframes
    stitched_pysero_df = pd.concat(df_list)
    stitched_pysero_df.reset_index(drop=True, inplace=True)
    # remove empty xkappa-biotin spots, round off dilution
    stitched_pysero_df = stitched_pysero_df[(stitched_pysero_df['antigen'] != 'xkappa-biotin') |
                            (stitched_pysero_df['antigen type'] == 'Fiducial')]
    stitched_pysero_df['serum dilution'] = stitched_pysero_df['serum dilution'].round(7)
    stitched_pysero_df.to_csv(os.path.join(output_dir, 'master_report.csv'))
else:
    stitched_pysero_df = pd.read_csv(os.path.join(output_dir, 'master_report.csv'), index_col=0, low_memory=False)

stitched_pysero_df.loc[stitched_pysero_df['antigen']=='xIgG Fc', 'antigen type'] = 'Positive'
#%% functions to compute ROC curves and AUC
# for plate_id in stitched_pysero_df['plate_id'].unique():
# for plate_id in ['plate_8']:
#     slice_cols = ['pipeline', 'serum ID', 'plate_id']
#     slice_keys = [['python'], sera_roc_list, [plate_id]]
#     scn_psr_slice_actions = ['keep', 'drop', 'keep']

norm_antigen = 'xIgG Fc'
# norm_antigen = 'xkappa-biotin'
# norm_antigen = None
norm_group = 'plate'
aggregate = 'mean'
# aggregate = None
# offset_antigen = 'GFP foldon'
offset_antigen = None
# norm_group = 'well'
offset_group = 'well'
antigen_list = ['SARS CoV2 N 50', 'SARS CoV2 RBD 250', 'SARS CoV2 spike 62.5']
# sample_type = 'Orasure'
sample_type = 'Serum'
# slice_cols = ['serum ID', 'antigen type', 'antigen']
# slice_keys = [sera_roc_list, ['Diagnostic'], antigen_list]
# scn_psr_slice_actions = ['drop', 'keep', 'keep']
slice_cols = ['serum ID', 'antigen type']
slice_keys = [sera_roc_list, ['Diagnostic']]
slice_actions = ['drop', 'keep', 'keep']
fpr = 0.05
# ci = 95
ci = None
hue = 'pipeline'
for pipeline in stitched_pysero_df['pipeline'].unique():
# for pipeline in ['nautilus']:
    df_norm = stitched_pysero_df.copy()
    df_norm = slice_df(df_norm, 'keep', 'pipeline', [pipeline])
    df_norm = slice_df(df_norm, 'keep', 'sample type', [sample_type])
    df_norm = normalize_od(df_norm, norm_antigen, group=norm_group)
    df_norm = offset_od(df_norm, offset_antigen, offset_group)
    suffix = '_'.join([pipeline, sample_type])
    if ci is not None:
        suffix = '_'.join([suffix, 'ci'])
    if norm_antigen is not None:
        suffix = '_'.join([suffix, norm_antigen, 'norm_per', norm_group])
    for col, action, key in zip(slice_cols, slice_actions, slice_keys):
        df_norm = slice_df(df_norm, action, col, key)
    roc_df = df_norm.copy()
    if aggregate is not None:
        roc_df = roc_df.groupby(['antigen', 'serum ID', 'well_id', 'plate_id', 'sample type',
                                 'serum type', 'serum dilution', 'pipeline', 'secondary ID',
                                 'secondary dilution'])['OD'].mean()
        roc_df = roc_df.reset_index()
        suffix = '_'.join([suffix, aggregate])
    # roc_df_split = roc_df.copy()
    # roc_df_split[['antigen', 'antigen conc']] = roc_df['antigen'].str.rsplit(n=1, expand=True)
    # roc_plot_grid(roc_df, fig_path, 'ROC')
    # roc_plot_grid(roc_df, fig_path, 'ROC_' + plate_id)
    roc_plot_grid(roc_df, fig_path, '_'.join(['ROC', suffix]), 'png', ci=ci, fpr=fpr, hue=hue)
    id_vars = [x for x in roc_df.columns if x not in ['False positive rate', 'True positive rate']]
    # roc_df = roc_df.melt(id_vars=id_vars,
    #                      var_name='category',
    #                      value_name='rate',
    #                      )
    # thr_plot_grid(roc_df, fig_path, '_'.join(['ROC_thr', suffix]), 'png')
#%%
df_serum = df_norm[['serum ID', 'serum type']].drop_duplicates()
print(len(df_serum))
print((df_serum['serum type']=='negative').sum())
#%% Plot categorical scatter plot for episurvey
for plate_id in stitched_pysero_df['plate_id'].unique():
# for plate_id in ['plate_4']:
    slice_cols = ['pipeline', 'serum ID', 'plate_id']
    slice_keys = [['python'], sera_cat_list, [plate_id]]
    scn_psr_slice_actions = ['keep', 'drop', 'keep']
    antigens = natsorted(stitched_pysero_df['antigen'].unique())
    serum_df = stitched_pysero_df.copy()
    for col, action, key in zip(slice_cols, scn_psr_slice_actions, slice_keys):
        serum_df = slice_df(serum_df, action, col, key)
    assert not serum_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
    # Draw a categorical scatterplot to show each observation
    g = sns.catplot(x="serum type", y="OD", hue="serum type", col_order=antigens, col="antigen",kind="swarm",
                    palette=["r", "c", "y"], data=serum_df, col_wrap=5)
    plt.savefig(os.path.join(fig_path, 'catplot_{}.jpg'.format(plate_id)),
                                  dpi=300, bbox_inches='tight')
    g.set(ylim=(-0.05, 0.4))
    plt.savefig(os.path.join(fig_path, 'catplot_zoom_{}.jpg'.format(plate_id)),
                                  dpi=300, bbox_inches='tight')
    plt.close('all')
#%%
def joint_plot(df_ori,
            x_col,
            y_col,
            hue,
            title,
            output_path,
            output_fname,
            bw='scott',
            n_levels=60,
            xlim=None,
            ylim=None,
            ):

    # g = sns.JointGrid(x=x_col, y=y_col, data=df)
    #                   # xlim=(0, 50), ylim=(0, 8))
    # g = g.plot_joint(sns.kdeplot, cmap="Purples_d")
    # g = g.plot_marginals(sns.kdeplot, color="m", shade=True)
    df = df_ori.dropna(subset=[x_col, y_col])
    diff_df = df[y_col] - df[x_col]
    me = diff_df.mean()
    mae = diff_df.abs().mean()
    # cmap = sns.cubehelix_palette(as_cmap=True, dark=0, light=1, reverse=False)
    cmap = 'Blues'
    fig = plt.figure()
    fig.set_size_inches((9, 9))
    g = sns.JointGrid(x_col, y_col, df,
                xlim=xlim, ylim=ylim)
    hue_vals = []
    for hue_val, hue_df in df.groupby(hue):
        hue_vals.append(hue_val)
        sns.kdeplot(hue_df[x_col], ax=g.ax_marg_x, legend=False, bw=bw)
        sns.kdeplot(hue_df[y_col], ax=g.ax_marg_y, vertical=True, legend=False, bw=bw)
        sns.kdeplot(hue_df[x_col], hue_df[y_col], ax=g.ax_joint,
                     legend=True, gridsize=400, bw=bw, n_levels=n_levels, shade=True, cmap=cmap)
    xfit = np.linspace(xlim[0], xlim[1], 2)
    g.ax_joint.plot(xfit, xfit, linewidth=5, color='k', linestyle='--', alpha=0.5)
    g.ax_joint.text(0.7 * xlim[1], 0.15 * ylim[1], 'Bias={:.3f}'.format(me), fontsize=16)  # add text
    g.ax_joint.text(0.7 * xlim[1], 0.1 * ylim[1], 'Noise={:.3f}'.format(mae), fontsize=16)  # add text
    plt.title(title)
    plt.legend(hue_vals, loc='upper left')
    plt.savefig(os.path.join(output_path, ''.join([output_fname, '.jpg'])),
                dpi=300, bbox_inches='tight')

#%% pivot the dataframe for xy scatter plot
norm_antigen = 'xIgG Fc'
# norm_antigen = None
norm_group = 'plate'
# norm_group = 'well'
# norm_antigen = 'xkappa-biotin'
limit = [0, 2.0]
neg_limit = [0, 0.1]
x_cols = ['plate_3', 'plate_7', 'plate_4']
y_cols = ['plate_9', 'plate_8', 'plate_10']
# for pipeline in stitched_pysero_df['pipeline'].unique():
for pipeline in ['nautilus']:
    df_norm = stitched_pysero_df.copy()
    df_norm = slice_df(df_norm, 'keep', 'pipeline', [pipeline])
    df_norm = normalize_od(df_norm, norm_antigen, group=norm_group)
    suffix = pipeline
    if norm_antigen is not None:
        suffix = '_'.join([pipeline, norm_antigen, 'norm_per', norm_group])

    pysero_df_pivot = pd.pivot_table(df_norm, values='OD',
                                 index=['well_id', 'antigen_row', 'antigen_col', 'serum ID', 'secondary ID', 'secondary dilution',
           'serum type', 'serum dilution', 'antigen', 'antigen type', 'pipeline'],
                                 columns=['plate_id'])
    pysero_df_pivot.reset_index(inplace=True)

    antigen_OD_df = slice_df(pysero_df_pivot, 'keep', 'antigen type', ['Diagnostic'])
    biotin_OD_df = slice_df(pysero_df_pivot, 'keep', 'antigen', ['xkappa-biotin'])
    # biotin_OD_df = slice_df(biotin_OD_df, 'keep', 'antigen type', ['Fiducial'])
    igg_OD_df = slice_df(pysero_df_pivot, 'keep', 'antigen', ['xIgG Fc'])
    antigen_pos_df = slice_df(antigen_OD_df, 'keep', 'serum type', ['positive'])
    antigen_neg_df = slice_df(antigen_OD_df, 'keep', 'serum type', ['negative'])

    for x_col, y_col in zip(x_cols, y_cols):
        # scatter_plot(pysero_df_pivot, x_col, y_col, fig_path, '_'.join(['OD_scatter', x_col, y_col, suffix]), xlim=limit, ylim=limit)
        # scatter_plot(antigen_OD_df, x_col, y_col, 'antigen', fig_path,
        #              '_'.join(['antigen_OD_scatter', x_col, y_col, suffix]), xlim=limit, ylim=limit)
        # scatter_plot(biotin_OD_df, x_col, y_col, 'biotin', fig_path,
        #              '_'.join(['biotin_OD_scatter', x_col, y_col, suffix]), xlim=limit, ylim=limit)
        # scatter_plot(igg_OD_df, x_col, y_col, 'igg', fig_path,
        #              '_'.join(['igg_OD_scatter', x_col, y_col, suffix]), xlim=limit, ylim=limit)
        joint_plot(antigen_OD_df, x_col, y_col, 'antigen type', 'antigen', fig_path,
                   '_'.join(['antigen_OD_joint', x_col, y_col, suffix]), bw='scott', n_levels=60, xlim=limit, ylim=limit)
        joint_plot(antigen_pos_df, x_col, y_col, 'antigen type', 'antigen', fig_path,
                   '_'.join(['antigen_pos_joint', x_col, y_col, suffix]), bw='scott', n_levels=60, xlim=limit, ylim=limit)
        joint_plot(antigen_neg_df, x_col, y_col, 'antigen type', 'antigen', fig_path,
                   '_'.join(['antigen_neg_joint', x_col, y_col, suffix]), bw='scott', n_levels=60, xlim=neg_limit, ylim=neg_limit)
    plt.close('all')

#%% 4PL fit
slice_cols = ['pipeline', 'serum ID']
slice_keys = [['python'], sera_fit_list]
scn_psr_slice_actions = ['keep', 'keep']
serum_df = stitched_pysero_df.copy()
serum_df_fit = stitched_pysero_df.copy()
for col, action, key in zip(slice_cols, scn_psr_slice_actions, slice_keys):
    serum_df = slice_df(serum_df, action, col, key)
    serum_df_fit = slice_df(serum_df_fit, action, col, key)
serum_df_fit = fit2df(serum_df_fit, fourPL, 'python')

#%% plot the ODs and fits
sera_4pl_list = [' '.join([x, 'fit']) for x in sera_fit_list]
markers = 'o'
hue = 'serum ID'
style = 'serum type'
antigens = natsorted(stitched_pysero_df['antigen'].unique())
assert not serum_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
palette = sns.color_palette(n_colors=len(serum_df[hue].unique()))
print('plotting standard curves...')
g = sns.lmplot(x="serum dilution", y="OD", col_order=antigens,
                hue=hue, hue_order=sera_fit_list, col="antigen", ci='sd', palette=palette, markers=markers,
                 data=serum_df, col_wrap=5, fit_reg=False, x_estimator=np.mean)
palette = sns.color_palette(n_colors=len(serum_df_fit[hue].unique()))
for antigen, ax in zip(antigens, g.axes.flat):
    df_fit = serum_df_fit[(serum_df_fit['antigen'] == antigen)]
    sns.lineplot(x="serum dilution", y="OD", hue=hue, hue_order=sera_4pl_list, data=df_fit,
                 style=style, palette=palette,
                 ax=ax, legend=False)
    ax.set(xscale="log")
plt.savefig(os.path.join(fig_path, '{}_fit.jpg'.format(sera_fit_list)),dpi=300, bbox_inches='tight')

for antigen, ax in zip(antigens, g.axes.flat):
    ax.set(ylim=[-0.05, 1.5])
plt.savefig(os.path.join(fig_path, '{}_fit_zoom.jpg'.format(sera_fit_list)),dpi=300, bbox_inches='tight')
#%%
# #%% Generate plots from pysero
#
# # 4 positive sera vs 4 negative sera for control antigens
#
# fig_path = os.path.join(data_folder, 'pysero_fit_OD_8_sera_per_ag')
#
#
# sera_list = ['pos 1', 'pos 2', 'pos 3', 'pos 4', 'neg 1', 'neg 2', 'neg 3', 'neg 4'] # All +ve sera are affected by comets at some antigens, first 4 negative sera are not affected by comets.
# control_ag = ['xkappa-biotin','xIgG Fc', 'GFP foldon']
#
# # selectsera_ctlag = python_df_fix[(python_df_fix['Sera ID'].isin(seralist) & python_df_fix['antigen'].isin(controlag))]
# # selectsera_diagag = python_df_fix[(python_df_fix['Sera ID'].isin(seralist) & ~python_df_fix['antigen'].isin(controlag))]
#
# python_sera_df = pysero_df[pysero_df['serum ID'].isin(sera_list) & (pysero_df['pipeline']==pipeline)]
# y_list = ["OD", "intensity", "background"]
#
# for y in y_list:
#     g = sns.relplot(x="dilution", y=y,
#                         hue="serum ID", style = "type", col="antigen", ci=None,
#                          data=python_sera_df, col_wrap=3, estimator=np.nanmedian, kind='line')
#     g.set(xscale="log")
#     # g.set(xscale="log", ylim=[-0.05, 0.8], xlim=[1E-6, 2E-2])
#     plt.savefig(os.path.join(fig_path, '_'.join(['pysero', y + '.jpg'])), dpi=300, bbox_inches='tight')
#
#
#
# #%% Compare pysero and Scienion aross all serum and antigen combinations
# fig_path = os.path.join(data_folder, 'pysero_plots', 'pysero_metadata_vs_scienion')
# os.makedirs(fig_path, exist_ok=True)
# y_list = ["OD", "intensity", "background"]
# for serum in pysero_df['serum ID'].unique():
#     serum_df = pysero_df[pysero_df['serum ID']==serum]
#     for y in y_list:
#         g = sns.relplot(x="serum dilution", y=y,
#                     hue="pipeline", col="antigen", ci='sd', style='pipeline',
#                      data=serum_df, col_wrap=3, estimator=np.mean, kind='line')
#         g.set(xscale="log", ylim=[-0.05, 1.2])
#         plt.savefig(os.path.join(fig_path, '_'.join([serum, y + '.jpg'])), dpi=300, bbox_inches='tight')
#         plt.close()
#
#

