import pytest
import xmltodict
import os
import pandas as pd
import cv2 as cv
import numpy as np


@pytest.fixture(scope="session")
def create_good_xml(tmp_path_factory):
    input_dir = tmp_path_factory.mktemp("input_dir")
    output_dir = tmp_path_factory.mktemp("output_dir")

    fiducials = [{'@row': 0,
                  '@col': 0,
                  '@spot_type': 'Reference, Diagnostic'},
                 {'@row': 1,
                  '@col': 0,
                  '@spot_type': 'Reference, Diagnostic'}
                 ]
    spots = [{'@row': 0,
              '@col': 1,
              '@id': 'spot-1-2',
              '@spot_type': 'Diagnostic'},
             {'@row': 1,
              '@col': 2,
              '@id': 'spot-2-3',
              '@spot_type': 'Diagnostic'},
             {'@row': 3,
              '@col': 4,
              '@id': 'spot-4-5',
              '@spot_type': 'Diagnostic'}
             ]
    repl = [{'@row': 0,
             '@col': 1,
             '@id': 'H1 HA',
             'id': ['spot-1-2', 'spot-2-3']},
            {'@row': 0,
             '@col': 1,
             '@id': 'H3 HA',
             'id': ['spot-4-5']}
            ]
    params = {'@rows': 6,
              '@cols': 6,
              '@vspace': 0.4,
              '@hspace': 0.4,
              '@expected_diameter': 0.2,
              '@background_offset': 1,
              '@background_thickness': 1,
              '@max_diameter': 1,
              '@min_diameter': 1,
              }

    doc = {'configuration': {'well_configurations': {'configuration': {'array': {}}}}}
    # set the hardware parameters
    doc['configuration']['well_configurations']['configuration']['array']['layout'] = params

    # set the fiducials
    doc['configuration']['well_configurations']['configuration']['array']['layout']['marker'] = fiducials

    # set the spot IDs
    doc['configuration']['well_configurations']['configuration']['array']['spots'] = {}
    doc['configuration']['well_configurations']['configuration']['array']['spots']['spot'] = spots

    # set the number of replicates
    doc['configuration']['well_configurations']['configuration']['array']['spots']['multiplet'] = repl

    with open(os.path.join(str(input_dir), 'temp.xml'), 'w', encoding='utf-8') as temp_xml:
        temp_xml.write(xmltodict.unparse(doc))

    return input_dir, output_dir


@pytest.fixture(scope="session")
def create_good_xlsx(tmp_path_factory):
    input_dir = tmp_path_factory.mktemp("input_dir")
    output_dir = tmp_path_factory.mktemp("pysero_output_dir")

    # make a dummy worksheet with realistic parameters
    params_worksheet = {'': '',
                        'rows': '6',
                        'columns': '6',
                        'v_pitch': '0.4',
                        'h_pitch': '0.45',
                        'spot_width': '0.2',
                        'pixel_size': '0.0049'
                        }
    keys = dict()
    vals = dict()
    for idx, value in enumerate(params_worksheet.keys()):
        keys[idx] = value
    for idx, value in enumerate(params_worksheet.values()):
        vals[idx] = value
    p = pd.Series(keys, name="Parameter")
    v = pd.Series(vals, name="Value")
    params_df = pd.DataFrame([p, v]).T

    # make a dummy antigen array layout with realistic fiducials
    fiducials = {0: {0: 'Fiducial', 1: '', 2: '', 3: '', 4: '', 5: 'Fiducial'},
                 1: {0: 'Fiducial', 1: '', 2: '', 3: '', 4: '', 5: ''},
                 2: {0: 'Positive Control', 1: '', 2: '', 3: '', 4: '', 5: 'Negative Control'},
                 3: {0: 'Positive Control', 1: '', 2: '', 3: '', 4: '', 5: 'Negative Control'},
                 4: {0: 'Positive Control', 1: '', 2: '', 3: '', 4: '', 5: 'Negative Control'},
                 5: {0: 'Fiducial', 1: '', 2: '', 3: '', 4: '', 5: 'Fiducial'}
                 }
    fiducials_df = pd.DataFrame(fiducials).T

    # make a dummy antigen array layout with realistic antigens
    antigens = {0: {0: 'xkappa-biotin', 1: 'Flu vaccine 2018-2019', 2: 'Flu vaccine 2018-2019',
                    3: 'Flu vaccine 2018-2019', 4: 'Flu vaccine 2018-2019', 5: 'xkappa-biotin'},
                1: {0: 'xkappa-biotin', 1: 'H1 HA', 2: 'H1 HA', 3: 'H1 HA', 4: 'H1 HA', 5: ''},
                2: {0: 'xlgG Fc', 1: 'H3 HA', 2: 'H3 HA', 3: 'H3 HA', 4: 'H3 HA', 5: 'GFP foldon'},
                3: {0: 'xlgG Fc', 1: 'H7 HA', 2: 'H7 HA', 3: 'H7 HA', 4: 'H7 HA', 5: 'GFP foldon'},
                4: {0: 'xlgG Fc', 1: 'HA FluB I', 2: 'HA FluB I', 3: 'HA FluB I', 4: 'HA FluB I', 5: 'GFP foldon'},
                5: {0: 'xkappa-biotin', 1: 'HA FluB II', 2: 'HA FluB II', 3: 'HA FluB II', 4: 'HA FluB II',
                    5: 'xkappa-biotin'}
                }
    antigens_df = pd.DataFrame(antigens).T

    # Add names of wells to rerun
    rerun_df = pd.DataFrame(data={'well_name': ['A3', 'B7']})

    # writing a two-worksheet excel file
    writer = pd.ExcelWriter(
        os.path.join(str(input_dir), 'pysero_output_data_metadata.xlsx'),
        index=False,
        engine='openpyxl',
    )
    params_df.to_excel(writer, sheet_name='imaging_and_array_parameters')
    fiducials_df.to_excel(writer, sheet_name='antigen_type')
    antigens_df.to_excel(writer, sheet_name='antigen_array')
    rerun_df.to_excel(writer, sheet_name='rerun_wells')
    writer.save()
    writer.close()

    return input_dir, output_dir


@pytest.fixture(scope="session")
def image_dir(tmpdir_factory):
    """
    Creates one directory with a number of well images named A1.png, ...

    :param tmpdir_factory: PyTest factory for temporary directories
    :return input_dir: Temporary directory with well images
    """
    im = np.zeros((5, 10), dtype=np.uint8)
    im[2:3, 3:8] = 128
    input_dir = tmpdir_factory.mktemp("input_dir")
    # Save a couple of images
    cv.imwrite(os.path.join(input_dir, 'A1.png'), im)
    cv.imwrite(os.path.join(input_dir, 'A2.png'), im)
    cv.imwrite(os.path.join(input_dir, 'B11.png'), im + 50)
    cv.imwrite(os.path.join(input_dir, 'B12.png'), im + 100)
    return input_dir


@pytest.fixture(scope="session")
def micromanager_dir(tmpdir_factory):
    """
    Creates a directory where well images are in one subdirectory each.
    Well names are encoded in the subdirectory name.

    :param tmpdir_factory: PyTest factory for temporary directories
    :return input_dir: Directory with well images in subdirectories
    """
    im = np.zeros((5, 10), dtype=np.uint8)
    im[2:3, 3:8] = 128
    input_dir = tmpdir_factory.mktemp("input_dir")
    well_names = ['A1', 'A2', 'B11', 'B12']
    # Save a couple of images
    for well in well_names:
        sub_dir = input_dir / well + "-what-ever"
        sub_dir.mkdir()
        cv.imwrite(os.path.join(sub_dir, 'micromanager_name.tif'), im)
    return input_dir
