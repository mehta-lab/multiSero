import pytest
import xmltodict
import os
import pandas as pd
import cv2 as cv
import numpy as np

from google_drive_downloader import GoogleDriveDownloader as gdd
import shutil


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture(scope="session")
def create_good_xml(tmp_path_factory):
    input_dir = tmp_path_factory.mktemp("input_dir")

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

    return str(input_dir)


@pytest.fixture(scope="session")
def create_good_xlsx(tmp_path_factory):
    input_dir = tmp_path_factory.mktemp("input_dir")

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

    # writing a two-worksheet excel file
    writer = pd.ExcelWriter(os.path.join(str(input_dir), 'pysero_output_data_metadata.xlsx'),
                            index=False,
                            engine='openpyxl')
    params_df.to_excel(writer, sheet_name='imaging_and_array_parameters')
    fiducials_df.to_excel(writer, sheet_name='antigen_type')
    antigens_df.to_excel(writer, sheet_name='antigen_array')
    writer.save()
    writer.close()

    return input_dir


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


# ================  fixtures for registration_workflow unit tests =========

@pytest.fixture(scope="session")
def scienion_clean(tmpdir_factory):
    """
    downloads and provides "clean" data acquired on scienion
    zip includes:
    - (raw).png
    - im_well.npy (cropped well)
    - grid_coords.npy (create_reference_grid)
    - spot_coords.npy (get_spot_coords)
    - reg_coords.npy (output of full registration)
    :param tmpdir_factory:
    :return:
    """

    input_dir = tmpdir_factory.mktemp("scienion_clean")

    clean_march_25_c4_zip = '13D2EHjg--dLTx3OkyqJCRacoEJHTg-9r'
    target_dir = input_dir + "/srczip.zip"

    gdd.download_file_from_google_drive(file_id=clean_march_25_c4_zip,
                                        dest_path=target_dir,
                                        unzip=True,
                                        showsize=True,
                                        overwrite=True)
    return input_dir


@pytest.fixture(scope="session")
def scienion_dir(tmpdir_factory):
    """
    Creates a directory with scienion well images.

    :param tmpdir_factory: PyTest factory for temporary directories
    :return input_dir: Directory with well images in subdirectories
    :return output_dir: Directory for writing output
    """
    input_dir = tmpdir_factory.mktemp("input_dir")
    output_dir = input_dir / "integration_dir"
    output_dir.mkdir()

    # Get images and metadata
    # clean, comet, missing fiducial, particles and metadata
    image_ids = [
        '116oiRJSX4FNWWxtqxuu-YsfAC10iBF7a',
        '1-VCwS5PvJnGID-3wQEiPVDoGH4LOkF0P',
        '11DlsQjJJurDQ960oEc3mgZnjGO1vyLL2',
        '1-XjOog2pzI8U9JyOV7orn1fRuYq9cNcS',
    ]
    meta_id = '10x7LIWJa_IBTGkKC8ZUq_3JiqBGd3uui'
    # Download images
    well_names = ['A1', 'A2', 'B11', 'B12']
    for i, gdrive_id in enumerate(image_ids):
        target_path = os.path.join(input_dir, well_names[i] + '.png')
        gdd.download_file_from_google_drive(
            file_id=gdrive_id,
            dest_path=target_path,
        )
    # Download metadata
    target_path = os.path.join(input_dir, 'pysero_output_data_metadata.xlsx')
    gdd.download_file_from_google_drive(
        file_id=meta_id,
        dest_path=target_path,
    )
    return input_dir, output_dir
