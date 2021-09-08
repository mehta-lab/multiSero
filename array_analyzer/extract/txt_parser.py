# bchhun, {2020-03-22}
import csv
import natsort
import numpy as np
import os
import xmltodict
from xml.parsers.expat import ExpatError
import xml.etree.ElementTree as ET
import pandas as pd
import math

import array_analyzer.extract.constants as constants

"""
functions like "create_<extension>_dict" parse files of <extension> and return:
fiduc: list of dict describing fiducial positions
spots: list of dict, other spot info
repl: list of dict describing 'replicates' AKA antigens
params: dict containing hardware and array parameters


functions like "populate_array_<type>" take <type> from above (like fiduc, spots, repl, param) and 
    populate np.ndarrays indices correspond to array positions
    
The values of the arrays depend on the function call
- populate_array_id : Cell id like "spot-6-2", "spot-5-5-" etc..
- populate_array_spots_type : Type like "Diagnostic", "Positive Control"
- populate_array_antigen : Antigen


*** NOTE ***
populating antigens is more complicated for .xml parsing than for .csv or .xlsx
    .xml files have keys:"antigen", values: multiple "spot_ID"
    .csv or .xlsx can map (col, row) directly to "antigen"
"""


def create_xml_dict(path_):
    """
    receives an .xml file generated by the Scienion sciReader software
    and returns dictionaries containing info.

    :param str path_: Full path of xml file
    :return dict fiduc: Fiducials and control info
    :return dict spots: Spot info
    :return dict repl: Replicate info
    :return dict params: Additional parameters
    """
    try:
        with open(path_) as fd:
            doc = xmltodict.parse(fd.read())
    except ExpatError:
        tree = ET.parse(path_)
        xml_data = tree.getroot()
        # here you can change the encoding type to be able to set it to the one you need
        xmlstr = ET.tostring(xml_data, encoding='utf-8', method='xml')
        doc = xmltodict.parse(xmlstr)

    try:
        # layout of array
        layout = doc['configuration']['well_configurations']['configuration']['array']['layout']

        # fiducials
        fiduc = layout['marker']

        # spot IDs
        spots = doc['configuration']['well_configurations']['configuration']['array']['spots']['spot']

        # replicates
        repl = doc['configuration']['well_configurations']['configuration']['array']['spots']['multiplet']

        array_params = dict()
        array_params['rows'] = int(layout['@rows'])
        array_params['columns'] = int(layout['@cols'])
        array_params['v_pitch'] = float(layout['@vspace'])
        array_params['h_pitch'] = float(layout['@hspace'])
        array_params['spot_width'] = float(layout['@expected_diameter'])
        array_params['bg_offset'] = float(layout['@background_offset'])
        array_params['bg_thickness'] = float(layout['@background_thickness'])
        array_params['max_diam'] = float(layout['@max_diameter'])
        array_params['min_diam'] = float(layout['@min_diameter'])

    except Exception as ex:
        raise AttributeError(f"exception while parsing .xml : {ex}")

    return fiduc, spots, repl, array_params


def create_csv_dict(path_):
    """
    Looks for three .csv files:
        "array_parameters.csv" contains array printing parameters as well as hardware parameters
        "array_format_type.csv" contains fiducial and control names, locations
        "array_format_antigen.csv" contains specific antigen names for only diagnostic spots
    Then, parses the .csvs and creates the four dictionaries:
    :param path_: list
        of strings to the .csv paths
    :return:
    """
    # assign names to each of the types == params, spot types, antigens
    fiduc = list()
    csv_antigens = list()
    array_params = dict()

    for meta_csv in path_:
        with open(meta_csv, newline='') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                # header row
                if "Parameter" in str(row[0]) or '' in str(row[0]):
                    continue

                # parse params
                elif "array_format_parameters" in meta_csv:
                    array_params[row[0]] = row[1]

                # parse fiducials
                elif "array_format_type" in meta_csv:
                    for col, value in enumerate(row[1:]):
                        pos = {'@row': str(row[0]),
                               '@col': str(col - 1),
                               '@spot_type': str(value)}
                        fiduc.append(pos)

                # parse antigens
                elif "array_format_antigen" in meta_csv:
                    for col, value in enumerate(row[1:]):
                        pos = {'@row': str(row[0]),
                               '@col': str(col - 1),
                               '@antigen': str(value)}
                        csv_antigens.append(pos)

    return fiduc, None, csv_antigens, array_params


def create_xlsx_dict(xlsx):
    """
    extracts fiducial, antigen, and array parameter metadata from .xlsx sheets
    then populates dictionaries or lists with appropriate information
    The output dictionaries and lists conform the .xml-style parsing.  This is for consistency

    :param dict xlsx: Opened xlsx sheets
    :return list fiduc: Fiducials and control info
    :return list spots: None.  spots IDs not needed for .xlsx
    :return list repl: Replicate (antigen)
    :return dict params: Additional parameters about hardware and array
    """
    fiduc = list()
    xlsx_antigens = list()
    array_params = dict()

    # populate array parameters
    for idx, value in enumerate(xlsx['imaging_and_array_parameters']['Parameter']):
        array_params[value] = xlsx['imaging_and_array_parameters']['Value'][idx]

    # Unless specified, run analysis with fiducials only
    # Otherwise run with all non-negative spots
    fiducials_only = True
    if 'fiducials_only' in array_params:
        if array_params['fiducials_only'] != 1:
            fiducials_only = False

    # populate fiduc list
    for col in xlsx['antigen_type'].keys()[1:]:
        for row, value in enumerate(xlsx['antigen_type'][col]):
            if type(value) is float:
                if math.isnan(value):
                    continue
            else:
                if not fiducials_only and "Negative" not in value:
                    pos = {'@row': row,
                           '@col': col,
                           '@spot_type': "Fiducial"}
                    fiduc.append(pos)
                elif "Fiducial" in value or "xkappa-biotin" in value or "Fiducial, Diagnostic" in value:
                    pos = {'@row': row,
                           '@col': col,
                           '@spot_type': "Fiducial"}
                    fiduc.append(pos)

    # find and populate antigen list
    for col in xlsx['antigen_array'].keys()[1:]:
        for row, value in enumerate(xlsx['antigen_array'][col]):
            if type(value) is float:
                if math.isnan(value):
                    continue
            else:
                pos = {'@row': row,
                       '@col': col,
                       '@antigen': str(value)}
                xlsx_antigens.append(pos)

    return fiduc, xlsx_antigens, array_params


def create_xlsx_array(path_):
    """
    (unfinished attempt to convert arrays in xlsx to np.ndarrays directly, using pandas)

    :param path_:
    :return:
    """

    array_params = dict()

    # populate array parameters
    params = pd.read_excel(path_, sheet_name='imaging_and_array_parameters')
    for idx, value in enumerate(params['Parameter']):
        array_params[value] = params['Value'][idx]

    # populate fiducials array
    fiduc_df = pd.read_excel(path_, sheet_name='antigen_type')
    fiduc = fiduc_df.to_numpy(dtype='U100')

    # populate xlsx
    xlsx_antigens_df = pd.read_excel(path_, sheet_name='antigen_array')
    xlsx_antigens = xlsx_antigens_df.to_numpy(dtype='U100')

    return fiduc, None, xlsx_antigens, array_params


def create_array(rows_, cols_, dtype='U100'):
    """
    creates an empty numpy array whose elements are long strings
    :param rows_: int
        provided by params
    :param cols_: int
        provided by params
    :param dtype: str
        type of element in this np array
    :return: np.ndarray
    """

    xml_numpy_array = np.empty(shape=(rows_, cols_), dtype=np.dtype(dtype))

    return xml_numpy_array


def populate_array_id(arr, spots):
    """
    receives an empty array
    populates the array with values corresponding to spot ID:
        ID like : "spot-1-2", "spot-2-4", "spot-3-3-", etc...

    :param arr: np.ndarray
        numpy array generated from "create_array"
    :param spots: dict
        dict from "create_xml_dict"
    :return: np.ndarray
        populated array
    """
    for spot in spots:
        r = int(spot['@row'])
        c = int(spot['@col'])
        ID = spot['@id']

        arr[r, c] = ID

    return arr


def populate_array_spots_type(arr, spots, fiduc):
    """
    receives an empty array
    populates the array with values corresponding to "spot_type":
        spot_type like : "Diagnostic", "PositiveControl", "NegativeControl"

    :param arr: np.ndarray
        numpy array generated from "create_array"
    :param spots: dict
        list of dict from "create_xml_dict"
    :param fiduc: list
        list of dict from "create_xml_dict"
    :return: np.ndarray
        populated array
    """

    for spot in spots:
        r = int(spot['@row'])
        c = int(spot['@col'])
        v = spot['@spot_type']

        arr[r, c] = v

    for f in fiduc:
        r = int(f['@row'])
        c = int(f['@col'])
        v = f['@spot_type']

        arr[r, c] = v

    return arr


def populate_array_fiduc(arr, fiduc):
    """
    assigns only fiducials to the positions of the array
    :param arr: np.ndarray
        target array
    :param fiduc: list
        list of dict output of "create_xml_dict"
    :return: np.ndarray
        modified target array
    """

    for f in fiduc:
        r = int(f['@row'])
        c = int(f['@col'])
        v = f['@spot_type']

        arr[r, c] = v

    return arr


def populate_array_antigen_xml(arr, id_arr_, repl):
    """
    populates an array with the antigen
    scans through "replicate" in the .xml and assigns all spots the appropriate antigen

    :param arr: np.ndarray
        numpy array generated from "create_array"
    :param id_arr_: np.ndarray
        array from populate_array_id
    :param repl: dict
        dict from "create_xml_dict"
    :return: np.ndarray
        populated array
    """
    for rep in repl:
        antigen = rep['@id']
        all_spots = rep['id']  # list of IDs
        for spot in all_spots:
            arr[np.where(id_arr_ == spot)] = antigen

    return arr


def populate_array_antigen(arr, csv_antigens_):
    """
    populates an array with antigen
        used for metadata obtained through .csv files
    :param arr:
    :param csv_antigens_: list
        list of dictionaries whose keys define array coordinates and values
    :return:
    """
    for antigen in csv_antigens_:
        r = antigen['@row']
        c = int(antigen['@col'])
        v = antigen['@antigen']
        arr[r, c] = v

    return arr


def rerun_xl_od(well_names, well_xlsx_path, rerun_names, xlsx_writer):
    """
    Load stats_per_well excel file and copy over existing well sheets
    before rerunning some of the wells.

    :param list well_names: Well names (e.g. ['B12', 'C2'])
    :param str well_xlsx_path: Full path to well stats xlsx sheet
    :param list rerun_names: Names of wells to be rerun
    :param pd.ExcelWriter xlsx_writer: Pandas excel writer
    """
    rerun_set = set(rerun_names)
    assert rerun_set.issubset(well_names), \
        "All rerun wells can't be found in input directory"
    assert os.path.isfile(well_xlsx_path),\
        "Can't find stats_per_well excel: {}".format(well_xlsx_path)
    ordered_dict = pd.read_excel(well_xlsx_path, sheet_name=None)
    written_wells = list(ordered_dict.keys())
    written_wells.remove('antigens')
    # Find the difference between the sets
    existing_wells = natsort.natsorted(
        list(set(written_wells) - rerun_set),
    )
    # Write existing wells to well stats
    for well_name in existing_wells:
        well_df = pd.DataFrame(ordered_dict[well_name])
        well_df.to_excel(xlsx_writer, sheet_name=well_name)
