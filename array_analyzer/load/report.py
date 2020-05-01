import array_analyzer.extract.constants as c


def write_well_arrays(data, well_name, array_type):
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

    if well_name in c.IMAGE_TO_WELL:
        (row, col) = c.IMAGE_TO_WELL[well_name]
    else:
        raise AttributeError(f"well name {well_name} is not recognized")

    if array_type == 'od':
        c.WELL_OD_ARRAY[row-1, col-1] = data
    if array_type == 'int':
        c.WELL_INT_ARRAY[row-1, col-1] = data
    if array_type == 'bg':
        c.WELL_BG_ARRAY[row-1, col-1] = data

