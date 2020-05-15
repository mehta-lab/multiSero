import pytest

import array_analyzer.utils.io_utils as io_utils
import array_analyzer.extract.image_parser as image_parser
import array_analyzer.extract.img_processing as img_processing
import warnings
import numpy as np


def test_scienion_clean_well_crop(scienion_clean):
    """
    test the function image_parser.find_well_border on "good" data

    :param scienion_clean:
    :return:
    """
    gdd_path = scienion_clean

    raw = gdd_path+'/Clean_March25_C4.png'
    crop = gdd_path+'/march25_c4_im_well.npy'

    image = io_utils.read_gray_im(raw.strpath)
    target = np.load(crop.strpath)

    try:
        well_center, well_radi, well_mask = image_parser.find_well_border(
            image,
            detmethod='region',
            segmethod='otsu',
        )
        im_well, _ = img_processing.crop_image_at_center(
            image,
            well_center,
            2 * well_radi,
            2 * well_radi,
        )
    except IndexError:
        warnings.warn("Couldn't find well in {}".format(raw))
        im_well = image

    assert(im_well.all() == target.all())
    # todo: add more assertions


def test_scienion_clean_spot_coords(scienion_clean):
    """
    test
    :param scienion_clean:
    :return:
    """
    gdd_path = scienion_clean

    im_well_path = gdd_path+'/march25_c4_im_well.npy'
    target_spot_coords_path = gdd_path+'/march25_c4_spot_coords.npy'

    im_well = np.load(im_well_path.strpath)
    target = np.load(target_spot_coords_path.strpath)
    params = {
        'rows': 6,
        'columns': 6,
        'v_pitch': 0.4,
        'h_pitch': 0.4,
        'spot_width': 0.2,
        'bg_offset': None,
        'bg_thickness': None,
        'max_diam': None,
        'min_diam': None,
        'pixel_size_scienion': 0.0049,
        'pixel_size_octopi': 0.00185,
        'pixel_size': 0.0049
    }

    spot_detector = img_processing.SpotDetector(
        imaging_params=params,
    )

    spot_coords = spot_detector.get_spot_coords(im_well)

    assert(target.all() == spot_coords.all())

