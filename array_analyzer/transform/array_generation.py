# bchhun, {2020-04-06}

import numpy as np


def build_centroid_binary_blocks(cent_list, image_, params_, return_type='region'):
    """
    places blocks at the center of each centroid supplied by "cent_list".  Returns this mask

    :param cent_list: list
        of centroids (x, y)
    :param image_: ndarray
        underlying image whose shape is output
    :param params_: dict
        parameters parsed from .txt metadata
    :param return_type: str
        to return the region mask or the product of mask*image
    :return: ndarray
    """

    # build_block_array uses values in the params_ to motivate array dimensions, spacings

    # values in mm
    v_pitch = params_['v_pitch']
    h_pitch = params_['h_pitch']
    spot_width = params_['spot_width']
    PIX_SIZE = params_['pixel_size_scienion']

    # values in pixels
    v_pix = v_pitch/PIX_SIZE
    h_pix = h_pitch/PIX_SIZE
    spot_pix = spot_width/PIX_SIZE

    # make the box 1.3x the size of the spot, unless it will cause overlap
    side = int(1.3*spot_pix if 1.3*spot_pix < v_pix-1 and 1.3*spot_pix < h_pix-1 else spot_pix)

    target = np.zeros(image_.shape)
    for centroid in cent_list:
        (center_y, center_x) = centroid

        # check that the blank fits within the bounds of the target array
        x_min = int(center_x - (side/2)) if int(center_x - (side/2)) > 0 else 0
        x_max = int(center_x + (side/2)) if int(center_x + (side/2)) < target.shape[0] else target.shape[0]
        y_min = int(center_y - (side/2)) if int(center_y - (side/2)) > 0 else 0
        y_max = int(center_y + (side/2)) if int(center_y + (side/2)) < target.shape[1] else target.shape[1]

        blank = np.ones((x_max - x_min, y_max - y_min))

        target[x_min:x_max, y_min:y_max] = blank

    if return_type == 'product':
        return target*image_
    elif return_type == 'region':
        return target
