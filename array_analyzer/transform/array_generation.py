# bchhun, {2020-04-06}

import numpy as np
import array_analyzer.extract.img_processing as img_processing
from array_analyzer.extract.image_parser import generate_props
from array_analyzer.utils.mock_regionprop import MockRegionprop


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


def get_spot_intensity(coords, im_int, background, params, search_range=2):
    """
    Extract signal and background intensity at each spot given the spot coordinate
    with the following steps:
    1. crop image around each grid point
    2. Segment 1 single spot from each image
    3. Get median intensity within the spot mask
    4. If segmentation in 2. returns no mask, use a circular mask with average spot size as the spot mask and do 3.

    :param coords: list or tuple
        [row, col] coordinates of spots
    :param im_int: ndarray
        intensity image of the spots (signals)
    :param background: ndarray
        background image without spots
    :param params: dict
        parameters parsed from the metadata
    :param float search_range: Factor of bounding box size in which to search for
        spots. E.g. 2 searches 2 * 2 * bbox width * bbox height
    :return: dict
        dictionary with format (row index, column index): prop
    """
    # values in mm
    spot_width = params['spot_width']
    pix_size = params['pixel_size']
    n_rows = params['rows']
    n_cols = params['columns']
    # make spot size always odd
    spot_size = 2 * int(0.3 * spot_width / pix_size) + 1
    bbox_width = bbox_height = spot_size

    row_min = np.min(coords[:, 0])
    row_max = np.max(coords[:, 0])
    col_min = np.min(coords[:, 1])
    col_max = np.max(coords[:, 1])

    # scaled max-col, max-row
    row_range = row_max - row_min
    col_range = col_max - col_min
    cent_map = {}
    cent_map_bg = {}
    for label, coord in enumerate(coords):
        # convert the centroid position to an integer that maps to array indices
        grid_row_idx = int(round((n_rows - 1) * ((coord[0] - row_min) / row_range)))
        grid_col_idx = int(round((n_cols - 1) * ((coord[1] - col_min) / col_range)))
        grid_id = (grid_row_idx, grid_col_idx)
        # make bounding boxes larger to account for interpolation errors
        spot_height = int(np.round(search_range * bbox_height))
        spot_width = int(np.round(search_range * bbox_width))
        im_spot_lg, bbox_lg = img_processing.crop_image_at_center(
            im=im_int,
            center=coord,
            height=spot_height,
            width=spot_width,
        )
        bg_spot_lg, _ = img_processing.crop_image_at_center(
            im=background,
            center=coord,
            height=spot_height,
            width=spot_width,
        )
        im_spot, bbox = img_processing.crop_image_at_center(
            im_int,
            coord,
            bbox_height,
            bbox_width,
        )
        bg_spot, _ = img_processing.crop_image_at_center(
            background,
            coord,
            bbox_height,
            bbox_width,
        )
        prop_int = MockRegionprop(
            intensity_image=im_spot,
            centroid=coord,
            label=label,
            bbox=bbox,
        )
        prop_bg = MockRegionprop(
            intensity_image=bg_spot,
            centroid=coord,
            label=label,
            bbox=bbox,
        )
        mask_spot = img_processing.thresh_and_binarize(
            image=im_spot_lg,
            method='bimodal',
        )
        if np.any(mask_spot):
            for prop, im in zip([prop_int, prop_bg], [im_spot_lg, bg_spot_lg]):
                prop_df = generate_props(
                    mask_spot,
                    intensity_image=im,
                    dataframe=True,
                )
                prop.mean_intensity = prop_df.at[0, 'mean_intensity']
                prop.intensity_image = prop_df.at[0, 'intensity_image']
                prop.image = prop_df.at[0, 'image']
                prop.median_intensity = np.median(prop.intensity_image[prop.image])
                prop.centroid = [bbox_lg[0] + prop_df.at[0, 'centroid-0'],
                                 bbox_lg[1] + prop_df.at[0, 'centroid-1']]
                prop.bbox = [bbox_lg[0] + prop_df.at[0, 'bbox-0'],
                             bbox_lg[1] + prop_df.at[0, 'bbox-1'],
                             bbox_lg[0] + prop_df.at[0, 'bbox-2'],
                             bbox_lg[1] + prop_df.at[0, 'bbox-3']
                             ]

        cent_map[grid_id] = prop_int
        cent_map_bg[grid_id] = prop_bg

    return cent_map, cent_map_bg
