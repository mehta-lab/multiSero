# bchhun, {2020-04-06}

import numpy as np
from ..extract.img_processing import crop_image_at_center, thresh_and_binarize
from ..extract.image_parser import generate_props
from ..utils.mock_regionprop import MockRegionprop


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


def get_spot_intensity(coords, im_int, background, params):
    """
    Extract signal and background intensity at each spot given the spot coordinate
    with the followling steps:
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

    y_min = np.min(coords[:, 0])
    y_max = np.max(coords[:, 0])
    x_min = np.min(coords[:, 1])
    x_max = np.max(coords[:, 1])

    # scaled max-x, max-y
    y_range = y_max - y_min
    x_range = x_max - x_min
    cent_map = {}
    cent_map_bg = {}
    for coord in coords:
        cen_y, cen_x = coord
        # convert the centroid position to an integer that maps to array indices
        grid_y_idx = int(round((n_rows - 1) * ((cen_y - y_min) / y_range)))
        grid_x_idx = int(round((n_cols - 1) * ((cen_x - x_min) / x_range)))
        grid_id = (grid_y_idx, grid_x_idx)
        label = 0
        # make bounding boxes larger to account for interpolation errors
        im_1spot_lg, bbox_lg = crop_image_at_center(im_int,
                                                    coord,
                                                    4 * bbox_height,
                                                    4 * bbox_width,
                                                    )

        bg_1spot_lg, _ = crop_image_at_center(background,
                                              coord,
                                              4 * bbox_height,
                                              4 * bbox_width,
                                              )
        im_1spot, bbox = crop_image_at_center(im_int,
                                              coord,
                                              bbox_height,
                                              bbox_width,
                                              )

        bg_1spot, _ = crop_image_at_center(background,
                                           coord,
                                           bbox_height,
                                           bbox_width,
                                           )

        prop_int = MockRegionprop(intensity_image=im_1spot,
                              centroid=coord,
                              label=label,
                              bbox=bbox)

        prop_bg = MockRegionprop(intensity_image=bg_1spot,
                                  centroid=coord,
                                  label=label,
                                 bbox=bbox)

        mask_1spot = thresh_and_binarize(im_1spot_lg, method='bright_spots', invert=True, min_size=10, thr_percent=92)
        if np.any(mask_1spot):
            for prop, im in zip([prop_int, prop_bg], [im_1spot_lg, bg_1spot_lg]):
                prop_df = generate_props(mask_1spot, intensity_image_=im, dataframe=True)
                # select the object with max area
                if len(prop_df.index) > 1:
                    prop_df = prop_df.loc[[prop_df['area'].idxmax()]]
                    prop_df.reset_index(drop=True, inplace=True)

                if len(prop_df.index) == 1:
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
        label += 1
    return cent_map, cent_map_bg
