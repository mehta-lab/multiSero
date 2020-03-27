# bchhun, {2020-03-26}

import skimage as si
import skimage.io
import os
import numpy as np


def save_all_wells(region_props_array, spot_ids_, output_folder, well_name):

    for row in range(region_props_array.shape[0]):
        for col in range(region_props_array.shape[1]):

            cell = spot_ids_[row][col]
            if cell == '':
                continue

            prop = region_props_array[row][col]
            if prop is not None:
                si.io.imsave(output_folder + os.sep + well_name + f"_{cell}.png",
                             (255 * prop.intensity_image).astype('uint8'))
            else:
                si.io.imsave(output_folder + os.sep + well_name + f"_{cell}.png",
                             (255 * np.ones((32, 32)).astype('uint8')))


def assign_region(target_, props_, intensity_image_=None):
    """
    put underlying intensity image from props_ into the target_ array
    assumes props_ image exists and will fit into target_

    :param target_:
    :param props_:
    :param intensity_image_:
    :return:
    """

    min_row, min_col, max_row, max_col = props_.bbox
    if intensity_image_ is None:
        target_[min_row:max_row, min_col:max_col] = props_.intensity_image
    else:
        target_[min_row:max_row, min_col:max_col] = intensity_image_[min_row:max_row, min_col:max_col]

    return target_


def create_composite_spots(target_array_, region_props_array_, intensity_image_=None):
    """
    insert all intensity images for each region prop in "region_props_array_" into the "target_array_"

    :param target_array_:
    :param region_props_array_:
    :param intensity_image_:
    :return:
    """
    for row in range(region_props_array_.shape[0]):
        for col in range(region_props_array_.shape[1]):
            props = region_props_array_[row, col]
            if props is not None:
                if intensity_image_ is None:
                    # print(f"\t shape = {props.intensity_image.shape}")
                    target_array_ = assign_region(target_array_, props)
                else:
                    target_array_ = assign_region(target_array_, props, intensity_image_)
    return target_array_


def save_composite_spots(source_array_, region_props_array_, output_folder, well_name, from_source=False):
    """

    :param source_array_: np.ndarray
        image representing original image data
    :param region_props_array_: np.ndarray
        region props describing segmented spots from image data
    :param output_folder: str
    :param well_name: str
        well name of format "A1, A2 ... C2, C3"
    :param from_source: bool
        True : images are extracted from source array
        False : images are pulled from regionprops.intensity_image
    :return:
    """

    t = np.mean(source_array_) * np.ones(shape=(source_array_.shape[0], source_array_.shape[1]))

    if from_source:
        t = create_composite_spots(t, region_props_array_, source_array_)
        si.io.imsave(output_folder + os.sep + well_name + f"_composite_spots_img.png",
                     (255 * t).astype('uint8'))
    else:
        t = create_composite_spots(t, region_props_array_)
        si.io.imsave(output_folder + os.sep + well_name + f"_composite_spots_prop.png",
                     (255 * t).astype('uint8'))




