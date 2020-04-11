# bchhun, {2020-03-26}

import skimage as si
import skimage.io
import os
import numpy as np
import matplotlib.pyplot as plt


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


def plot_spot_assignment(od_well, i_well, bg_well,
                         im_crop, props_by_loc,
                         bgprops_by_loc, image_name,
                         output_name, params):
    plt.imshow(im_crop,cmap='gray')
    plt.colorbar()
    for r in np.arange(params['rows']):
        for c in np.arange(params['columns']):
            try:
                ceny, cenx = props_by_loc[(r, c)].centroid
            except:
                spot_text = '(' + str(r) + ',' + str(c) + ')'
                print(spot_text + 'not found')
            else:
                cenybg, cenxbg = bgprops_by_loc[(r, c)].centroid
                plt.plot(cenx, ceny, 'm+',ms=10)
                plt.plot(cenxbg, cenybg, 'gx',ms=10)
                spot_text = '(' + str(r) + ',' + str(c) + ')'
                plt.text(cenx, ceny - 5, spot_text, va='bottom', ha='center', color='w')
                plt.text(0, 0, image_name[:-4] + ',spot count=' + str(len(props_by_loc)))
    figcentroid = plt.gcf()
    centroids_debug = output_name + '_overlayCentroids.png'
    figcentroid.savefig(centroids_debug, bbox_inches='tight')
    plt.close(figcentroid)

    plt.figure(figsize=(6, 1.5))
    plt.subplot(131)
    plt.imshow(i_well, cmap='gray')
    plt.colorbar()
    plt.title('intensity')

    plt.subplot(132)
    plt.imshow(bg_well, cmap='gray')
    plt.colorbar()
    plt.title('background')

    plt.subplot(133)
    plt.imshow(od_well, cmap='gray')
    plt.colorbar()
    plt.title('OD')

    figOD = plt.gcf()
    od_debug = output_name + '_od.png'
    figOD.savefig(od_debug)
    plt.close(figOD)
