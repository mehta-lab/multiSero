# bchhun, {2020-03-26}

import cv2 as cv
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


def save_composite_spots(source_array_,
                         region_props_array_,
                         output_name,
                         from_source=False):
    """

    :param source_array_: np.ndarray
        image representing original image data
    :param region_props_array_: np.ndarray
        region props describing segmented spots from image data
    :param str output_name: Path plus well name, no extension
        well name of format "A1, A2 ... C2, C3"
    :param from_source: bool
        True : images are extracted from source array
        False : images are pulled from regionprops.intensity_image
    :return:
    """
    t = np.mean(source_array_) * np.ones(source_array_.shape)

    if from_source:
        t = create_composite_spots(t, region_props_array_, source_array_)
        cv.imwrite(
            output_name + "_composite_spots_img.png",
            (255 * t).astype('uint8'),
        )
    else:
        t = create_composite_spots(t, region_props_array_)
        cv.imwrite(
            output_name + "_composite_spots_prop.png",
            (255 * t).astype('uint8'),
        )


def plot_centroid_overlay(im_crop,
                          params,
                          props_by_loc,
                          bgprops_by_loc,
                          output_name):

    plt.imshow(im_crop, cmap='gray')
    plt.colorbar()
    im_name = os.path.basename(output_name)
    for r in np.arange(params['rows']):
        for c in np.arange(params['columns']):
            try:
                ceny, cenx = props_by_loc[(r, c)].centroid
            except:
                spot_text = '(' + str(r) + ',' + str(c) + ')'
                print(spot_text + 'not found')
            else:
                cenybg, cenxbg = bgprops_by_loc[(r, c)].centroid
                plt.plot(cenx, ceny, 'm+', ms=10)
                plt.plot(cenxbg, cenybg, 'gx', ms=10)
                spot_text = '(' + str(r) + ',' + str(c) + ')'
                plt.text(cenx, ceny - 5, spot_text, va='bottom', ha='center', color='w')
                plt.text(0, 0, im_name + ',spot count=' + str(len(props_by_loc)))

    figcentroid = plt.gcf()
    centroids_debug = output_name + '_overlay_centroids.png'
    figcentroid.savefig(centroids_debug, bbox_inches='tight')
    plt.close(figcentroid)


def plot_od(od_well,
            i_well,
            bg_well,
            output_name):

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


def plot_background_overlay(im, im_background, output_name):
    """
    Writes color image with background overlaid.

    :param np.array im: 2D grayscale image
    :param np.array im_background: 2D grayscale image
    :param str output_name: Path and image name minus extension
    """
    im_stack = np.stack([im_background, im, im_background], axis=2)
    skimage.io.imsave(
        output_name + "_crop_bg_overlay.png",
        (255 * im_stack).astype('uint8'),
    )


def plot_registration(image,
                      spot_coords,
                      grid_coords,
                      reg_coords,
                      output_name,
                      max_intensity=255,
                      margin=100):
    """
    Plots all detected spots, initial fiducial coordinates and registered grid.

    :param np.array image: Input image
    :param np.array spot_coords: Detected spot coordinates (nbr spots x rows, cols)
    :param np.array grid_coords: Initial estimate of fiducial coordinates
    :param np.array reg_coords: Registered coordinates
    :param str output_name: Path + well name, _registration.png will be added
    :param int max_intensity: Maximum image intensity (expecting uint8 or 16)
    :param int margin: Margin around spots to crop image before plotting
    """
    im = (image / max_intensity * 255).astype(np.uint8)

    all_coords = np.vstack([spot_coords, grid_coords, reg_coords])
    im_shape = image.shape
    row_min = int(max(0, np.min(all_coords[:, 0]) - margin))
    row_max = int(min(im_shape[0], np.max(all_coords[:, 0]) + margin))
    col_min = int(max(0, np.min(all_coords[:, 1]) - margin))
    col_max = int(min(im_shape[1], np.max(all_coords[:, 1]) + margin))
    im_roi = im[row_min:row_max, col_min:col_max]

    im_roi = cv.cvtColor(im_roi, cv.COLOR_GRAY2RGB)
    plt.imshow(im_roi)
    plt.plot(spot_coords[:, 1] - col_min + 1, spot_coords[:, 0] - row_min + 1, 'rx', ms=8)
    plt.plot(grid_coords[:, 1] - col_min + 1, grid_coords[:, 0] - row_min + 1, 'b+', ms=8)
    plt.plot(reg_coords[:, 1] - col_min + 1, reg_coords[:, 0] - row_min + 1, 'g.', ms=8)
    plt.axis('off')
    fig_save = plt.gcf()
    fig_save.savefig(output_name + '_registration.png', bbox_inches='tight')
    plt.close(fig_save)
