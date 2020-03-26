# bchhun, {2020-03-22}

import os
import numpy as np
import re

import skimage.io as io
import skimage.util as u

from skimage.color import rgb2grey
from skimage.filters import threshold_minimum, median, gaussian, threshold_local
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.feature import canny
from skimage.morphology import binary_closing, binary_dilation, selem, disk, binary_opening
from skimage import measure

from .img_processing import create_unimodal_mask, create_otsu_mask

"""
method is
1) read_to_grey(supplied images)
# find center of well
2) thresh and binarize from 1
3) find well border from 2
4) crop image from 3

# find center of spots from crop
5) thresh and binarize from 4
6) clean spot binary from 5
7) generate props from 6
8) generate props dict from 7
9) assign props dict to array from 8
"""


def read_to_grey(path_):
    """
    a generator that receives file path and returns the next rgb image as greyscale and its name

    :param path_: path to folder with all images
    :return: next image as greyscale np.ndarray, filename
    """

    images = [file for file in os.listdir(path_) if '.png' in file or '.tif' in file or '.jpg' in file]
    # remove any images that are not images of wells.
    wellimages = [file for file in images if re.match(r'[A-P][0-9]{1,2}', file)]
    # sort by letter, then by number (with '10' coming AFTER '9')
    wellimages.sort(key=lambda x: (x[0], int(x[1:-4])))

    for image_base_path in wellimages:
        image_path = path_+os.sep+image_base_path
        im = io.imread(image_path)
        i = rgb2grey(im)
        yield i, os.path.basename(image_path)


def thresh_and_binarize(image_, method='rosin'):
    """
    receives greyscale np.ndarray image
        inverts the intensities
        thresholds on the minimum peak
        converts the image into binary about that threshold

    :param image_: np.ndarray
    :param method: str
        'bimodal' or 'unimodal'
    :return: spots threshold_min on this image
    """
    inv = u.invert(image_)
    if method == 'bimodal':
        thresh = threshold_minimum(inv)
        inv = inv > thresh
        spots = inv.astype(int)

    elif method == 'rosin':
        spots = create_unimodal_mask(inv, str_elem_size=3)

    else:
        raise ModuleNotFoundError("not a supported method for thresh_and_binarize")

    return spots


def ivan_adaptive_threshold(image_):

    im_inv_crop = u.invert(image_)

    # median
    region = np.ones((10, 10))
    im1_filt = median(im_inv_crop, selem=region)

    # gaussian
    im1_filt = gaussian(im1_filt, sigma=10)

    # adaptive
    block = 15
    offset = -0.05
    method = 'gaussian'
    thr = threshold_local(im1_filt, block, method, offset)

    # binary
    im1_bw = im_inv_crop > thr
    str_elem = disk(5)
    im1_bw = binary_opening(im1_bw, selem=str_elem)

    return im1_bw


def find_well_center(image, method='otsu'):
    """
    finds the border of the well to motivate future cropping around spots
        hough_radii are potential radii of the well in pixels
            this should be motivated based on magnification
        edge filter
        fit hough circle
        find the peak of the SINGULAR hough circle

    :param image: np.ndarray
        raw image, not inverted
    :param method: str
        'otsu' or 'hough'
    :return: center x, center y, radius of the one hough circle
    """
    if method == 'otsu':
        well_mask = create_otsu_mask(image, str_elem_size=10)
        # plt.imshow(well_mask, cmap='gray')

        labels = measure.label(well_mask)
        props = measure.regionprops(labels)

        # let's assume ONE circle for now (take only props[0])
        cx, cy = props[0].centroid
        radii = int(props[0].minor_axis_length / 2 / np.sqrt(2))

    elif method == 'hough':
        hough_radii = [300, 400, 500, 600]

        binary_ = thresh_and_binarize(image, method='bimodal')

        edges = canny(binary_, sigma=3)
        hough_res = hough_circle(edges, hough_radii)
        aaccums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=1)
        cx, cy = cx[0], cy[0]
        radii = radii[0]
    else:
        cx, cy, radii = None, None, None

    return int(cx), int(cy), int(radii)


def crop_image(arr, cx_, cy_, radius_, border_=200):
    """
    crop the supplied image to include only the well and its spots

    :param arr: image
    :param cx_: float
    :param cy_: float
    :param radius_:
    :param border_:
    :return:
    """
    cx_ = int(np.rint(cx_))
    cy_ = int(np.rint(cy_))
    crop = arr[
           cx_ - (radius_ - border_): cx_ + (radius_ - border_),
           cy_ - (radius_ - border_): cy_ + (radius_ - border_)
           ]

    return crop


def clean_spot_binary(arr, kx=10, ky=10):
    return binary_closing(arr, selem=np.ones((kx, ky)))


def generate_spot_background(spotmask, distance=3, annulus=5):
    """
    
    compute an annulus around each spot to estimate background.
    
    Parameters
    ----------
    spotmask : binary mask of spots
    distance : distance from the edge of segmented spot.
    annulus : width of the annulus

    Returns
    -------
    spotbackgroundmask: binary mask of annuli around spots.
    
    TODO: 'comets' should be ignored, and this approach may not be robust to it.
    """
    se_inner = selem.disk(distance, dtype=bool)
    se_outer = selem.disk(distance+annulus, dtype=bool)
    inner = binary_dilation(spotmask, se_inner)
    outer = binary_dilation(spotmask, se_outer)
    spot_background = np.bitwise_xor(inner, outer)

    return spot_background


def generate_props(mask, intensity_image_=None):
    """
    converts binarized image into a list of region-properties using scikit-image
        first generates labels for the cleaned (binary_closing) binary image
        then generates regionprops on the remaining

    :param mask: np.ndarray
        binary version of cropped image
    :param intensity_image_: np.ndarray
        intensity image corresponding to this binary
    :return: list
        of skimage region-props object
    """
    labels = measure.label(mask)
    props = measure.regionprops(labels, intensity_image=intensity_image_)
    return props


def filter_props(props_, attribute, condition, condition_value):
    """

    :param props_: RegionProps
    :param attribute: str
        a regionprop attribute
        https://scikit-image.org/docs/dev/api/skimage.measure.html#regionprops
    :param condition: str
        one of "greater_than", "equals", "less_than"
    :param condition_value: int, float
        the value to evaluate around
    :return:
    """

    if condition == 'greater_than':
        props = [p for p in props_ if getattr(p, attribute) > condition_value]
    elif condition == 'equals':
        props = [p for p in props_ if getattr(p, attribute) == condition_value]
    elif condition == 'less_than':
        props = [p for p in props_ if getattr(p, attribute) < condition_value]
    else:
        props = props_

    return props


def generate_props_dict(props_, rows, cols, min_area=100, img_x_max=2048, img_y_max=2048):
    """
    based on the region props, creates a dictionary of format:
        key = (centroid_x, centroid_y)
        value = region_prop object

    :param props_: list of region props
        approximately 36-48 of these, depending on quality of the image
    :param rows: int
    :param cols: int
    :param min_area: int
    :param img_x_max: int
    :param img_y_max: int
    :return: dict
        of format (cent_x, cent_y): prop
    """

    # find minx, miny to "zero center" the array
    minx = img_x_max
    miny = img_y_max
    # find maxx, maxy to scale to array index values
    maxx = 0
    maxy = 0
    for prop in props_:
        if prop.area > min_area:
            if prop.centroid[0] < minx:
                minx = prop.centroid[0]
            if prop.centroid[1] < miny:
                miny = prop.centroid[1]
            if prop.centroid[0] > maxx:
                maxx = prop.centroid[0]
            if prop.centroid[1] > maxy:
                maxy = prop.centroid[1]

    # scaled max-x, max-y
    smaxx = maxx - minx
    smaxy = maxy - miny

    chk_list = []
    cent_map = {}
    for prop in props_:
        if prop.area > min_area:
            cx, cy = prop.centroid
            csx = cx - minx
            csy = cy - miny

            # convert the centroid position to an integer that maps to array indices
            norm_cent_x = int(round((rows-1) * (csx / smaxx)))
            norm_cent_y = int(round((cols-1) * (csy / smaxy)))

            # print(f"\ncentroid = {prop.centroid}\n\tnorm_cent = {norm_cent_x, norm_cent_y}")

            chk_list.append((norm_cent_x, norm_cent_y))
            cent_map[(norm_cent_x, norm_cent_y)] = prop

    if len(chk_list) != len(set(chk_list)):
        print("ERROR, DUPLICATE ENTRIES")
        raise AttributeError("generate props array failed\n"
                             "duplicate spots found in one position\n")

    return cent_map


def generate_region_array(props_):
    # props contains bounding box info
    # this will be useful to fill in NONE array values that have extremely low signal
    pass


def assign_props_to_array(arr, cent_map_):
    """
    takes an empty array and assigns region_property objects to each position, based on print array position

    :param arr: np.ndarray
        of shape = print array shape
    :param cent_map_: dict
        generated by "generate_props_array"
    :return:
    """

    for key, value in cent_map_.items():
        arr[key[0], key[1]] = value

    return arr


def assign_region(target_, props_, source_array_=None):
    """
    put underlying intensity image from props_ into the target_ array
    assumes props_ image exists and will fit into target_

    :param target_:
    :param props_:
    :return:
    """

    min_row, min_col, max_row, max_col = props_.bbox
    if source_array_ is None:
        target_[min_row:max_row, min_col:max_col] = props_.intensity_image
    else:
        target_[min_row:max_row, min_col:max_col] = source_array_[min_row:max_row, min_col:max_col]

    return target_


def create_composite_spots(target_array_, region_props_array_, source_array_=None):
    """
    insert all intensity images for each region prop in "region_props_array_" into the "target_array_"

    :param target_array_:
    :param region_props_array_:
    :return:
    """
    for row in range(region_props_array_.shape[0]):
        for col in range(region_props_array_.shape[1]):
            props = region_props_array_[row, col]
            if props is not None:
                if source_array_ is None:
                    # print(f"\t shape = {props.intensity_image.shape}")
                    target_array_ = assign_region(target_array_, props)
                else:
                    target_array_ = assign_region(target_array_, props, source_array_)
    return target_array_


