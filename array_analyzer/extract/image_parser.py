# bchhun, {2020-03-22}

import os
from copy import copy
import numpy as np
import re

import skimage.io as io
import skimage.util as u

from skimage.color import rgb2grey
from skimage.filters import threshold_minimum, threshold_otsu
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.feature import canny
from skimage.morphology import binary_closing, binary_dilation, selem, disk, binary_opening
from scipy.ndimage import binary_fill_holes
from skimage import measure

from .img_processing import  create_unimodal_mask

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


def read_to_grey(path_, wellimage_):
    """
    a generator that receives file path and returns the next rgb image as greyscale and its name

    :param path_: path to folder with all images
    :param wellimage_: name of the file with image of the well.
    :return: next image as greyscale np.ndarray, filename
    """

    image_path = path_+os.sep+wellimage_
    im = io.imread(image_path)
    i = rgb2grey(im)
    return i, os.path.basename(image_path)


def thresh_and_binarize(image_, method='rosin', invert=True):
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

    if invert:
        image_ = u.invert(image_)

    if method == 'bimodal':
        thresh = threshold_minimum(image_, nbins=512)

        spots = copy(image_)
        spots[image_ < thresh] = 0
        spots[image_ >= thresh] = 1

    elif method == 'otsu':
        thresh = threshold_otsu(image_, nbins=512)

        spots = copy(image_)
        spots[image_ < thresh] = 0
        spots[image_ >= thresh] = 1

    elif method == 'rosin':
        spots = create_unimodal_mask(image_, str_elem_size=3)
    else:
        raise ModuleNotFoundError("not a supported method for thresh_and_binarize")

    return spots


def find_well_border(image, segmethod='bimodal', detmethod='region'):
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
    segmented_img = thresh_and_binarize(image, method=segmethod, invert=True)
    well_mask = segmented_img == 0
    # Now remove small objects.
    str_elem_size=10
    str_elem = disk(str_elem_size)
    well_mask = binary_opening(well_mask, str_elem)
    well_mask = binary_fill_holes(well_mask)


    if detmethod == 'region':
        labels = measure.label(well_mask)
        props = measure.regionprops(labels)

        # let's assume ONE circle for now (take only props[0])
        cy, cx = props[0].centroid # notice that the coordinate order is different from hough.
        radii = int((props[0].minor_axis_length + props[0].major_axis_length)/ 4 / np.sqrt(2))
        # Otsu threshold fails occasionally and leads to asymmetric region. Averaging both axes makes the segmentation robust.
        # If above files, try bounding box.

    elif detmethod == 'hough':
        hough_radii = [300, 400, 500, 600]

        well_mask = thresh_and_binarize(image, method='bimodal')

        edges = canny(well_mask, sigma=3)
        hough_res = hough_circle(edges, hough_radii)
        aaccums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=1)
        cx, cy = cx[0], cy[0]
        radii = radii[0]
    else:
        cx, cy, radii = None, None, None

    return cx, cy, radii, well_mask


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
    cx_=int(np.rint(cx_))
    cy_=int(np.rint(cy_))
    crop = arr[
           cy_ - (radius_ - border_): cy_ + (radius_ - border_),
           cx_ - (radius_ - border_): cx_ + (radius_ - border_)
           ]

    return crop


def clean_spot_binary(arr, kx=10, ky=10):
    return binary_closing(arr, selem=np.ones((kx, ky)))


def generate_spot_background(spotmask,distance=3,annulus=5):
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
    se_inner=selem.disk(distance,dtype=bool)
    se_outer=selem.disk(distance+annulus,dtype=bool)
    inner=binary_dilation(spotmask,se_inner)
    outer=binary_dilation(spotmask,se_outer)
    spot_background=np.bitwise_xor(inner,outer)

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


def select_props(props_, attribute, condition, condition_value):
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
    elif condition == 'is_in':
        props = [p for p in props_ if getattr(p, attribute) in condition_value]
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

def compute_od(props_array,bgprops_array):
    """
    
    Parameters
    ----------
    props_array: object: 
     2D array of regionprops objects at the spots over data.
    bgprops_array: object: 
     2D array of regionprops objects at the spots over background.

    Returns
    -------
    od_norm
    i_spot
    i_bg
    """
    assert props_array.shape == bgprops_array.shape, 'regionprops arrays representing sample and background are not the same.'
    rows=props_array.shape[0]
    cols=props_array.shape[1]
    i_spot=np.empty((rows,cols))
    i_bg=np.empty((rows,cols))
    od_norm=np.empty((rows,cols))

    i_spot[:]=np.NaN
    i_bg[:]=np.NaN
    od_norm[:]=np.NaN

    for r in np.arange(rows):
        for c in np.arange(cols):
            if props_array[r,c] is not None:
                i_spot[r,c]=props_array[r,c].mean_intensity
                i_bg[r,c]=bgprops_array[r,c].mean_intensity
    od_norm=i_bg/i_spot

    return od_norm, i_spot, i_bg
