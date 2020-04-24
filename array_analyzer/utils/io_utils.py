import cv2 as cv
import glob
import natsort
import os
import re
import skimage.io as io
from skimage.color import rgb2grey


def read_to_grey(path_, wellimage_):
    """
    a generator that receives file path and returns the next rgb image as greyscale and its name

    :param str path_: path to folder with all images
    :param str wellimage_: name of the file with image of the well.
    :return: next image as greyscale np.ndarray
    """
    image_path = os.path.join(path_, wellimage_)
    im = io.imread(image_path)
    im = rgb2grey(im)
    return im, os.path.basename(image_path)


def read_gray_im(im_path):
    """
    Read image from full path to file location.
    :param str im_path: Path to image
    :return np.array im: Grayscale image
    """
    try:
        im = cv.imread(im_path, cv.IMREAD_GRAYSCALE | cv.IMREAD_ANYDEPTH)
    except IOError as e:
        raise("Can't read image", e)
    return im


def get_image_paths(input_dir):
    """
    Searches input directory and its subdirectory for all images
    and returns sorted list of paths.

    :param str input_dir: Input directory, may contain images or subdirectories
        with one image each
    :return dict well_images: Well name key, path to found image value
    """
    extensions = ('*.png', '*.tif', '*.jpg')

    image_names = []
    for ext in extensions:
        search_str = os.path.join(input_dir, '**', ext)
        image_names.extend(glob.glob(search_str, recursive=True))

    # Sort images
    image_names = natsort.natsorted(image_names)
    # Find well names from image paths
    well_images = {}
    if len(os.path.basename(image_names[0])) > 6:
        # Micromanager naming convention, find well from subdir name
        for im in image_names:
            well_name = im.split('/')[-2]
            # split again for well name, assume - separation
            well_name = well_name.split('-')[0]
            well_images[well_name] = im
    else:
        # Assume images are named e.g. A0.png
        for im in image_names:
            well_name = os.path.basename(im)[:-4]
            well_images[well_name] = im

    return well_images
