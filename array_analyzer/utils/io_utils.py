import cv2 as cv
from datetime import datetime
import glob
import natsort
import numpy as np
import os
import skimage.io as io
from skimage.color import rgb2grey
import re


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
    if not isinstance(im, np.ndarray):
        raise IOError('Invalid image path')
    return im


def get_image_paths(input_dir):
    """
    Searches input directory and its subdirectory for all images
    and returns sorted list of paths.

    :param str input_dir: Input directory, may contain images or subdirectories
        with one image each
    :return dict well_images: Well name key, path to found image value
    """
    extensions = ('.png', '.tif')

    image_names = []
    for ext in extensions:
        search_str = os.path.join(input_dir,'[A-P][0-9]*'+ext)
        image_names.extend(glob.glob(search_str, recursive=False))

    # Sort images
    image_names = natsort.natsorted(image_names)
    # Find well names from image paths
    well_images = {}
    if len(image_names) > 0:
        # Assume images are named e.g. A0.png
        for im_name in image_names:
            well_name = os.path.basename(im_name)[:-4]
            if re.match(r'[A-P][0-9]{1,2}',well_name): #  double-check that the file represents a well.
                well_images[well_name] = im_name
    else:
        # Micromanager naming convention, find well from subdir name
        image_names = []
        for ext in extensions:
            search_str = os.path.join(input_dir,'[A-P][0-9]*','*'+ext)
            image_names.extend(glob.glob(search_str, recursive=False))

        # Sort images
        image_names = natsort.natsorted(image_names)
        for im_name in image_names:
            well_name = im_name.split('/')[-2]
            # split again for well name, assume - separation
            well_name = well_name.split('-')[0]
            if re.match(r'[A-P][0-9]{1,2}',well_name): #  double-check that the file represents a well.
                well_images[well_name] = im_name

    return well_images


def make_run_dir(input_dir, output_dir):
    """
    For a specific processing run, create a subdirectory in the output directory
    which specifies which input directory was used and when the processing took
    place.

    :param str input_dir: Path to input directory, to be processed
    :param str output_dir: Path to main output directory
    :return str run_dir: Path to directory where processed data is stored
    """
    run_dir = os.path.join(
        output_dir,
        '_'.join([os.path.basename(os.path.normpath(input_dir)),
                  str(datetime.now().month),
                  str(datetime.now().day),
                  str(datetime.now().hour),
                  str(datetime.now().minute),
                  str(datetime.now().second)]),
    )
    os.makedirs(run_dir, exist_ok=True)
    return run_dir

