import cv2 as cv
import numpy as np
import os
import pytest


@pytest.fixture(scope="session")
def image_dir(tmpdir_factory):
    im = np.zeros((5, 10), dtype=np.uint8)
    im[2:3, 3:8] = 128
    input_dir = tmpdir_factory.mktemp("input_dir")
    # Save a couple of images
    cv.imwrite(os.path.join(input_dir, 'A1.png'), im)
    cv.imwrite(os.path.join(input_dir, 'A2.png'), im)
    cv.imwrite(os.path.join(input_dir, 'B11.png'), im + 50)
    cv.imwrite(os.path.join(input_dir, 'B12.png'), im + 100)
    return input_dir


# @pytest.fixture(scope="session")
# def micromanager_dir(tmpdir_factory):
#     im = np.zeros((5, 10), dtype=np.uint8)
#     im[2:3, 3:8] = 128
#     input_dir = tmpdir_factory.mktemp("input_dir")
#     well_names = ['A1', 'A2', 'B11', 'B12']
#     # Save a couple of images
#     for well in well_names:
#         sub_dir_name = os.path.join("input_dir", well + "-what-ever")
#         sub_dir = tmpdir_factory.mktemp(sub_dir_name)
#         cv.imwrite(os.path.join(sub_dir, 'mm_name.tif'), im)
#     return input_dir
