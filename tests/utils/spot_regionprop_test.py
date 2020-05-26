import numpy as np
import pytest

import array_analyzer.utils.spot_regionprop as regionprop


@pytest.fixture
def spot_and_mask():
    im_size = 51
    sigma = 10
    mu = 25
    row, col = np.meshgrid(
        np.linspace(0, im_size, im_size),
        np.linspace(0, im_size, im_size),
    )
    im_spot = np.exp(-((row - mu) / sigma) ** 2 / 2 - ((col - mu) / sigma) ** 2 / 2)
    # Invert for dark spot
    im_spot = 1. - im_spot
    mask_spot = (im_spot < .5).astype(np.uint8)
    return im_spot, mask_spot


def test_regionprop_init():
    prop = regionprop.SpotRegionprop(label=5)
    assert prop.label == 5
    assert prop.image is None
    assert prop.centroid is None
    assert prop.bbox is None
    assert prop.mask is None
    assert prop.mean_intensity is None
    assert prop.median_intensity is None
    assert prop.masked_image is None


def test_assign_median_masked(spot_and_mask):
    im_spot, mask_spot = spot_and_mask
    prop = regionprop.SpotRegionprop(label=1)
    prop.image = im_spot
    prop.mask = mask_spot
    prop.assign_median_masked()
    assert prop.label == 1
    assert prop.median_intensity < .3
    assert prop.masked_image.shape == im_spot.shape
    # Inside should have im spot values, edges should be zero
    assert prop.masked_image[0, 0] == 0.
    assert prop.masked_image[25, 25] == im_spot[25, 25]


def test_make_mask():
    prop = regionprop.SpotRegionprop(label=0)
    mask = prop.make_mask(51)
    assert mask.shape == (51, 51)
    assert mask[0, 0] == 0
    assert mask[25, 25] == 1


def test_generate_props_from_disk(spot_and_mask):
    im_spot, _ = spot_and_mask
    prop = regionprop.SpotRegionprop(label=3)
    mask = prop.make_mask(im_spot.shape[0])
    bbox = [5, 10, 15, 20]
    centroid = [25, 25]
    prop.generate_props_from_disk(im_spot, mask, bbox, centroid)
    assert prop.label == 3
    assert prop.image.all() == im_spot.all()
    assert prop.centroid == centroid
    assert prop.bbox == bbox
    assert prop.mask.all() == mask.all()
    assert prop.mean_intensity == np.mean(im_spot[mask > 0])
    assert prop.median_intensity == np.median(im_spot[mask > 0])
    assert prop.masked_image.all() == (im_spot * mask).all()


def test_generate_props_from_mask(spot_and_mask):
    im_spot, mask_spot = spot_and_mask
    prop = regionprop.SpotRegionprop(label=4)
    bbox = [5, 10, 15, 20]
    prop.generate_props_from_mask(im_spot, mask_spot, bbox)
    assert prop.label == 4
    assert prop.image.all() == im_spot.all()
    assert prop.centroid == [24 + bbox[0], 24 + bbox[1]]
    # There are the coordinates for the mask bounding box
    assert prop.bbox == [18, 23, 42, 47]
    assert prop.mask.all() == mask_spot.all()
    assert prop.mean_intensity == np.mean(im_spot[mask_spot > 0])
    assert prop.median_intensity == np.median(im_spot[mask_spot > 0])
    assert prop.masked_image.all() == (im_spot * mask_spot).all()
