# bchhun, {2020-03-29}
from ..utils.mock_regionprop import MockRegionprop
from skimage.feature import peak_local_max
import numpy as np
import math

from ..extract.image_parser import *


def fix_comets_thresh(props_):
    """
    attempt to
    :param props_:
    :return:
    """

    for idx, p in enumerate(props_):
        if p.eccentricity > 0.5:
            # find peak intensity pixel
            int_im = p.intensity_image
            centroid = p.local_centroid

            # perform threshold otsu on this
            # determine new centroid based on old centroid (old+new, old+new)
            # create mock prop with new int image, area, centroid

            region_thresh = thresh_and_binarize(int_im)
            region_props = generate_props(region_thresh)
            (newx, newy) = region_props.local_centroid

            # regionprops are iteratable, so if we want to copy all but replace centroid, we can do that
            mockprop = MockRegionprop(centroid=(newx+centroid[0], newy+centroid[1]))
            mockprop.area = 25*25

            props_[idx] = mockprop

    return props_


def fix_comets_min(props_):
    # for oblong objects, move the centroid the max intensity

    for idx, p in enumerate(props_):
        if p.eccentricity > 0.75:
            minint_coord = np.where(p.intensity_image == p.min_intensity)
            x_min, y_min, _, _ = p.bbox

            new_x, new_y = x_min+int(np.mean(minint_coord[0])), int(np.mean(y_min+minint_coord[1]))

            mockprop = MockRegionprop(centroid=(new_x, new_y))
            mockprop.area = 25 * 25

            props_[idx] = mockprop
        else:
            pass

    return props_
