import numpy as np
import pandas as pd
from skimage import measure
from skimage.morphology import disk


class SpotRegionprop:

    def __init__(self, label=None):
        """
        Object holding spot images, masks, and their properties:
        centroid, bounding box, mean and median intensity.

        :param int label: Spot label
        """
        self.image = None
        self.centroid = None
        self.label = label
        self.bbox = None
        self.mask = None
        # Values to be computed
        self.mean_intensity = None
        self.median_intensity = None
        self.masked_image = None

    def assign_median_masked(self):
        self.median_intensity = np.median(self.image[self.mask > 0])
        self.masked_image = self.image * self.mask

    @staticmethod
    def make_mask(im_size):
        """
        Creates disk shaped mask the size of image.

        :param int im_size: Image size (assume square shape)
        :return np.array mask: Binary disk shaped mask
        """
        mask = disk(int(im_size / 2), dtype=np.uint8)
        return mask

    def generate_props_from_disk(self, image, mask, bbox, centroid):
        """
        Assign properties from disk shaped mask.

        :param np.ndarray image: Large single spot ROI intensity image
        :param np.ndarray mask: Binary disk mask corresponding to image
        :param list bbox: Bounding box of single spot image
        :param list/tuple centroid: Center coordinate
        """
        self.image = image
        self.mask = mask
        self.centroid = centroid
        self.bbox = bbox
        self.mean_intensity = np.mean(self.image[self.mask > 0])
        self.assign_median_masked()

    def generate_props_from_mask(self, image, mask, bbox):
        """
        converts binarized image into region-properties using
        scikit-image.

        :param np.ndarray image: Large single spot ROI intensity image
        :param np.ndarray mask: Binary mask corresponding to image
        :param list bbox: Bounding box of single spot image
        """
        properties = ('label', 'centroid', 'mean_intensity',
                      'intensity_image', 'image', 'area', 'bbox')
        skimage_props = measure.regionprops_table(
            mask,
            intensity_image=image,
            properties=properties,
        )
        skimage_props = pd.DataFrame(skimage_props)

        self.mean_intensity = skimage_props.at[0, 'mean_intensity']

        self.centroid = [bbox[0] + skimage_props.at[0, 'centroid-0'],
                         bbox[1] + skimage_props.at[0, 'centroid-1']]

        min_row = skimage_props.at[0, 'bbox-0']
        min_col = skimage_props.at[0, 'bbox-1']
        max_row = skimage_props.at[0, 'bbox-2']
        max_col = skimage_props.at[0, 'bbox-3']

        self.bbox = [bbox[0] + min_row,
                     bbox[1] + min_col,
                     bbox[0] + max_row,
                     bbox[1] + max_col,
                     ]
        self.image = image[min_row:max_row, min_col:max_col]
        self.mask = mask[min_row:max_row, min_col:max_col]

        self.assign_median_masked()
