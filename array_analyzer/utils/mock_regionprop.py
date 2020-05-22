import numpy as np
import warnings
import pandas as pd
from skimage import measure
from skimage.morphology import disk


class MockRegionprop:

    def __init__(self,
                 image,
                 centroid,
                 label=None,
                 ):
        self.image = image
        self.centroid = centroid
        self.label = label
        self.bbox = None
        # Values to be computed in assign props
        self.mean_intensity = None
        self.mask = None
        self.mean_intensity = None
        self.median_intensity = None
        self.masked_image = None

    def _set_mean_intensity(self):
        self.mean_intensity = np.mean(self.intensity_image[self.image])

    def _set_median_intensity(self):
        self.median_intensity = np.median(self.intensity_image[self.image])

    def _set_intensity_image(self):
        self.intensity_image = self.intensity_image * self.image

    def _check_image_shape(self):
        if self.intensity_image == [] or self.intensity_image is None:
            warnings.warn("MockRegionProp received an empty intensity image")
            return False
        if self.image == [] or self.image is None:
            warnings.warn("MockRegionProp received an empty boolean mask")
            return False
        if self.intensity_image.shape != self.image.shape:
            warnings.warn("MockRegionProp received shape-mismatched mask and image")
            return False

    def assign_props(self, mask=None):
        if mask is None:
            self.mask = disk(int(self.image.shape[0] / 2), dtype=np.bool_)

        self.mean_intensity = np.mean(self.image[self.mask])
        self.median_intensity = np.median(self.image[self.mask])
        self.masked_image = self.image * self.mask

    def generate_props(self, image, mask, bbox):
        """
        converts binarized image into a list of region-properties using
        scikit-image first generates labels for the cleaned (binary_closing)
        binary image then generates regionprops on the remaining

        :param np.ndarray image: Large single spot ROI intensity image
        :param np.ndarray mask: Binary mask corresponding to image
        """
        self.mask = mask
        properties = ('label', 'centroid', 'mean_intensity',
                      'intensity_image', 'image', 'area', 'bbox')
        skimage_props = measure.regionprops_table(
            self.mask,
            intensity_image=image,
            properties=properties,
        )
        skimage_props = pd.DataFrame(skimage_props)

        self.mean_intensity = skimage_props.at[0, 'mean_intensity']
        self.image = image
        self.mask = mask
        self.median_intensity = np.median(self.image[self.mask])
        self.centroid = [bbox[0] + skimage_props.at[0, 'centroid-0'],
                         bbox[1] + skimage_props.at[0, 'centroid-1']]
        self.bbox = [bbox[0] + skimage_props.at[0, 'bbox-0'],
                     bbox[1] + skimage_props.at[0, 'bbox-1'],
                     bbox[0] + skimage_props.at[0, 'bbox-2'],
                     bbox[1] + skimage_props.at[0, 'bbox-3']
                     ]
