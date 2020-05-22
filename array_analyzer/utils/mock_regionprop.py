import numpy as np
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
