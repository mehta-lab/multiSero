from skimage.morphology import disk
import numpy as np
import warnings


class MockRegionprop:

    def __init__(self,
                 intensity_image,
                 centroid,
                 image=None,
                 label=None,
                 bbox=None,
                 ):
        self.centroid = centroid
        self.label = label
        self.intensity_image = intensity_image
        self.bbox = bbox
        if image is None:
            self.image = disk(int(intensity_image.shape[0] / 2), dtype=np.bool_)
        else:
            self.image = image

        if not self._check_image_shape():
            self.median_intensity = -1
        else:
            self._set_median_intensity()
            self._set_mean_intensity()
            self._set_intensity_image()

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
