from skimage.morphology import disk
import numpy as np
class MockRegionprop:

    def __init__(self,
                 intensity_image,
                 centroid,
                 image=None,
                 mean_intensity=None,
                 label=None,
                 bbox=None,
                 ):
        self.centroid = centroid
        self.mean_intensity = mean_intensity
        self.label = label
        self.intensity_image = intensity_image
        if image is None:
            image = disk(int(intensity_image.shape[0] / 2), dtype=np.bool_)
        self.image = image
        self.mean_intensity = np.mean(intensity_image[image])
        self.median_intensity = np.median(intensity_image[image])
        self.intensity_image = intensity_image * image
        self.bbox = bbox
