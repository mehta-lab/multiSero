class MockRegionprop:

    def __init__(self,
                 centroid=None,
                 mean_intensity=None,
                 label=None,
                 intensity_image=None,
                 ):
        self.centroid = centroid
        self.mean_intensity = mean_intensity
        self.label = label
        self.intensity_image = intensity_image
