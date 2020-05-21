from copy import copy

import cv2 as cv
import numpy as np

from scipy.signal import find_peaks
from scipy.ndimage.filters import gaussian_filter1d
from skimage.measure import label
from skimage import util as u
from skimage.morphology import disk, ball, binary_opening, binary_erosion
from skimage.filters import threshold_otsu, threshold_multiotsu, threshold_minimum
from scipy.ndimage import binary_fill_holes
from skimage.segmentation import clear_border


def get_unimodal_threshold(input_image):
    """Determines optimal unimodal threshold

    https://users.cs.cf.ac.uk/Paul.Rosin/resources/papers/unimodal2.pdf
    https://www.mathworks.com/matlabcentral/fileexchange/45443-rosin-thresholding

    :param np.array input_image: generate mask for this image
    :return float best_threshold: optimal lower threshold for the foreground
     hist
    """

    hist_counts, bin_edges = np.histogram(
        input_image,
        bins=256,
        range=(input_image.min(), np.percentile(input_image, 99.5))
    )
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # assuming that background has the max count
    max_idx = np.argmax(hist_counts)
    int_with_max_count = bin_centers[max_idx]
    p1 = [int_with_max_count, hist_counts[max_idx]]

    # find last non-empty bin
    pos_counts_idx = np.where(hist_counts > 0)[0]
    last_binedge = pos_counts_idx[-1]
    p2 = [bin_centers[last_binedge], hist_counts[last_binedge]]

    best_threshold = -np.inf
    max_dist = -np.inf
    for idx in range(max_idx, last_binedge, 1):
        x0 = bin_centers[idx]
        y0 = hist_counts[idx]
        a = [p1[0] - p2[0], p1[1] - p2[1]]
        b = [x0 - p2[0], y0 - p2[1]]
        cross_ab = a[0] * b[1] - b[0] * a[1]
        per_dist = np.linalg.norm(cross_ab) / np.linalg.norm(a)
        if per_dist > max_dist:
            best_threshold = x0
            max_dist = per_dist
    assert best_threshold > -np.inf, 'Error in unimodal thresholding'
    return best_threshold


def create_unimodal_mask(input_image, str_elem_size=3):
    """Create a mask with unimodal thresholding and morphological operations

    unimodal thresholding seems to oversegment, erode it by a fraction

    :param np.array input_image: generate masks from this image
    :param int str_elem_size: size of the structuring element. typically 3, 5
    :return: mask of input_image, np.array
    """

    if np.min(input_image) == np.max(input_image):
        thr = np.unique(input_image)
    else:
        thr = get_unimodal_threshold(input_image)
    if len(input_image.shape) == 2:
        str_elem = disk(str_elem_size)
    else:
        str_elem = ball(str_elem_size)
    # remove small objects in mask
    thr_image = binary_opening(input_image > thr, str_elem)
    mask = binary_erosion(thr_image, str_elem)
    return mask


def create_otsu_mask(input_image, scale=1):
    """Create a binary mask using morphological operations

    Opening removes small objects in the foreground.

    :param np.array input_image: generate masks from this image
    :param float scale: Scale the threshold
    :return: mask of input_image, np.array
    """
    if np.min(input_image) == np.max(input_image):
        return np.ones(input_image.shape)
    else:
        thr = threshold_otsu(input_image, nbins=512)
    return input_image > (scale * thr)


def create_multiotsu_mask(input_image, n_class, fg_class, str_elem_size=3):
    """Create a binary mask using morphological operations

    Opening removes small objects in the foreground.

    :param np.array input_image: generate masks from this image
    :param int str_elem_size: size of the structuring element. typically 3, 5
    :return: mask of input_image, np.array
    """
    if np.min(input_image) == np.max(input_image):
        return np.ones(input_image.shape)
    else:
        thr = threshold_multiotsu(input_image, classes=n_class, nbins=512)
        im_label = np.digitize(input_image, bins=thr)
        mask = im_label == fg_class

    if len(input_image.shape) == 2:
        str_elem = disk(str_elem_size)
    else:
        str_elem = ball(str_elem_size)
    # remove small objects in mask
    mask = binary_opening(mask, str_elem)
    return mask


def find_profile_peaks(profile, margin, prominence):
    # invert because black spots
    profile = profile.max() - profile
    max_pos = int(np.where(profile == profile.max())[0][0])
    # Make sure max is not due to leaving the center
    add_margin = 0
    half_margin = int(margin / 2)
    if max_pos > len(profile) - half_margin:
        profile = profile[:-half_margin]
    elif max_pos < half_margin:
        profile = profile[half_margin:]
        add_margin = half_margin
    profile = gaussian_filter1d(profile, 3)
    min_prom = profile.max() * prominence
    peaks, _ = find_peaks(profile, prominence=min_prom, distance=50)
    if len(peaks) >= 4:
        spot_dists = peaks[1:] - peaks[:-1]
    else:
        spot_dists = None
    mean_pos = peaks[0] + (peaks[-1] - peaks[0]) / 2 + add_margin
    return mean_pos, spot_dists


def grid_estimation(im,
                    spot_coords,
                    margin=50,
                    prominence=.15):
    """
    Based on images intensities and detected spots, make an estimation
    of grid location so that ICP algorithm is initialized close enough for convergence.
    TODO: This assumes that you always detect the first peaks
    this may be unstable so think of other ways to initialize...
    :param np.array im: Grayscale image
    :param np.array spot_coords: Spot x,y coordinates (nbr spots x 2)
    :param int margin: Margin for cropping outside all detected spots
    :param float prominence: Fraction of max intensity to filter out insignificant peaks
    :return tuple start_point: Min x, y coordinates for initial grid estimate
    :return float spot_dist: Estimated distance between spots
    """
    im_shape = im.shape
    x_min = int(max(margin, np.min(spot_coords[:, 0]) - margin))
    x_max = int(min(im_shape[1] - margin, np.max(spot_coords[:, 0]) + margin))
    y_min = int(max(margin, np.min(spot_coords[:, 1]) - margin))
    y_max = int(min(im_shape[0] - margin, np.max(spot_coords[:, 1]) + margin))
    im_roi = im[y_min:y_max, x_min:x_max]
    # Create intensity profiles along x and y and find peaks
    profile_x = np.mean(im_roi, axis=0)
    mean_x, dists_x = find_profile_peaks(profile_x, margin, prominence)
    profile_y = np.mean(im_roi, axis=1)
    mean_y, dists_y = find_profile_peaks(profile_y, margin, prominence)

    mean_point = (x_min + mean_x, y_min + mean_y)
    spot_dist = np.hstack([dists_x, dists_y])
    # Remove invalid distances
    spot_dist = spot_dist[np.where(spot_dist != None)]
    if spot_dist.size == 0:
        # Failed at estimating spot dist. Return default or error out?
        spot_dist = None
    else:
        spot_dist = np.median(spot_dist)

    return mean_point, spot_dist


def crop_image_from_coords(im, grid_coords, margin=200):
    """
    Given image coordinates, crop image around them with a margin.

    :param np.array im: 2D image
    :param np.array grid_coords: Fitted grid coordinates which should be contained
        in the cropped image (nbr points x 2)
    :param int margin: How much margin around the coordinates
    :return np.array im_roi: Cropped image
    :return np.array grid_coords: Grid coordinates with new origin
    """
    im_shape = im.shape
    x_min = int(max(0, np.min(grid_coords[:, 0]) - margin))
    x_max = int(min(im_shape[1], np.max(grid_coords[:, 0]) + margin))
    y_min = int(max(0, np.min(grid_coords[:, 1]) - margin))
    y_max = int(min(im_shape[0], np.max(grid_coords[:, 1]) + margin))
    im_crop = im[y_min:y_max, x_min:x_max]

    # Update coordinates with new origin
    crop_coords = grid_coords.copy()
    crop_coords[:, 0] = crop_coords[:, 0] - x_min + 1
    crop_coords[:, 1] = crop_coords[:, 1] - y_min + 1
    return im_crop, crop_coords


def crop_image_at_center(im, center, height, width):
    """
    Crop the supplied image to include only the well and its spots

    :param im: image
    :param float center: Center (row, col) of the crop box
    :param float height: height of the crop box
    :param float width: width of the crop box
    :return np.array crop: Cropped image
    :return list bbox: Bounding box coordinates [row min, col min, row max, col max]
    """
    c_row, c_col = center
    im_h, im_w = im.shape
    # truncate the bounding box when it exceeds image size
    bbox = np.rint([max(c_row - height / 2, 0),
                   max(c_col - width / 2, 0),
                   min(c_row + height / 2, im_h),
                   min(c_col + width / 2, im_w)]).astype(np.int32)
    crop = im[bbox[0]:bbox[2], bbox[1]:bbox[3]]
    return crop, bbox


def crop_image(arr, cx_, cy_, radius_, border_=200):
    """
    crop the supplied image to include only the well and its spots

    :param arr: image
    :param float cx_: Center x coordinate
    :param float cy_: Center y coordinate
    :param float radius_: Crop radius
    :param int border_: Margin on each side in pixels
    :return np.array crop: Cropped image
    """
    cx_ = int(np.rint(cx_))
    cy_ = int(np.rint(cy_))
    crop = arr[
           cy_ - (radius_ - border_): cy_ + (radius_ - border_),
           cx_ - (radius_ - border_): cx_ + (radius_ - border_)
           ]

    return crop


def get_largest_component(spot_segm):
    """
    Remove everything but the largest connected component from segmented image.

    :param np.array spot_segm: Binary segmented 2D image
    :return np.array largest_component: Largest connected component in image
    """
    labels = label(spot_segm)
    largest_component = labels.copy()
    if labels.max() > 0:
        largest_component = labels == np.argmax(np.bincount(labels.flat)[1:]) + 1
    return largest_component.astype(np.uint8)


def thresh_and_binarize(image,
                        method='rosin',
                        invert=True,
                        min_size=10,
                        thr_percent=95):
    """
    receives greyscale np.ndarray image
        inverts the intensities
        thresholds on the minimum peak
        converts the image into binary about that threshold

    :param np.ndarray image: 2D grayscale image
    :param str method: Trheshold type: 'bimodal', 'otsu', 'rosin' or 'bright_spots'
    :param bool invert: Invert image if spots are dark
    :param int min_size: Minimum structuring element disk size
    :param int thr_percent: Thresholding percentile
    :return: spots threshold_min on this image
    """
    image_ = image.copy()
    if invert:
        image_ = u.invert(image_)

    if method == 'bimodal':
        thresh = threshold_minimum(image_, nbins=512)
        spots = (image_ > thresh).astype(np.uint8)

    elif method == 'otsu':
        spots = create_otsu_mask(image_, scale=1)

    elif method == 'rosin':
        spots = create_unimodal_mask(image_, str_elem_size=3)

    elif method == 'bright_spots':
        spots = image_ > np.percentile(image_, thr_percent)
        str_elem = disk(min_size)
        spots = binary_opening(spots, str_elem)
        spots = clear_border(spots)

    else:
        raise ModuleNotFoundError("not a supported method for thresh_and_binarize")
    # Keep only largest connected component
    spots = get_largest_component(spots)

    return spots


class SpotDetector:
    """
    Detects spots in well image using a Laplacian of Gaussian filter
    followed by blob detection
    """

    def __init__(self,
                 imaging_params,
                 min_thresh=100,
                 max_thresh=255,
                 max_area=10000,
                 min_circularity=.1,
                 min_convexity=.5,
                 min_dist_between_blobs=10,
                 min_repeatability=2,
                 max_intensity=255):
        """
        :param int min_thresh: Minimum threshold
        :param int max_thresh: Maximum threshold
        :param int min_area: Minimum spot area in pixels
        :param int max_area: Maximum spot area in pixels
        :param float min_circularity: Minimum circularity of spots
        :param float min_convexity: Minimum convexity of spots
        :param float min_dist_between_blobs: minimal distance in pixels between two
            spots for them to be called as different spots
        :param int min_repeatability: minimal number of times the same spot has to be
            detected at different thresholds
        :param int max_intensity: Maximum image intensity (default uint8)
        """

        self.min_thresh = min_thresh
        self.max_thresh = max_thresh
        self.min_dist_between_blobs = min_dist_between_blobs
        self.min_repeatability = min_repeatability
        self.min_circularity = min_circularity
        self.min_convexity = min_convexity
        self.max_intensity = max_intensity
        self.sigma_gauss = int(np.round(imaging_params['spot_width'] /
                                imaging_params['pixel_size'] / 4))
        self.min_area = 4 * self.sigma_gauss ** 2
        self.max_area = max_area
        self.nbr_expected_spots = imaging_params['rows'] * imaging_params['columns']

        self.blob_detector = self._make_blob_detector()
        self.log_filter = self._make_log_filter()

    def _make_blob_detector(self):
        # Set spot detection parameters
        blob_params = cv.SimpleBlobDetector_Params()
        # Change thresholds
        blob_params.minThreshold = self.min_thresh
        blob_params.maxThreshold = self.max_thresh
        # Filter by Area
        blob_params.filterByArea = True
        blob_params.minArea = self.min_area
        blob_params.maxArea = self.max_area
        # Filter by Circularity
        blob_params.filterByCircularity = True
        blob_params.minCircularity = self.min_circularity
        # Filter by Convexity
        blob_params.filterByConvexity = True
        blob_params.minConvexity = self.min_convexity
        blob_params.minDistBetweenBlobs = self.min_dist_between_blobs
        blob_params.minRepeatability = self.min_repeatability
        # This detects bright spots, which they are after top hat
        blob_params.blobColor = self.max_intensity
        detector = cv.SimpleBlobDetector_create(blob_params)
        return detector

    def _make_log_filter(self):
        """
        Creates a uniform 2D Laplacian of Gaussian filter with given sigma.

        :param int sigma: Standard deviation of Gaussian
        :return np.array log_filter: 2D LoG filter
        """
        n = np.ceil(self.sigma_gauss * 6)
        y, x = np.ogrid[-n // 2:n // 2 + 1, -n // 2:n // 2 + 1]
        sigma_sq = 2 * self.sigma_gauss ** 2
        y_filter = np.exp(-(y ** 2 / sigma_sq))
        x_filter = np.exp(-(x ** 2 / sigma_sq))
        log_filter = (-sigma_sq + x ** 2 + y ** 2) * \
                     (x_filter * y_filter) * \
                     (1 / (np.pi * sigma_sq * self.sigma_gauss ** 2))
        # Total filter should sum to 1 to not alter mean intensity
        log_filter = log_filter / sum(sum(log_filter))
        return log_filter

    def get_spot_coords(self, im, margin=0, im_mean=100, im_std=25):
        """
        Use OpenCVs simple blob detector (thresholdings and grouping by properties)
        to detect all dark spots in the image. First filter with a Laplacian of
        Gaussian with sigma matching spots to enhance spots in image.

        :param np.array im: uint8 mage containing spots
        :param int margin: Pixel margin around image edged where spots should be
            ignored (to ignore boundary effects)
        :param float im_mean: Set normalized image to fixed mean
        :param float im_std: Set normalized image to fixed std
        :return np.array spot_coords: x, y coordinates of spot centroids
            (nbr spots x 2)
        """
        # First invert image to detect peaks
        im_norm = (im.max() - im) / self.max_intensity
        # Filter with Laplacian of Gaussian
        im_norm = cv.filter2D(im_norm, -1, self.log_filter)
        # Normalize
        im_norm = im_norm / im_norm.std() * im_std
        im_norm = im_norm - im_norm.mean() + im_mean
        im_norm[im_norm < 0] = 0
        im_norm[im_norm > self.max_intensity] = self.max_intensity
        im_norm = im_norm.astype(np.uint8)

        # Detect peaks in filtered image
        keypoints = self.blob_detector.detect(im_norm)

        spot_coords = np.zeros((len(keypoints), 2))
        # Remove outliers and convert to np.array
        x_max, y_max = im.shape
        idx = 0
        for keypoint in range(len(keypoints)):
            pt = keypoints[keypoint].pt
            if margin < pt[0] < x_max - margin and margin < pt[1] < y_max - margin:
                spot_coords[idx, 0] = pt[0]
                spot_coords[idx, 1] = pt[1]
                idx += 1
        spot_coords = spot_coords[:idx, :]
        return spot_coords
