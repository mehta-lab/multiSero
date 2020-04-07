import cv2 as cv
import numpy as np

from scipy.signal import find_peaks
from scipy.ndimage.filters import gaussian_filter1d
from skimage.morphology import disk, ball, binary_opening, binary_erosion
from skimage.filters import threshold_otsu, threshold_multiotsu
from scipy.ndimage import binary_fill_holes

from .background_estimator import BackgroundEstimator2D


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


def get_background(img, fit_order):
    bg_estimator = BackgroundEstimator2D(block_size=128)
    background = bg_estimator.get_background(img, order=fit_order, normalize=False)
    return background


def get_spot_coords(im,
                    min_thresh=0,
                    max_thresh=255,
                    min_area=50,
                    max_area=10000,
                    min_circularity=0.1,
                    min_convexity=0.5):
    """
    Use OpenCVs simple blob detector (thresholdings and grouping by properties)
    to detect all dark spots in the image
    :param np.array im: uint8 mage containing spots
    :param int min_thresh: Minimum threshold
    :param int max_thresh: Maximum threshold
    :param int min_area: Minimum spot area in pixels
    :param int max_area: Maximum spot area in pixels
    :param float min_circularity: Minimum circularity of spots
    :param float min_convexity: Minimum convexity of spots
    :return np.array spot_coords: x, y coordinates of spot centroids (nbr spots x 2)
    """
    params = cv.SimpleBlobDetector_Params()

    # Change thresholds
    params.minThreshold = min_thresh
    params.maxThreshold = max_thresh
    # Filter by Area
    params.filterByArea = True
    params.minArea = min_area
    params.maxArea = max_area
    # Filter by Circularity
    params.filterByCircularity = True
    params.minCircularity = min_circularity
    # Filter by Convexity
    params.filterByConvexity = True
    params.minConvexity = min_convexity

    detector = cv.SimpleBlobDetector_create(params)

    # Normalize image
    im_norm = ((im - im.min()) / (im.max() - im.min()) * 255).astype(np.uint8)
    # Detect blobs
    keypoints = detector.detect(im_norm)

    spot_coords = np.zeros((len(keypoints), 2))
    # Convert to np.arrays
    for c in range(len(keypoints)):
        pt = keypoints[c].pt
        spot_coords[c, 0] = pt[0]
        spot_coords[c, 1] = pt[1]

    return spot_coords


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
    x_min = int(max(margin, np.min(grid_coords[:, 0]) - margin))
    x_max = int(min(im_shape[1] - margin, np.max(grid_coords[:, 0]) + margin))
    y_min = int(max(margin, np.min(grid_coords[:, 1]) - margin))
    y_max = int(min(im_shape[0] - margin, np.max(grid_coords[:, 1]) + margin))
    im_crop = im[y_min:y_max, x_min:x_max]

    # Update coordinates with new origin
    crop_coords = grid_coords.copy()
    crop_coords[:, 0] = crop_coords[:, 0] - x_min + 1
    crop_coords[:, 1] = crop_coords[:, 1] - y_min + 1
    return im_crop, crop_coords
