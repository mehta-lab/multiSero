import cv2 as cv
import numpy as np


def create_reference_grid(mean_point,
                          nbr_grid_rows=6,
                          nbr_grid_cols=6,
                          spot_dist=83):
    """
    Generate initial spot grid based on image scale, center point, spot distance
    and spot layout (nbr rows and cols).

    :param tuple start_point: (x,y) coordinates of center of grid
    :param int nbr_grid_rows: Number of spot rows
    :param int nbr_grid_cols: Number of spot columns
    :param int spot_dist: Distance between spots
    :return np.array grid_coords: (x, y) coordinates for reference spots (nbr x 2)
    """
    start_x = mean_point[0] - spot_dist * (nbr_grid_cols - 1) / 2
    start_y = mean_point[1] - spot_dist * (nbr_grid_rows - 1) / 2
    x_vals = np.linspace(start_x, start_x + (nbr_grid_cols - 1) * spot_dist, nbr_grid_cols)
    y_vals = np.linspace(start_y, start_y + (nbr_grid_rows - 1) * spot_dist, nbr_grid_rows)
    grid_x, grid_y = np.meshgrid(x_vals, y_vals)
    grid_x = grid_x.flatten()
    grid_y = grid_y.flatten()
    grid_coords = np.vstack([grid_x.T, grid_y.T]).T

    return grid_coords


def icp(source, target, max_iterate=50, matrix_diff=1.):
    """
    Iterative closest point. Expects x, y coordinates of source and target in
    an array with shape: nbr of points x 2

    :param np.array source: Source spot coordinates
    :param np.array target: Target spot coordinates
    :param int max_iterate: Maximum number of registration iterations
    :param float matrix_diff: Sum of absolute differences between transformation
        matrices after one iteration
    :return np.array t_matrix: 2D transformation matrix (2 x 3)
    """
    src = source.copy().astype(np.float32)
    dst = target.copy().astype(np.float32)

    src = np.expand_dims(src, 0)
    dst = np.expand_dims(dst, 0)

    # Initialize kNN module
    knn = cv.ml.KNearest_create()
    labels = np.array(range(dst.shape[1])).astype(np.float32)
    knn.train(dst[0], cv.ml.ROW_SAMPLE, labels)
    # Initialize transformation matrix
    t_matrix = np.eye(3)
    t_temp = np.eye(3)
    t_old = t_matrix

    # Iterate while matrix difference > threshold
    for i in range(max_iterate):

        # Find closest points
        ret, results, neighbors, dist = knn.findNearest(src[0], 1)
        # Outlier removal
        idxs = np.squeeze(neighbors.astype(np.uint8))
        dist_max = 2 * np.median(dist)
        normal_idxs = np.where(dist < dist_max)[0]
        idxs = idxs[normal_idxs]
        # Find rigid transform
        t_iter = cv.estimateRigidTransform(
            src[0, normal_idxs, :],
            dst[0, idxs, :],
            fullAffine=False,
        )
        if t_iter is None:
            print("ICP optimization failed.")
            return None
        t_temp[:2] = t_iter
        src = cv.transform(src, t_iter)
        t_matrix = np.dot(t_temp, t_matrix)
        # Estimate diff
        t_diff = sum(sum(abs(t_matrix[:2] - t_old[:2])))
        t_old = t_matrix
        print(t_diff)
        if t_diff < matrix_diff:
            break

    return t_matrix[:2]
