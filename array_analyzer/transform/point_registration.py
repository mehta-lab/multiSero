import cv2 as cv
import numpy as np


def create_reference_grid(center_point,
                          nbr_grid_rows=6,
                          nbr_grid_cols=6,
                          spot_dist=83):
    """
    Generate initial spot grid based on image scale, center point, spot distance
    and spot layout (nbr rows and cols).

    :param tuple center_point: (row, col) coordinates of center of grid
    :param int nbr_grid_rows: Number of spot rows
    :param int nbr_grid_cols: Number of spot columns
    :param int spot_dist: Distance between spots
    :return np.array grid_coords: (row, col) coordinates for reference spots (nbr x 2)
    """
    start_row = center_point[0] - spot_dist * (nbr_grid_rows - 1) / 2
    start_col = center_point[1] - spot_dist * (nbr_grid_cols - 1) / 2
    row_vals = np.linspace(
        start_row,
        start_row + (nbr_grid_rows - 1) * spot_dist,
        nbr_grid_rows,
    )
    col_vals = np.linspace(
        start_col,
        start_col + (nbr_grid_cols - 1) * spot_dist,
        nbr_grid_cols,
    )
    grid_cols, grid_rows = np.meshgrid(col_vals, row_vals)
    grid_cols = grid_cols.flatten()
    grid_rows = grid_rows.flatten()
    grid_coords = np.vstack([grid_rows.T, grid_cols.T]).T

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

        if t_diff < matrix_diff:
            break

    return t_matrix[:2]


def create_gaussian_particles(stds,
                              mean_point=(0, 0),
                              scale_mean=1.,
                              angle_mean=0.,
                              nbr_particles=100):
    """
    Create particles from parameters x, y, scale and angle given mean and std.
    A particle is considered one set of parameters for a 2D translation matrix.

    :param np.array stds: Standard deviations of x, y, angle and scale
    :param tuple mean_point: Mean offset in x and y
    :param float scale_mean: Mean scale of translation
    :param float angle_mean: Mean angle of translation estimation
    :param int nbr_particles: Number of particles
    :return np.array particles: Set of particle coordinates (nbr particles x 4)
    """
    particles = np.empty((nbr_particles, 4))
    particles[:, 0] = mean_point[0] + (np.random.randn(nbr_particles) * stds[0])
    particles[:, 1] = mean_point[1] + (np.random.randn(nbr_particles) * stds[1])
    particles[:, 2] = angle_mean + (np.random.randn(nbr_particles) * stds[2])
    particles[:, 3] = scale_mean + (np.random.randn(nbr_particles) * stds[3])
    return particles


def get_translation_matrix(particle):
    """
    Create a 2D translation matrix from x, y, scale and angle.

    :param np.array particle: The four parameters x, y, scale and angle
    :return np.array t_matrix: 2D translation matrix (3 x 2)
    """
    a = particle[3] * np.cos(particle[2] * np.pi / 180)
    b = particle[3] * np.sin(particle[2] * np.pi / 180)
    t_matrix = np.array([[a, b, particle[0]],
                        [-b, a, particle[1]]])
    return t_matrix


def particle_filter(fiducial_coords,
                    spot_coords,
                    particles,
                    stds,
                    max_iter=100,
                    stop_criteria=.1,
                    iter_decrease=.8,
                    remove_outlier=False,
                    debug=False):
    """
    Particle filtering to determine best grid location.
    Start with a number of randomly placed particles. Compute distances
    to nearest neighbors among detected spots. Do importance sampling of
    the particles and slightly distort them while iterating until convergence.

    :param np.array fiducial_coords: Initial estimate of fiducial coordinates
    :param np.array spot_coords: Coordinates of detected spots (nbr spots x 2)
    :param np.array particles: Translation matrix parameters (nbr particles x 2)
    :param np.array stds: Standard deviations of x, y, angle, scale
    :param int max_iter: Maximum number of iterations
    :param float stop_criteria: Absolute difference of distance between iterations
    :param float iter_decrease: Reduce standard deviations each iterations to slow
        down permutations
    :param bool remove_outlier: If registration hasn't converged, remove worst fitted
        spot when running particle filter
    :param bool debug: Print total distance in each iteration if true
    :return np.array t_matrix: Estimated 2D translation matrix
    :return float min_dist: Minimum total distance from fiducials to spots
    """
    # Pretrain spot coords
    dst = spot_coords.copy().astype(np.float32)
    knn = cv.ml.KNearest_create()
    labels = np.array(range(dst.shape[0])).astype(np.float32)
    knn.train(dst, cv.ml.ROW_SAMPLE, labels)

    nbr_particles = particles.shape[0]
    dists = np.zeros(nbr_particles)
    temp_stds = stds.copy()

    # Iterate until min dist doesn't change
    min_dist_old = 10 ** 6
    for i in range(max_iter):

        for p in range(nbr_particles):
            particle = particles[p]
            # Generate transformation matrix
            t_matrix = get_translation_matrix(particle)
            trans_coords = cv.transform(np.array([fiducial_coords]), t_matrix)
            trans_coords = trans_coords[0].astype(np.float32)
            # Find nearest spots
            ret, results, neighbors, dist = knn.findNearest(trans_coords, 1)
            if remove_outlier:
                # Remove worst fitted spot
                dist = dist[dist != np.amax(dist)]

            dists[p] = sum(dist)

        min_dist = np.min(dists)
        if debug:
            print(min_dist)
        # See if min dist is not decreasing anymore
        if abs(min_dist_old - min_dist) < stop_criteria:
            break
        min_dist_old = min_dist

        # Low distance should correspond to high probability
        weights = 1 / dists
        # Make weights sum to 1
        weights = weights / sum(weights)

        # Importance sampling
        idxs = np.random.choice(nbr_particles, nbr_particles, p=weights)
        particles = particles[idxs, :]

        # Reduce standard deviations a little every iteration
        temp_stds = temp_stds * iter_decrease ** i
        # Distort particles
        for c in range(4):
            distort = np.random.randn(nbr_particles)
            particles[:, c] = particles[:, c] + distort * temp_stds[c]

    # Return best particle
    particle = particles[dists == dists.min(), :][0]

    # Generate transformation matrix
    t_matrix = get_translation_matrix(particle)
    return t_matrix, min_dist
