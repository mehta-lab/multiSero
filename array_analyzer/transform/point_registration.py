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


def create_gaussian_particles(mean_point,
                              stds,
                              scale_mean=1,
                              angle_mean=0,
                              nbr_particles=100):
    """
    Create particles from parameters x, y, scale and angle given mean and std.
    :param mean_point:
    :param stds:
    :param scale_mean:
    :param angle_mean:
    :param nbr_particles:
    :return:
    """
    # x, y, scale, angle
    particles = np.empty((nbr_particles, 4))
    particles[:, 0] = mean_point[0] + (np.random.randn(nbr_particles) * stds[0])
    particles[:, 1] = mean_point[1] + (np.random.randn(nbr_particles) * stds[1])
    particles[:, 2] = angle_mean + (np.random.randn(nbr_particles) * stds[2])
    particles[:, 3] = scale_mean + (np.random.randn(nbr_particles) * stds[3])
    return particles


def get_translation_matrix(particle):

    a = particle[3] * np.cos(particle[2] * np.pi / 180)
    b = particle[3] * np.sin(particle[2] * np.pi / 180)
    t_matrix = np.array([[a, b, particle[0]],
                        [-b, a, particle[1]]])
    return t_matrix


def particle_filter(fiducial_coords,
                    spot_coords,
                    particles,
                    stds,
                    max_iter=50,
                    stop_criteria=.01):
    """
    Particle filtering to determine best grid location
    :param fiducial_coords:
    :param spot_coords:
    :param particles:
    :param stds:
    :param max_iter:
    :param stop_criteria:
    :return:
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
        # Reduce standard deviations a little every iteration
        temp_stds = temp_stds * 0.8 ** i

        for p in range(nbr_particles):
            particle = particles[p]
            # Generate transformation matrix
            t_matrix = get_translation_matrix(particle)
            trans_coords = cv.transform(np.array([fiducial_coords]), t_matrix)
            trans_coords = trans_coords[0].astype(np.float32)

            # Find nearest spots
            ret, results, neighbors, dist = knn.findNearest(trans_coords, 1)
            dists[p] = sum(dist)  # np.linalg.norm(dist)

        min_dist = np.min(dists)
        # print(min_dist)
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

        # Distort particles
        for c in range(4):
            distort = np.random.randn(nbr_particles)
            particles[:, c] = particles[:, c] + distort * temp_stds[c]


    # Return best particle
    particle = particles[dists == dists.min(), :][0]
    # Generate transformation matrix
    t_matrix = get_translation_matrix(particle)
    return t_matrix
