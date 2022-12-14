import numpy as np
import itertools


class BackgroundEstimator2D:
    """Estimates flat field image"""

    def __init__(self,
                 block_size=128,
                 order=2,
                 normalize=True):
        """
        Background images are estimated once per channel for 2D data
        :param int block_size: Size of blocks image will be divided into
        :param int order: Order of polynomial (default 2)
        :param bool normalize: Normalize surface by dividing by its mean
            for background correction (default True)
        """

        if block_size is None:
            block_size = 128
        self.block_size = block_size
        self.order = order
        self.normalize = normalize

    def sample_block_medians(self, im):
        """Subdivide a 2D image in smaller blocks of size block_size and
        compute the median intensity value for each block. Any incomplete
        blocks (remainders of modulo operation) will be ignored.

        :param np.array im:         2D image
        :return np.array(float) sample_coords: Image coordinates for block
                                               centers
        :return np.array(float) sample_values: Median intensity values for
                                               blocks
        """
        im_shape = im.shape
        assert self.block_size < im_shape[0], "Block size larger than image height"
        assert self.block_size < im_shape[1], "Block size larger than image width"

        nbr_blocks_x = im_shape[0] // self.block_size
        nbr_blocks_y = im_shape[1] // self.block_size
        sample_coords = np.zeros((nbr_blocks_x * nbr_blocks_y, 2),
                                 dtype=np.float64)
        sample_values = np.zeros((nbr_blocks_x * nbr_blocks_y, ),
                                 dtype=np.float64)
        for x in range(nbr_blocks_x):
            for y in range(nbr_blocks_y):
                idx = y * nbr_blocks_x + x
                sample_coords[idx, :] = [x * self.block_size + (self.block_size - 1) / 2,
                                         y * self.block_size + (self.block_size - 1) / 2]
                sample_values[idx] = np.median(
                    im[x * self.block_size:(x + 1) * self.block_size,
                       y * self.block_size:(y + 1) * self.block_size]
                )
        return sample_coords, sample_values

    def fit_polynomial_surface_2d(self,
                                  sample_coords,
                                  sample_values,
                                  im_shape):
        """
        Given coordinates and corresponding values, this function will fit a
        2D polynomial of given order, then create a surface of given shape.

        :param np.array sample_coords: 2D sample coords (nbr of points, 2)
        :param np.array sample_values: Corresponding intensity values (nbr points,)
        :param tuple im_shape:         Shape of desired output surface (height, width)

        :return np.array poly_surface: 2D surface of shape im_shape
        """
        assert (self.order + 1)*(self.order + 2)/2 <= len(sample_values), \
            "Can't fit a higher degree polynomial than there are sampled values"
        # Number of coefficients is determined by (order + 1)*(order + 2)/2
        orders = np.arange(self.order + 1)
        variable_matrix = np.zeros(
            (sample_coords.shape[0], int((self.order + 1)*(self.order + 2)/2)),
        )
        order_pairs = list(itertools.product(orders, orders))
        # sum of orders of x,y <= order of the polynomial
        variable_iterator = itertools.filterfalse(lambda x: sum(x) > self.order, order_pairs)
        for idx, (m, n) in enumerate(variable_iterator):
            variable_matrix[:, idx] = sample_coords[:, 0] ** n * sample_coords[:, 1] ** m
        # Least squares fit of the points to the polynomial
        coeffs, _, _, _ = np.linalg.lstsq(variable_matrix, sample_values, rcond=-1)
        # Create a grid of image (x, y) coordinates
        x_mesh, y_mesh = np.meshgrid(np.linspace(0, im_shape[1] - 1, im_shape[1]),
                                     np.linspace(0, im_shape[0] - 1, im_shape[0]))
        # Reconstruct the surface from the coefficients
        poly_surface = np.zeros(im_shape, np.float)
        order_pairs = list(itertools.product(orders, orders))
        # sum of orders of x,y <= order of the polynomial
        variable_iterator = itertools.filterfalse(lambda x: sum(x) > self.order, order_pairs)
        for coeff, (m, n) in zip(coeffs, variable_iterator):
            poly_surface += coeff * x_mesh ** m * y_mesh ** n

        return poly_surface

    def get_background(self, im):
        """
        Combine sampling and polynomial surface fit for background estimation.
        To background correct an image, divide it by background.

        :param np.array im: 2D grayscale image
        :return np.array background: Background image
        """
        # Get grid of coordinates with median intensity values
        coords, values = self.sample_block_medians(im=im)
        # Estimate background from grid
        background = self.fit_polynomial_surface_2d(
            sample_coords=coords,
            sample_values=values,
            im_shape=im.shape,
        )
        # Normalize by mean
        if self.normalize:
            background /= np.mean(background)

        return background
