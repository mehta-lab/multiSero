import numpy as np

import array_analyzer.transform.point_registration as registration

im_shape = (50, 75)
reg_coords = np.ones((5, 2), dtype=np.float32)


def test_create_reference_grid():
    grid_coords = registration.create_reference_grid((50, 100), 2, 4, 10)
    row_vals =[45, 55]
    col_vals = np.arange(85, 125, 10)
    count = 0
    for row in row_vals:
        for col in col_vals:
            assert grid_coords[count, 0] == row
            assert grid_coords[count, 1] == col
            count += 1


def test_check_reg_coords():
    reg_coords[2, :] = [24., 72.1]
    reg_ok = registration.check_reg_coords(reg_coords, im_shape, True)
    assert reg_ok is True


def test_check_reg_coords_large_col():
    reg_coords[2, :] = [24., 75.1]
    reg_ok = registration.check_reg_coords(reg_coords, im_shape, True)
    assert reg_ok is False


def test_check_reg_coords_large_row():
    reg_coords[2, :] = [50., 72.1]
    reg_ok = registration.check_reg_coords(reg_coords, im_shape, True)
    assert reg_ok is False


def test_check_reg_coords_negatve():
    reg_coords[2, :] = [-.1, 72.1]
    reg_ok = registration.check_reg_coords(reg_coords, im_shape, True)
    assert reg_ok is False


def test_check_reg_coords_already_false():
    reg_ok = registration.check_reg_coords(np.ones((5, 2)), (5, 5), False)
    assert reg_ok is False
