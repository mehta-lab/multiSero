import os
import pytest
import array_analyzer.utils.io_utils as io_utils


def test_read_gray_im(image_dir):
    im = io_utils.read_gray_im(os.path.join(image_dir, 'A1.png'))
    im_shape = im.shape
    assert im_shape[0] == 5
    assert im_shape[1] == 10


def test_no_im(image_dir):
    with pytest.raises(IOError):
        io_utils.read_gray_im(os.path.join(image_dir, 'no_im.png'))


def test_get_image_paths(image_dir):
    well_images = io_utils.get_image_paths(image_dir)
    assert len(well_images) == 4
    expected_keys = ['A1', 'A2', 'B11', 'B12']
    for key in expected_keys:
        im_name = os.path.basename(well_images[key])
        assert im_name == key + '.png'


def test_get_mm_image_paths(micromanager_dir):
    well_images = io_utils.get_image_paths(micromanager_dir)
    assert len(well_images) == 4
    expected_keys = ['A1', 'A2', 'B11', 'B12']
    for key in expected_keys:
        im_name = os.path.basename(well_images[key])
        assert im_name == 'micromanager_name.tif'
        well_name = well_images[key].split('/')[-2]
        # split again for well name, assume - separation
        well_name = well_name.split('-')[0]
        assert key == well_name


def test_make_run_dir(tmp_path):
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    input_dir = 'input_dir_name'
    run_dir = io_utils.make_run_dir(input_dir, output_dir)
    output_subdir = os.listdir(output_dir)
    assert len(output_subdir) == 1
    # Get last part of subdir and compare
    subdir_name = run_dir.split('/')[-1]
    assert subdir_name == output_subdir[0]
