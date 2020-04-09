"""
The constants.py namespace contains constants used throughout the repo

"""

# constants parsed from .xml file
params = {
    'rows': None,
    'columns': None,
    'v_pitch': None,
    'h_pitch': None,
    'spot_width': None,
    'bg_offset': None,
    'bg_thickness': None,
    'max_diam': None,
    'min_diam': None,
    'pixel_size_scienion': 0.0049,
    'pixel_size_octopi': 0.00185,
}

layout = dict()
fiducials = dict()
spots = dict()
replicates = dict()

# array-constants
SPOT_ID_ARRAY = None
SPOT_TYPE_ARRAY = None
ANTIGEN_ARRAY = None
FIDUCIAL_ARRAY = None

# constants needed for workflows

FIDUCIALS = []
FIDUCIALS_IDX = []
FIDUCIALS_IDX_8COLS = []
SCENION_SPOT_DIST = []
STDS = []