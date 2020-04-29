"""
The constants.py namespace contains constants used throughout the repo
    These values are assigned at the start of analysis using the "metadata.MetaData" class
    Set the values by populating the .xml or .csv file as per documentation
"""

METADATA_EXTENSION = None

# === constants parsed from metadata ===
#   the constants below are all dictionaries
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
    'pixel_size': None
}

layout = dict()
fiducials = dict()
spots = dict()
replicates = dict()

# === array-constants ===
#   the constants below are all np.ndarrays whose elements are "U100" strings
SPOT_ID_ARRAY = None
SPOT_TYPE_ARRAY = None
FIDUCIAL_ARRAY = None
ANTIGEN_ARRAY = None

# === constants needed for workflows ===
#   these values are extracted from the above ARRAYs

# list of (column, row) fiducial locations
FIDUCIALS = []

# list of int fiducial location in array.  Example for 6x6 array: (0,0) = 0, (0,1) = 6, (5,5) = 35
FIDUCIALS_IDX = []

# values used by point_registration.py
SPOT_DIST_PIX = int()
SPOT_DIST_UM = int()
STDS = [500, 500, .1, .001]  # x, y, angle, scale
REG_DIST_THRESH = 1000

# constants for saving
RUN_PATH = ''
