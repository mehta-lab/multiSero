"""
The constants.py namespace contains constants used throughout the repo
    These values are assigned at the start of analysis using the "metadata.MetaData" class
    Set the values by populating the .xml or .csv file as per documentation
"""

METADATA_EXTENSION = None

# constants parsed from metadata
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
}

layout = dict()
fiducials = dict()
spots = dict()
replicates = dict()

# array-constants
#   the constants below are all np.ndarrays whose elements are "U100" strings
SPOT_ID_ARRAY = None
SPOT_TYPE_ARRAY = None
FIDUCIAL_ARRAY = None
ANTIGEN_ARRAY = None

# constants needed for workflows
#   the constants below are lists of tuples, lists of integers, or simply integers
FIDUCIALS = []
FIDUCIALS_IDX = []
FIDUCIALS_IDX_8COLS = []
SCENION_SPOT_DIST = int()
STDS = []

# constants for saving
RUN_PATH = ''

# hardware parameters
