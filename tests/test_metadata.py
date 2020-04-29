import pytest

from array_analyzer.extract import metadata
import array_analyzer.extract.constants as c

"""
Tests to add:
METADATA:
normal operation:
*requires:
    - input folder with good xlsx, good 
    - dummy output folder
1) extension .xml, .xlsx are found
2) no extension is found (not implemented error)
3) c.params are populated correctly for each of
    3a) .xml, .xlsx
4) arrays are populated correctly for each of 
    4a) .xml, .xlsx
5) fiducials are calculated correctly for each of
    5a) 6x6 array metadata, .xml and .xlsx
    5b) 6x8 array metadata, .xml and .xlsx
    
bad operation:
* requires:
    - input folder with bad xml:
        - multiple .xml
        - .xml with correct dictionaries
1) catch when multiple .xml
2) catch when incorrectly labeled .xlsx
3)      txtparser xlsx: catch when incorrect sheets
2) c.params has rows, columns, pitch, size defined
3) 
"""


def test_bad_xlsx(create_bad_xlsx, create_good_xlsx):
    bad = create_bad_xlsx
    good = create_good_xlsx
    c.METADATA_EXTENSION = 'xlsx'


    pass




import pytest
from testfixtures import TempDirectory
import pandas as pd


@pytest.fixture()
def create_bad_xlsx():
    # setup temp folder for fake .xlsx file
    temp = TempDirectory()
    temp_dir = temp.path

    # create dummy .xlsx
    # build dictionaries, use pandas to write .xlsx to temp
    # yield temp folder


    # breakdown temp folder
    #   delete temp and .xlsx

@pytest.fixture()
def create_good_xlsx():
    pass
