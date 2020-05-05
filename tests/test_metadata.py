import pytest

from array_analyzer.extract.metadata import MetaData
import array_analyzer.extract.constants as constants
import os

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


def test_one_xml(create_good_xml):
    constants.METADATA_EXTENSION = 'xml'

    output_dir = create_good_xml
    MetaData(output_dir, output_dir)


def test_two_xml(create_good_xml):
    constants.METADATA_EXTENSION = 'xml'
    output_dir = create_good_xml

    with open(os.path.join(output_dir, 'second_xml.xml'), 'w') as fp:
        pass
    with pytest.raises(IOError):
        MetaData(output_dir, output_dir)


def test_xlsx(create_good_xlsx):
    constants.METADATA_EXTENSION = 'xlsx'
    output_dir = create_good_xlsx
    MetaData(output_dir, output_dir)

