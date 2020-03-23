# bchhun, {2020-03-23}

"""
Here we will call the methods from the ETL folders

todo: include CLI capabilities with flags

"""

"""
txt_parser workflow:
--------------------
1) xml_to_dict the xml file
2) create ID-array
3) create antigen-array

image_parser workflow:
----------------------
4) read_to_grey(supplied images)

# find center of well
5) thresh and binarize from 1
6) find well border from 2
7) crop image from 3

# find center of spots from crop
8) thresh and binarize from 4
9) clean spot binary from 5
10) generate props from 6
11) generate props dict from 7
12) assign props dict to array from 8

xlsx report generation workflow:
--------------------------------
13) "create base template"
14) "populate main tab" using :
        - workbook from 13
        - "ID-array" from 2
        - "props-array" from 12
        - "well" from "read_to_grey" from 4
15) *repeat 13-14* calling next image and well name

"""