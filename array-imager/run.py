# bchhun, {2020-03-23}

"""
Here we will call the methods from the ETL folders

============================
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
15) "populate main replictes" :
        - workbook from 13
        - "props-array" from 12
        - "antigen-array" from 3
        - "well" from "read_to_grey" from 4
16) (populate well-tab)
17) (populate well-replicate-tab)

18) *repeat 14-17* using next image and well name
19) save .xlsx
==============================

==============================
FULL WORKFLOW

cli
---
- input folder
- output folder

extract
-------
A) search folder for all images, all .xmls to list
B) xlsx_report.create_base_template() step 13

C) txt_parse workflow above to create ID-array, antigen-array
D) image_parser workflow above to loop 4 (read_to_grey)
    within loop:
    E) image_parser steps 5-12

    transform
    ---------
    F) (ANY "transform" methods that will further analyze the data from E)
        (this is set aside as a place to make diagnosis calls and for downstream calculations using spot properties)

    load
    ----
    G) xlsx_report generation workflow steps 14-17

"""

