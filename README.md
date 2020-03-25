# serology-COVID19
analysis of ELISA serological data 


### array_imager

the script "run_array_imager.py" can be run from command line if you comment out

```python
    path = '/Users/bryant.chhun/PycharmProjects/array-imager/Plates_given_to_manu/2020-01-15_plate4_AEP_Feb3_6mousesera'
    input = ['-i', path, '-o', path]
    main(input)
```

and uncomment:
```python
    main(sys.argv[1:])
```

then at cli, you can pass:

```bash

    python run_array_imager.py -i <input_folder> -o <output_folder>

```

This will look for .xml file in the "input_folder" (must be exactly 1) and grab all .png, .jpg, and .tiff images there.

Next it will extract the spots and create a subfolder for the specific processing run named "run_<hour>_<min>_<sec>".
Within that "run folder", a subfolder is created for for each image (well: A1, A2, B1, C2...).  This "well" subfolder contains
separate images for each spot like such:

- <well_name>_crop.png
- <well_name>_crop_binary.png
- <well_name>_spot-1-2.png
- <well_name>_spot-2-2.png
- etc...

Finally, the excel spreadsheet report is written to the "run folder".  Currently it contains only the summary and replicate summary.

It does not contain the fine detail parameters per-well (like Eric's reference data)