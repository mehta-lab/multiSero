# serology-COVID19
analysis of ELISA serological data 


## array_analyzer

the script "run_array_analyzer.py" can be run from command line if you comment out

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

### performance of array_imager

Using %timeit in the jupyter notebook, I find:

```python
    %timeit edges = canny(binary, sigma=3)
    932 ms ± 37.1 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    
    %timeit hough_res = hough_circle(edges, hough_radii)
    643 ms ± 21.3 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    
    %timeit aaccums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=1)
    777 ms ± 12.5 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

```
While the processing reports (at runtime) total well processing time as such:

```bash
    A6.png
        time to process=2.6693918704986572
    A7.png
        time to process=2.677278995513916
    A5.png
        time to process=2.785794973373413
    B9.png
        time to process=2.57763409614563
    B8.png
        time to process=2.824164867401123
    A4.png
        time to process=2.5805811882019043

```
From this, one can see how nearly 2.3 of the 2.5 seconds are taken by the "canny" and "hough" transforms

### improvements, todos:

- if the Rosin thresholding cuts out a spot, the spot will be “None” in the array (instead of a regionprop).  In this case, the summary in the Excel sheet is also empty.
We can solve this by “guessing” the well position based on neighboring regionprop.bbox values.

- it looks like the summary stats depend a lot on the rosin-thresholded image (as in, the stats are based on pixels inside the mask).  We can solve this by re-labeling the image based on the regionprop.bbox, then feeding this again through regionprop.

- we can explore replacements to canny (which is probably overkill), and maybe even Hough circle that may speed this up