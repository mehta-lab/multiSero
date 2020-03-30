# serology-COVID19
Serological survey of antibody response to COVID-19 using ELISA and ELISA-arrays.

This open-source repository provides tools to analyze antibody responses from the data acquired using ELISA assays in [conventional multi-well format](https://doi.org/10.1101/2020.03.17.20037713) and [ELISA-array format](https://doi.org/10.1101/2019.12.20.885285).
The key goal is to enable rapid serological surveys of COVID-19 immunity. 

# Status
The code is still in infancy and is being rapidly developed. Present code is being written to analyze data from ELISA-arrays imaged with a variety of plate readers.
* [SciReader CL2](https://www.scienion.com/products/scireaders/).
* Open and configurable platform [Octopi](https://www.biorxiv.org/content/10.1101/684423v1) adapted for imaging multi-well plates.

The code is structured to be broadly useful for other serological analyses, and imagers.

# Structure

Current version is written to analyze data acquired using ELISA-arrays. The code is divided in two major parts:

**array_analyzer**:  python module to analyze the images from ELISA-arrays acquired with variety of plate readers.
* Inputs:
    * One image per well of a plate named by the well-index (`A1,A2,A3,...`).

    Well A1 (Flu Experiment)
    ![Well A1](https://drive.google.com/uc?export=view&id=1umSuvXz85mLqXzSVF_cm2QXCAFq0bXxL)
    
    Well E12 (Flu Experiment)
    ![Well E12](https://drive.google.com/uc?export=view&id=1umaFN5Ojrv_A0YGu3zVgKLCBJDhGpNiY)
    
    * .xml file that provides metadata for arrangement of antigens in 2D array. 
    
       [Link to .xml metadata for Flu Experiment](https://drive.google.com/file/d/1FoYHN28hAeBhkrGcikenEjfG9bzeZBMW/view?usp=sharing)

* Outputs.
    * Excel file (`OD...xlsx`): sheet named `antigens` shows arrangement of antigen spots, sheets named `A1,A2,A3,...` reports background corrected optical densities (ODs) at those spots.
    * Several debug plots to assess image analysis. 
    

**notebooks_interpretation**: collection of jupyter notebooks that show how to use output of `array_analyzer` to evaluate antibody binding. 

# Usage

the script "run_array_analyzer.py" can be run from command line if you comment out

```python
    main(input)
```

and uncomment:
```python
    main(sys.argv[1:])
```

then at cli, you can type:

```bash

    python run_array_imager.py -i <input_folder> -o <output_folder>

    (Optionally, you can add a flag "-d" for debug, which writes diagnostic images.)

```

This will look for .xml file in the "input_folder" (must be exactly 1) and grab all .png, .jpg, and .tiff images there.

Next it will extract the spots and create a subfolder for the specific processing run named "run_hour_min_sec".

Finally, within that run folder an excel workbook "OD.xlsx" is written, summarizing the Optical Density measurements of all spots
and all wells.  Individual spot and spot-composite images are written if -d debug mode is on.  This can be useful to see
how well the algorithm identified spots.

- <well_name>_crop.png
- <well_name>_crop_binary.png
- <well_name>_spot-1-2.png
- <well_name>_spot-2-2.png
- etc...

# Program workflow

1) Read array metadata
2) Identify the center of the well, crop it out
3) Identify the spots in the array.
4) Read each spot's region and convert that to an OD measurement
5) Generate a report of spots and ODs

##### 1) Read array metadata
- the .xml is read and parsed into separate dictionries
- contains print array info such as rows, columns, pitch, spot diameter
- contains locations of antigens, controls and fiducials

Example print run metadata
```xml
<layout rows="6" cols="8" vspace="0.41" hspace="0.405" expected_diameter="0.2" background_offset="0.05" background_thickness="0.05" max_diameter="0.3" min_diameter="0.1">
  <marker row="0" col="0" spot_type="Reference, Diagnostic" />
  <marker row="0" col="7" spot_type="Reference, Diagnostic" />
  <marker row="1" col="0" spot_type="Reference, Diagnostic" />
  <marker row="1" col="7" spot_type="Reference, PositiveControl" />
  <marker row="4" col="7" spot_type="Reference, NegativeControl" />
  <marker row="5" col="0" spot_type="Reference, Diagnostic" />
  <marker row="5" col="7" spot_type="Reference, Diagnostic" />
</layout>
<spots>
  <spot id="spot-1-2" row="0" col="1" spot_type="Diagnostic" />
  <spot id="spot-1-3" row="0" col="2" spot_type="Diagnostic" />
  .
  .
  .
```

##### 2) Identifying the center of the well
- as you can see from the images above, the well position can vary from image to image.
- we threshold to identify the boundary using OTSU

    ![Well A1 boundary](https://drive.google.com/uc?export=view&id=1uF7GiQRk0Agjrz3tiZ3Fls0sNkLNuTno)
    
    ![Well A1_cropped](https://drive.google.com/uc?export=view&id=1uMZLHPh3cFnQiJZLZIm2pvRPaYluEi5k)

##### 3) Identifying the spots in the array
- we use a unimodal threshold described by [Rosin](https://users.cs.cf.ac.uk/Paul.Rosin/resources/papers/unimodal2.pdf) 
to segment spots

    ![Well A1_cropped](https://drive.google.com/uc?export=view&id=1uUf775-vCuRmkxc0q52wr4Qgiofwr3C-)
    
- next we gather properties of each segmented region to apply some filters

- At this point, we assume the filters leave only real spots.  From these spots, we use the fiducials to interpolate 
a block onto each spot

    ![Well_A1_block_placed](https://drive.google.com/uc?export=view&id=1um6yui3yFX8BJQRCOxFAsqrL6IlLUytU)

The above steps can be repeated for estimated background levels

##### 4) Generating a report
- The information from each spot above is summarized as an Optical Density (OD) in the following images and report

    ![A1_OD](https://drive.google.com/uc?export=view&id=1ChcSAJeCkkT4PBBMezBOa60jwkbDqpyZ)
    
    Antigen worksheet of report
    ```text
        | 0            | 1     | 2     | 3     | 4      | 5     | 6     | 7    |
    |---|--------------|-------|-------|-------|--------|-------|-------|------| 
    | 0 |              | 114   | 100   | KZ52  | c13C6  | c2G4  | c4G7  |      | 
    | 1 |              | 114   | 100   | KZ52  | c13C6  | c2G4  | c4G7  |      | 
    | 2 | anti-HIS tag | 15731 | 15742 | 15750 | 15878  | 15946 | 15960 | Q411 | 
    | 3 | anti-HIS tag | 15731 | 15742 | 15750 | 15878  | 15946 | 15960 | Q411 | 
    | 4 | anti-HIS tag | 15974 | 16061 | FVM04 | VIC122 | Q206  | Q314  |      | 
    | 5 |              | 15974 | 16061 | FVM04 | VIC122 | Q206  | Q314  |      | 
    ```
    
    Optical Density worksheet of report
    ```text
        | 0           | 1           | 2           | 3           | 4           | 5           | 6           | 7           |
    |---|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------| 
    | 0 | 1.09348586  | 1.274936469 | 1.130797958 | 1.225837517 | 1.042772052 | 1.20056763  | 1.072387244 | 1.100760893 | 
    | 1 | 1.088949505 | 1.139443433 | 1.179673528 | 1.253505128 | 1.062039834 | 1.125543296 | 1.103111509 | 1.234257728 | 
    | 2 | 1.526547102 | 1.090505501 | 1.255221344 | 1.056073248 | 1.510661763 | 1.029962578 | 1.120706316 | 1.080969964 | 
    | 3 | 1.426604953 | 1.109832859 | 1.275447921 | 1.04740674  | 1.274042419 | 1.059029446 | 1.194087153 | 1.0730993   | 
    | 4 | 1.496396187 | 1.21556865  | 1.417112818 | 1.183898394 | 1.184774831 | 1.160777915 | 1.139684973 | 1.994962323 | 
    | 5 | 1.08574711  | 1.205482259 | 1.406633916 | 1.185217466 | 1.189582126 | 1.149181387 | 1.034703719 | 1.072368158 | 
    ```
    
    [link to report](https://drive.google.com/file/d/1usd1cVAJFzWANqR92SucaT6PpW581q8b/view?usp=sharing)
    
# Improvements
##### Iterative Closest Point fitting
- After segmenting the spots, we can fit a predicted grid array onto the segmentation via a rigid transform.
After multiple iterations the fit will naturally converge to filter out erroneous spots.

# Contributing

We welcome bug reports, feature requests, and contributions to the code. Please see issues on the repository for areas we need input on. 
The master branch is protected and meant to be always functional. Develop on fork of the repo and branches of this repo. Pull requests are welcome.
Please generate PRs after testing your code against real data and make sure that master branch is always functional.