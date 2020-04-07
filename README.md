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
    ![Well A1](https://drive.google.com/uc?export=view&id=1utiSZF_jnIDFAuDYZ2TvZS7BjwmBqOQh)

    Well E12 (Flu Experiment)
    ![Well E12](https://drive.google.com/uc?export=view&id=1uwtxcpIDsBDwET7IEvcdjwYn4Uxz4_mf)
    
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

    python run_array_imager.py -i <input_folder> -o <output_folder> -m <method>

    (where <method> is one of "fit" or "interp")
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



