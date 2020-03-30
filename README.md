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
**array_analyzer**: python module to analyze the images from ELISA-arrays acquired with variety of plate readers.
* Inputs:
    * One image per well of a plate named by the well-index (`A1,A2,A3,...`).
    #TODO: add few images and their names.
    * .xml file that provides metadata for arrangement of antigens in 2D array. 
    # TODO: link to an xml in the repo or google-drive.

* Outputs.
    * Excel file (`OD...xlsx`): sheet named `antigens` shows arrangement of antigen spots, sheets named `A1,A2,A3,...` reports background corrected optical densities (ODs) at those spots.
    * Several debug plots to assess image analysis. 

**notebooks_interpretation**: collection of jupyter notebooks that show how to use output of `array_analyzer` to evaluate antibody binding. 

# Usage

## array_analyzer

the script "run_array_analyzer.py" can be run from command line if you comment out

```python
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
#TODO: describe what is currently output. 

# Contributing

We welcome bug reports, feature requests, and contributions to the code. Please see issues on the repository for areas we need input on. 
The master branch is protected and meant to be always functional. Develop on fork of the repo and branches of this repo. Pull requests are welcome.
Please generate PRs after testing your code against real data and make sure that master branch is always functional.

## array_analyzer

#TODO: Following is the transformation of data and links to current implementation (date: ):
# [segmentation of well](permalink to current master)
# [identification of spot](permalink to current master)
# ...
