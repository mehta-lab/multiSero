![unittest](https://github.com/czbiohub/serology-COVID19/workflows/unittest/badge.svg)

# pysero

pysero enables serological measurements with multiplexed and standard ELISA assays.

The project automates measurement of antibody response from data collected with ELISA assays in [conventional multi-well format](https://doi.org/10.1101/2020.03.17.20037713) and [ELISA-array format](https://doi.org/10.1101/2019.12.20.885285).

## Features

The immediate goal is to enable rapid serological surveys for COVID-19. 

The project aims to implement serological analysis for all antigen multiplexing technologies that emerge to be useful: 
* classical ELISA.
* antigen arrays printed with [Scienion](https://www.scienion.com/products/sciflexarrayers/).
* antigen arrays printed with [Echo](https://www.labcyte.com/echo-liquid-handling).
* antigen multiplexing with [Luminex](https://www.luminexcorp.com/blog/multiplex-technologies-more-effective-than-elisa-for-antibody-detection/) beads. 

The image data can be acquired with diversity of imagers and plate readers:
 * turn-key plate imagers, e.g., [SciReader CL2](https://www.scienion.com/products/scireaders/).
 * [Octopi](https://www.biorxiv.org/content/10.1101/684423v1) platform with online analysis - in development in collaboration with [Prakash Lab](http://web.stanford.edu/group/prakash-lab/cgi-bin/labsite/).
 * any widefield microscope with motorized XY stage.

The project will also have tools for combining data from different types of assays for interpretation of antibody binding.

## Overview of pipeline

<img src="docs/Workflow%20Schematic.png" width="600">

## Status
The code is being rapidly developed. 
Current code is validated for analysis of ELISA-arrays imaged with Scienion reader and is being refined for antigen arrays imaged with Octopi.


## Usage
The script "run_array_analyzer.py" can be run from command line

```buildoutcfg
python run_array_analyzer.py --input <input dir> --output <output dir> --method <'interp' or 'fit'> --debug
```

This will look for .xml file in the input directory (must be exactly 1) and grab all .png, .jpg, and .tiff images there.

Next it will extract the spots and create a subfolder for the specific processing run named "run_hour_min_sec".

Finally, within that run folder an excel workbook "OD.xlsx" is written, summarizing the Optical Density measurements of all spots
and all wells.  Individual spot and spot-composite images are written if -d debug mode is on.  This can be useful to see
how well the algorithm identified spots.

- <well_name>_crop.png
- <well_name>_crop_binary.png
- <well_name>_spot-1-2.png
- <well_name>_spot-2-2.png
- etc...

Collection of jupyter notebooks in 
**notebooks_interpretation** show how to use ODs to evaluate antibody binding. 

## Workflow
Steps in the analysis workflow are described in [workflow](docs/workflow.md)

## Contributions
We welcome bug reports, feature requests, and contributions to the code. Please see  [this](docs/CONTRIBUTING.md) page for most fruitful ways to contribute.

