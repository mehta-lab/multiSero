
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
