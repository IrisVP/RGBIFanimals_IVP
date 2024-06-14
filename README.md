# RGBIFanimals_IVP
## Introduction
Species of maritime fauna all over the world are known to travel great distances in the oceans and seas. Some of these organisms are called ‘alien species’ when they travel from their natural habitat to another location that they didn’t inhabit before. With this workflow containing R scripts, we want to detect these alien species. This detection is realised by using occurrence data from GBIF and samplelocations to calculate sea distances between them.

## Origin
This is a revised workflow of RGBIFanimals: <br />
https://github.com/DvBMolEc/rgbifanimals

## General information
Usage of the scripts from RGBIFanimals original workflow is present in the user guide below. Example output can be found in the directory 'Example_Output'. Example input can be found in the directory 'Inputs'.

## Revised
In this Github page, the 2_Functions.R and 3_Main_script.R were revised and joined (+ adjusted) into one new big script called 3.2_New_Main_script.R. <br />
/!\ WARNING /!\ => The results wil not be correct because of the usage of a world map that is distorted when calculating sea distances!

## New workflow 'AlienDetective'
In the new workflow AlienDetective, there is a workflow of R scripts that does give good results. These scripts are also present in the directory 'new_workflow_Iris'. For more information on the new workflow, take a look at the AlienDetective Github page. <br />
https://github.com/IrisVP/AlienDetective

# Original README from RGBIFanimals
These R scripts will use BOLDigger output (or any other formatted sequence data containing taxonomic information), and determine the presence and diversity of alien species using GBIF occurrence data.

## User guide
1. Copy the first four R scripts in a new folder
2. Copy the folder _Inputs_ and its content in the just created folder
  - The _Inputs_ folder contains a tr.rdata file, which is used to calculate distances in 3_Main_script.R. Using this file saves some computing time.
3. Necessary input files have to be placed in the _Inputs_ filder and consist of
  - The _Boldigger_output.csv_, an example of which can be found in the repository
  - A _MetaData.csv_, of which the structure can be found in the example file
  - A _Synonyms.csv_ file, which contains taxonomic synonyms. This should be redundant as 3_Functions.R contains a function which checks the official name on WORMS, but making determined list can speed up processing.
4. Run the R scripts in order.
  - 1_Preperation.R will prepare the data for use in 3_Main_script.R and 4_Visualisation.R by creating several additional dataframes from the BOLDigger output
    All necessary files from 1_Preperation.R will be saved into a new folder _Output_.
  - 2_Functions.R contains all the functions used in 3_Main_script.R.
    - There is no need to actually run this script, everything will be run in the next script
  - 3_Main_script.R will calculate the shortest for each species found at each sampling location, to the nearest GBIF occurrence. 
    - Both the distance as the crow flies, as well as the distance over sea are calculated.
    - Depending on the amount if input data, the latter can take several hours of computing time.
  - 4_Visualisation.R will use the dataframes from 1_Preparation.R and the distances calculated in 3_Main_script.R to create several plots regarding native and alien species distribution, and diversity indices.
  - 5_Results.R is not necessary to run, but will provide insight in the ShortestPath versus DistanceOverSea differences
