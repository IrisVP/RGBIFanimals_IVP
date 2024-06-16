### Short description of R scripts

For detailed documentation: see the 'Documentation' directory for .Rmd files

#### 1_Calculation_sea_distances.R

1. Download and install libraries, then load data.
2. Retrieve occurrence data from GBIF, focusing on Europe. The settings for fetching data can be personalized here.
3. Calculate distances: both flying and sea distances are calculated and saved in files inside the 'Output_calculations' directory.

#### 2_Visualisation_distribution.R

Histograms are created showing the frequencies of species occurrences per sample location and species. These graphs are saved in the 'Output_calculations' directory, with a separate directory for each type of graph:

1. The first graph shows the sea distances per species and sample location.
2. The second graph is a combined graph of flying distances and sea distances.
3. The third graph is similar to the first graph but colored by country.
4. The fourth graph is also similar to the first graph but colored by year category.

Patterns can be observed in these histograms. The x-axis represents distances in km, and the y-axis represents frequency. An alien species can be detected when there is a spike in frequency at a large distance. Using the country plot, the origin region of the species can be identified. Using the year plot, the year in which the species migrated to that location can be determined.
