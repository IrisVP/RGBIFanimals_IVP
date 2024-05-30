
############################################################################################
# LOAD PACKAGES
############################################################################################
library("gdistance")
library("dplyr")
require("geosphere")
require("rgbif")
library("ggplot2")
library("tidyr")
library("worrms")
library("sf")
library("maps")
library("FRK")

species_name <- "Acartia bifilosa"

for (species_name in unique(long$Specieslist)){
  print(paste0("species_name: ", species_name))
  # if file exists:        put content into variable res
  # if file doesn't exist: make one, get data, and put into variable res
  occurrence_coord <- paste0("OccurrenceData_test_more/", species_name, ".csv")
  if (file.exists(occurrence_coord) == TRUE) {
    res <- read.csv(occurrence_coord, header = TRUE)
  } else {
    res1 <- occ_data(scientificName = species_name, 
                     hasCoordinate = TRUE,
                     limit = 10000,
                     basisOfRecord = "Observation",
    )
    res2 <- occ_data(scientificName = species_name, 
                     hasCoordinate = TRUE, 
                     limit = 10000,
                     basisOfRecord = "Machine observation",
    )
    res3 <- occ_data(scientificName = species_name, 
                     hasCoordinate = TRUE, 
                     limit = 10000,
                     basisOfRecord = "Human observation",
    )
    res4 <- occ_data(scientificName = species_name, 
                     hasCoordinate = TRUE, 
                     limit = 10000,
                     basisOfRecord = "Material sample",
    )
    res5 <- occ_data(scientificName = species_name, 
                     hasCoordinate = TRUE, 
                     limit = 10000,
                     basisOfRecord = "Material citation",
    )
    res6 <- occ_data(scientificName = species_name, 
                     hasCoordinate = TRUE, 
                     limit = 10000,
                     basisOfRecord = "Preserved specimen",
    )
    res7 <- occ_data(scientificName = species_name, 
                     hasCoordinate = TRUE, 
                     limit = 10000,
                     basisOfRecord = "Fossil specimen",
    )
    res8 <- occ_data(scientificName = species_name, 
                     hasCoordinate = TRUE, 
                     limit = 10000,
                     basisOfRecord = "Living specimen",
    )
    res9 <- occ_data(scientificName = species_name, 
                     hasCoordinate = TRUE, 
                     limit = 10000,
                     basisOfRecord = "Occurrence",
    )
    
    res1_coord <- res1$data[, c('decimalLongitude', 'decimalLatitude')]
    res2_coord <- res2$data[, c('decimalLongitude', 'decimalLatitude')]
    res3_coord <- res3$data[, c('decimalLongitude', 'decimalLatitude')]
    res4_coord <- res4$data[, c('decimalLongitude', 'decimalLatitude')]
    res5_coord <- res5$data[, c('decimalLongitude', 'decimalLatitude')]
    res6_coord <- res6$data[, c('decimalLongitude', 'decimalLatitude')]
    res7_coord <- res7$data[, c('decimalLongitude', 'decimalLatitude')]
    res8_coord <- res8$data[, c('decimalLongitude', 'decimalLatitude')]
    res9_coord <- res9$data[, c('decimalLongitude', 'decimalLatitude')]
    
    res_total <- rbind(res1_coord, res2_coord, res3_coord, res4_coord, res5_coord, res6_coord,
                       res7_coord, res8_coord, res9_coord)
    
    #rename the column names
    colnames(res_total) <- c('Longitude', 'Latitude')
    # Remove occurrences where longitude or latitude is NA
    res_total <- res_total[!is.na(res_total$Latitude) & !is.na(res_total$Longitude),]
    # check if there's no information for a species
    if (nrow(res_total) == 0) {
      error_message <- paste0("No information found on GBIF for ", species_name)
      
      # check if directory with error messages exists, if it doesn't: make one
      if (!dir.exists("test_outputs/errors")){
        dir.create("test_outputs/errors", recursive = TRUE)
        error_file_name <- paste0("test_outputs/errors/error_", species_name, ".csv")
        # error files written to test_outputs/errors/
        write.csv(error_message, file = error_file_name)
        return(FALSE)
        
        # write error file to the directory
      } else {
        error_file_name <- paste0("test_outputs/errors/error_", species_name, ".csv")
        # error files written to test_outputs/errors/
        write.csv(error_message, file = error_file_name)
        return(FALSE)
      }
    }
  }
  print("file has successfully been written")
  write.csv(res_total, occurrence_coord)
}
