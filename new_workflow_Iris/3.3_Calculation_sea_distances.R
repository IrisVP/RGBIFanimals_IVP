
##################################################################################
# This script is an alternative script for calculating sea distances and filtering
# on alien species
##################################################################################
# OUTPUTS:
# The DistanceOverSea.csv file is created in test_outputs/ directory instead of Output/
# The Occurrence data files are created in OccurrenceData_test/ instead of OccurrenceData/
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
library("sp")
library("rgdal")
library("raster")
library("rnaturalearth")
library("rnaturalearthdata")
library("maps")
library("maptools")
library("FRK")


############################################################################################
# LOAD DATA
############################################################################################
# Set working directory to directory where the R-script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Read a species list
df <- read.csv("Output/RealFirst10_Species_Location.csv")
# Read coordinates file
Coordinates <- read.csv("Inputs/Coordinates.csv")


############################################################################################
# RESHAPE DF WITH LOCATIONS AND SPECIES FROM WIDE TO LONG FORMAT
############################################################################################

long <- pivot_longer(df, !Specieslist)
# from tidyr package, reshape df from wide to long format
long <- long[long$value > 0, ]

#####################################################################################
# FIRST CHECK IF FILE WITH OCCURRENCE DATA IN OCCURRENCEDATA DIRECTORY EXISTS
# IF NOT: MAKE ONE AND GET DATA FROM GBIF
# limit in occ_data() function is changeable for personal preference
####################################################################################

# practice species name
#species_name <- "Acartia bifilosa"

# Function to ensure all required columns are present (used in coming for loop)
ensure_columns <- function(df, required_columns) {
  missing_columns <- setdiff(required_columns, colnames(df))
  for (col in missing_columns) {
    if (!col %in% colnames(df)) {
      df[[col]] <- NA
    }
  }
  return(df)
}
fetch_data_in_batches <- function(species_name, basisOfRecord, batch_size = 10000) {
  start <- 0
  combined_data <- data.frame()
    res <- occ_data(scientificName = species_name, 
                    hasCoordinate = TRUE, 
                    limit = batch_size,
                    basisOfRecord = basisOfRecord,
                    continent = "europe")
  combined_data <- rbind(combined_data, res$data)
  Sys.sleep(1)  # Adding a small delay to be polite to the server
  
  return(combined_data)
}

required_columns <- c('decimalLongitude', 'decimalLatitude', 'year', 'month', 'country')

for (species_name in unique(long$Specieslist)){
  print(paste0("species_name: ", species_name))
  # if file exists:        put content into variable res
  # if file doesn't exist: make one, get data, and put into variable res
  occurrence_coord <- paste0("OccurrenceData_test/", species_name, ".csv")
  #if (file.exists(occurrence_coord) == TRUE) {
  #  res <- read.csv(occurrence_coord, header = TRUE)
  #} else {
  data_list <- list(
    fetch_data_in_batches(species_name, "Observation"),
    fetch_data_in_batches(species_name, "Machine observation"),
    fetch_data_in_batches(species_name, "Human observation"),
    fetch_data_in_batches(species_name, "Material sample"),
    fetch_data_in_batches(species_name, "Material citation"),
    fetch_data_in_batches(species_name, "Living specimen"),
    fetch_data_in_batches(species_name, "Occurrence")
  )

  # Initialize an empty list to store processed data frames
  processed_data <- list()
  
  # Loop over each data frame in the list, ensure columns, and select required columns
  for (i in seq_along(data_list)) {
    temp_df <- data_list[[i]]
    
    if (!is.null(temp_df) && nrow(temp_df) > 0) { # check if temp_df is not NULL and not empty
      #print(paste0("Summary for data frame ", i, " before ensuring columns:"))
      #print(summary(temp_df))
      
      # Ensure the required columns are present
      temp_df <- ensure_columns(temp_df, required_columns)
      temp_df <- temp_df[, required_columns]
        
      processed_data[[i]] <- temp_df
      #print(paste0("Summary for data frame ", i, " after ensuring columns:"))
      #print(summary(processed_data[[i]]))
    } else {
      print(paste0("Data frame ", i, " is empty."))
    }
  }
    
  if (length(processed_data) > 0) {
    res_total <- do.call(rbind, processed_data)
    print(paste0("Total NA values in res_total: ", sum(is.na(res_total))))
    #print("Summary of res_total:")
    #print(summary(res_total))
  } else {
    print(paste0("No data to combine for species: ", species_name))
  }
  #rename the column names
  colnames(res_total) <- c('Longitude', 'Latitude', 'year', 'month', 'country')
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
  print("file has successfully been written")
  write.csv(res_total, occurrence_coord)
}


############################################################################################
# REVISION: CALCULATE DISTANCES
############################################################################################
# practice species
#species_name <- "Acartia bifilosa"
#species_location <- "Koster"

# Iterate over species_name and location_name
Calculation_seadistance <- function(species_name, species_location){
  
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  
  occurrence_coord <- paste0("OccurrenceData_test/", species_name, ".csv")
  OccurrenceData <- read.csv(occurrence_coord, header = TRUE)
  
  # getting coordinates from Coordinates dataframe to get longitude and latitude for ARMS location
  ################################################################################################
  
  # use grep to get the row from the Coordinates df where location_name is present
  location_row_index <- grep(species_location, Coordinates$Observatory.ID)
  # save longitude and latitude for the row that you selected with grep
  longitude <- Coordinates[location_row_index, "Longitude"]
  latitude <- Coordinates[location_row_index, "Latitude"]
  # make a dataframe out of the longitude and latitude called samplelocation
  samplelocation <- data.frame(Latitude = latitude, Longitude = longitude)
  
  # check if OccurrenceData has coordinates
  if(nrow(OccurrenceData) < 1) {
    stop("OccurrenceData has no coordinates (line 139)")
  }
  
  ##########################################################################
  # DISTANCES CALCULATION
  ##########################################################################

  # Load world data
  world <- ne_countries(scale = "medium", returnclass = "sf")
  r <- raster(extent(-180, 180, -90, 90), res = 0.1)
  r <- rasterize(world, r, field = 1, fun = max, na.rm = TRUE)
  costs <- reclassify(r, cbind(1, Inf))
  costs[is.na(costs)] <- 1  # Set water cells to a low cost (e.g., 1)
  
  # for loop to iterate over OccurrenceData
  # start a sea_distances and flying_distances list
  sea_distances <- c()
  flying_distances <- c()
  
  for (row in 1:nrow(OccurrenceData)) {
    print(paste0("Calculating longitude: ", OccurrenceData[row, 3], " and latitude: ", OccurrenceData[row, 2]))
    print(paste0("for samplelocation: ", samplelocation[,1], " ", samplelocation[,2]))
    
    # Define points using correct projection
    point1 <- SpatialPoints(cbind(samplelocation$Longitude, samplelocation$Latitude), proj4string = CRS(proj4string(r)))
    point2 <- SpatialPoints(cbind(OccurrenceData[row, 3], OccurrenceData[row, 2]), proj4string = CRS(proj4string(r)))
  
    #####################
    ## FLYING DISTANCE ##
    #####################
  
    # Calculate the straight-line distance (accounting for the Earth's curvature)
    straight_line_distance <- distHaversine(coordinates(point1), coordinates(point2))
    print(paste("Straight line distance:", straight_line_distance, "meters"))
    
    flying_distances <- append(flying_distances, straight_line_distance)
  
    #################
    ## SEA DISTANCE##
    #################
  
    transition_matrix <- "transitMatrix.rds"
    if (!file.exists(transition_matrix)) {
      # Create a transition object for adjacent cells
      transitMatrix <- transition(costs, transitionFunction = function(x) 1/mean(x), directions = 16)
      # Set infinite costs to NA to prevent travel through these cells
      transitMatrix <- geoCorrection(transitMatrix, scl = TRUE)
      # Save/Load transition matrix (this should be the same for all calculations)
      saveRDS(transitMatrix, file = "transitMatrix.rds")
      
    } else {
      transitMatrix <- readRDS(file = "transitMatrix.rds")
    }
    
    # Coerce points to SpatialPointsDataFrame for compatibility with gdistance
    point1_df <- SpatialPointsDataFrame(coords = point1, data = data.frame(id = 1), proj4string = CRS(proj4string(r)))
    point2_df <- SpatialPointsDataFrame(coords = point2, data = data.frame(id = 2), proj4string = CRS(proj4string(r)))
  
    # Compute cost distance
    cost_distance <- costDistance(transitMatrix, point1_df, point2_df)
  
    # Calculate the shortest path
    shortest_path <- shortestPath(transitMatrix, point1_df, point2_df, output = "SpatialLines")
  
    # Plotting the shortest path and the world map
    #plot(r, main = "Shortest Water Path")
    #plot(world, add = TRUE, col = "grey")
    #plot(shortest_path, add = TRUE, col = "blue", lwd = 2)
    #points(point1_df, col = "red", pch = 20)
    #points(point2_df, col = "green", pch = 20)
  
    # Assuming 'shortest_path' is your SpatialLines object from the shortestPath function
    # First, ensure the CRS is set on the original SpatialLines object
    crs_info <- proj4string(shortest_path)  # or use crs(shortest_path) if using `sp`
  
    # If it's not set, set it here, assuming the original data was in WGS 84 (EPSG:4326)
    if (is.na(crs_info)) {
      proj4string(shortest_path) <- CRS("+init=epsg:4326")  # Set to WGS 84
    }
  
    # Convert SpatialLines to sf object
    shortest_path_sf <- st_as_sf(shortest_path)
  
    # Confirm CRS is set for sf object, if not, set it:
    if (is.na(st_crs(shortest_path_sf))) {
      st_crs(shortest_path_sf) <- 4326  # EPSG code for WGS 84
    }
  
    # Transform to a suitable projected CRS for distance calculation (e.g., UTM zone 33N)
    shortest_path_utm <- st_transform(shortest_path_sf, 32633)  # UTM zone 33N
  
    # Calculate the length in meters
    path_length <- st_length(shortest_path_utm)
  
    # Print the length
    print(paste0("distance through sea: ", path_length))
    sea_distances <- append(sea_distances, path_length)
  }
  
  print(paste0("sea_distances: ", sea_distances))
  print(paste0("flying distances: ", flying_distances))
}
  # check if directory and/or csv files already exists with sea distances(files for each species/location)
  distance_file <- paste0("theoretical_data/sea_distances/", species_name, "_distancesTo_", species_location)
  content_file <- cbind(sea_distances, OccurrenceData$year, OccurrenceData$month, OccurrenceData$country)
  colnames(content_file) <- c("x", "year", "month", "country")
  
  if(!dir.exists("theoretical_data/sea_distances")){
    dir.create("theoretical_data/sea_distances")
    write.table(content_file, file = distance_file)
    colnames(content_file)
  } else {
    write.table(content_file, file = distance_file)
  }
  
  ###########################################################
  # CALCULATE SHORTEST SEA DISTANCE
  ###########################################################
  
  #shortest_fly_distance <- min(flying_distances, na.rm = TRUE)
  #print(shortest_fly_distance)
  #shortest_sea_distance <- min(sea_distances, na.rm = TRUE)
  #print(shortest_sea_distance)
  
  # make a list of all shortest distances per species/location
  #content_shortest <- data.frame(species_name = species_name,
  #                               species_location = species_location,
  #                               shortest_sea_distance = shortest_sea_distance,
  #                               shortest_fly_distance = shortest_fly_distance)
  
  #shortest_distances <- append(shortest_distances, content_shortest)
  #return(shortest_distances)

}

error_messages <- c()

# start a list to save shortest distances in
#colnames <- c("species_name", "species_location", "shortest_sea_distance", "shortest_fly_distance")
#shortest_distances <- data.frame(matrix(ncol = length(colnames), nrow = 0))
#colnames(shortest_distances) <- colnames

result <- Map(Calculation_seadistance, long$Specieslist, long$name)
#shortest_distances <- result

# write error file
if (length(error_messages) > 0) {
  error_df <- data.frame(error_message = unlist(error_messages))
  write.csv(error_df, "test_outputs/species_error.csv", row.names = FALSE)
}

# write shortest distances file
write.table(shortest_distances, file = "test_outputs/shortest_distances.csv")





