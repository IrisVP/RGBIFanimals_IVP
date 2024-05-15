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

############################################################################################
# LOAD DATA
############################################################################################
# Set working directory to directory where the R-script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # requires installation of package "rstudioapi"
# Read a species list
df <- read.csv("Output/Test_error_Species_Location.csv")
# Read coordinates file
Coordinates <- read.csv("Inputs/Coordinates.csv")

############################################################################################
# CREATE A RASTERED WORLD
############################################################################################

# if the file is already made, it will load the file and return it
if(file.exists("Inputs/new_tr.rdata")){
  load("Inputs/new_tr.rdata")
} else{
  # Load a map
  wrld_simpl2 <- map_data("world") # comes from ggplot2
  wrld_simpl <- df_to_SpatialPolygons(df = wrld_simpl2,
                                      keys = "region",
                                      coords = c("long", "lat"),
                                      proj = CRS("+proj=longlat +ellps=sphere"))
  # data(wrld_simpl)
  # Generate a scaffold for the raster file
  world_crs <- crs(wrld_simpl)
  worldshp <- spTransform(wrld_simpl, world_crs)
  ras <- raster(nrow=1200, ncol=1200)
  
  # Generate a raster file
  worldmask <- rasterize(worldshp, ras)
  worldras <- is.na(worldmask) # inverse water and land, so ocean becomes 1 and land 0
  worldras[worldras==0] <- 999 # set land to 999

  # Create a Transition object from the raster
  tr <- transition(worldras, function(x) 1/mean(x), 16)
  tr = geoCorrection(tr, scl=FALSE)
}
############################################################################################
# FIND SHORTEST PATH FOR EACH SAMPLING LOCATION TO THE SPECIES OBSERVATION
############################################################################################
sapply(df$Specieslist, function(species){
  ### FIND CLOSEST REGISTERED PLACE ###
  # Add headers to the output file if it is not made already
  #species <- "Halichondria bowerbanki"
  print(species)
  if(!file.exists("Output/ShortestPath.csv")){
    write.table(paste(c("Speciesname", "Total observations", "Unique locations/observations", 
                        as.character(Coordinates$Observatory.ID), "Errors"),collapse = ","),
                file = "Output/ShortestPath.csv", append = TRUE, quote = FALSE, 
                col.names = FALSE, row.names = FALSE)
  }
  ### CHECK IN FILE
  # Check the occurrence data, If there is an error there it catches it and writes an error in the output file
  file <- "Output/ShortestPath.csv"
  contents <- readChar(file, file.info(file)$size)
  length <- length(grep(paste(species, collapse = ","), contents))
  
  if(length > 0){
      warning(paste(species, "has already been written to the file"))
  }
  ### CHECK OCCURRENCE DATA ###
  filename <- paste0("OccurrenceData/", species, ".csv")
  
  if (file.exists(filename)) {
    res_new <- read.csv(filename, header = TRUE)
  } else {  
    ### GET OCCURRENCE DATA ###
    res <- occ_data(scientificName = species,
                   hasCoordinate = TRUE, limit = 10000) #changed limit, was 200000
    if (is.null(res) || !is.list(res) || !is.list(res$data)) {
      print(paste("res$data is NULL for species:", species))
      warning("No records found in GBIF for this species")
      # Error handling: write error message to CSV
      write.table(paste(c(species, rep("", 15), "ERROR while getting the occurrence data"),
                        collapse = ","),
                  file = "Output/ShortestPath.csv", append = TRUE, quote = FALSE, 
                  col.names = FALSE, row.names = FALSE)
      
    } else {
        ### MAKE COLUMN NAMES AND REMOVE NAs
        res_new <- res$data[, c('decimalLongitude', 'decimalLatitude')]
        #rename the column names
        colnames(res_new) <- c('Longitude', 'Latitude')
        # Remove occurrences where longitude or latitude is NA
        res_new <- res_new[!is.na(res_new$Latitude) & !is.na(res_new$Longitude),]
        # Only keep the longitudes between -180 & 180 and latitudes between -90 & 90
        res_new <- res_new[res_new$Longitude >= -180 & res_new$Longitude <= 180 & 
                             res_new$Latitude >= -180 & res_new$Latitude <= 180, ]
        # Save number of rows in obs
        obs <- nrow(res_new)
        # Remove duplicate coordinates
        res_new <- res_new[!duplicated(res_new[, c("Longitude", "Latitude")]),]
        # Save number of UNIQUE rows in uobs
        uobs <- nrow(res_new)
        
        # try(res <- get_occurrence_data(check_official_name(species))) ERROR
        if (nrow(res_new) == 0) stop("No information found for this species")
        write.csv(res_new, filename)
        
        # Plot the distribution
        plotname <- paste0("plots/", species, ".jpeg")
        if(!file.exists(plotname)){
          ggplot() +
            geom_polygon(aes(x = long, y = lat, group = group), data = map_data("world")) +
            geom_hex(aes(x= Longitude, y = Latitude), data = res) +
            # coord_cartesian(xlim = c(-30, 60), ylim = c(30,70)) +
            # geom_label_repel(aes(x= Longitude, y = Latitude, label = Observatory.ID), data = Coordinates) +
            geom_point(aes(x= Longitude, y = Latitude), data = Coordinates, col = "red") + 
            ggtitle(title)
          ggsave(plotname)
        }
        # Calculate the distances between all the points and all the locations
        distances <- distm(res_new[, c("Longitude", "Latitude")], Coordinates[, c("Longitude", "Latitude")], 
                           fun = distVincentyEllipsoid)
        # Find the closest point to every sampling location
        shortest <- round(apply(distances, 2, min), 0)
        # Writing results to the ShortestPath.csv file
        write.table(paste(c(species, obs, uobs, shortest, ""),collapse = ","), 
                    file = "Output/ShortestPath.csv", append = TRUE, quote = FALSE, 
                    col.names = FALSE, row.names = FALSE)
        
        }
    }
    
})

############################################################################################
# RESHAPE DF FROM WIDE TO LONG FORMAT
############################################################################################

long <- pivot_longer(df, !Specieslist)# from tidyr package, reshape df from wide to long format
# all columns in df except Specieslist should be pivoted => will be treated as value columns
# names will be stored in a new column called "name"
long <- long[long$value > 0, ]

############################################################################################
# FIND SHORTEST ROUTE IN SEA PART IN MAIN SCRIPT => UNTIL THE END OF SCRIPT
############################################################################################

#############################################
### CHECK OCCURRENCE DATA AND SAVE IT ###
#############################################

  # Iterate over the locations and species names in long dataframe
Map(function(species_name, location_name) {
  # Subset "long" dataframe for the current species and location
  species_data <- long %>% filter(Specieslist == species_name & name == location_name)
  # print location name as a check
  print(paste0("Location_name for species in df long= ", location_name))
  # use grep to get the row from the Coordinates df where location_name is present
  location_row_index <- grep(location_name, Coordinates$Observatory.ID)
  # save longitude and latitude for the row that you selected with grep
  longitude <- Coordinates[location_row_index, "Longitude"]
  latitude <- Coordinates[location_row_index, "Latitude"]
  # make a dataframe out of the longitude and latitude and print as a check
  samplelocation <- data.frame(Longitude = longitude, Latitude = latitude)
  print(paste0("Location => samplelocation$Longitude= ", samplelocation$Longitude))
  print(paste0("Location => samplelocation$Latitude= ", samplelocation$Latitude))
  
  #####################################################################################
  # CHECK IF FILE WITH OCCURRENCE DATA IN OCCURRENCEDATA DIRECTORY EXISTS
  # IF NOT: MAKE ONE AND GET DATA FROM GBIF
  # limit in occ_data() function is changeable for personal preference
  #####################################################################################
  # if file exists: put content into variable res
  # if file doesn't exist: make one, get data, and put into variable res
  filename_occ <- paste0("OccurrenceData/", species_name, ".csv")
  if (file.exists(filename_occ) == TRUE) {
    print(filename_occ)  # = check print
    res <- read.csv(filename_occ, header = TRUE)
  } else {
    res <- occ_data(scientificName = species_name, hasCoordinate = TRUE, limit = 100000) #changed limit, was 200000
    res <- res$data[, c('decimalLongitude', 'decimalLatitude')]
    #rename the column names
    colnames(res) <- c('Longitude', 'Latitude')
    # Remove occurrences where longitude or latitude is NA
    res <- res[!is.na(res$Latitude) & !is.na(res$Longitude),]
    # try(res <- get_occurrence_data(check_official_name(species)))
    if (nrow(res) == 0) stop("No information found for this species")
    write.csv(res, filename_occ)
      
  }
  OccurrenceData <- res   # working with OccurrenceData in the next steps
  
  ###########################################################################################
  ### FIND SHORTEST ROUTE IN SEA ###
  # Is the main function in this part: withholds the next few steps until the end of the loop
  ###########################################################################################
  
  filename <- "Output/DistanceOverSea.csv"
  if(!file.exists(filename)){
    write.table(paste(c(names(colnames(long)), "inrange", "pointscalculated","distance"),collapse = ","), 
                file = filename, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
  }
  ##############################################
  ### CHECK IN FILE ###
  ##############################################
  ### checks if combination of column 1 and 2 from df 'long' is already present in DistanceOverSea.csv
  ### if so, there is a warning
  contents <- readChar(filename, file.info(filename)$size)
  length <- length(grep(paste(long[, 1:2], collapse = ","), contents))>0
  if(length == TRUE){
    warning(paste(c(as.character(long[, 1:2]), "has already been written to DistanceOverSea.csv"), collapse=" "))
  }
  
  # Remove duplicate longitudes & latitudes  ==> from original function
  unique_file <- OccurrenceData[!duplicated(OccurrenceData[, c("Longitude", "Latitude")]), ]
  
  ############################################################
  ### FILTER ON DISTANCE ###
  # Remove samples taken further away than the closest point
  ############################################################
  ### step1 ### calculate all distances to samplelocation
  distances <- as.numeric(distm(samplelocation[,c("Longitude", "Latitude")], 
                                unique_file[,c("Longitude", "Latitude")], 
                                fun = distVincentyEllipsoid))
  print(paste0("distances have been calculated for", species_name, "and", location_name, "in line 220"))
  # 'fun =' which method is used = this method is for distances on earth (ellipsoid)
  
  ### step2 ### Calculate the length through sea for the closest point
  ### make matrices out of longitudes and latitudes of samplelocation and unique_file
  sampleloc_matrix <- matrix(as.numeric(c(samplelocation$Longitude, samplelocation$Latitude)), ncol = 2)
  uniquefile_matrix <- matrix(as.numeric(c(unique_file$Longitude, unique_file$Latitude)), ncol = 2)
  
  sea_dist <-  geosphere::lengthLine(gdistance::shortestPath(tr, sampleloc_matrix, #function sp_format!
                                                             uniquefile_matrix, 
                                                             output = "SpatialLines"))
  print(paste0("sea distances have been calculated for ", species_name, " and ", location_name, " in line 230"))
  ### calculated sea_dist => shortest path between sampleloc and occurrence 
  
  ### step3 ### filter out the datapoints further away than the sea_dist
  filtered <- unique_file[distances <= sea_dist,]
  print("distances further than sea_dist have been filtered out")
  ### retain distances that are equal to or smaller than sea_dist calculation
  ### why do I get NA values here?
  ### removing NA values
  cleaned_filtered <- na.omit(filtered)

  if(nrow(cleaned_filtered) > 0){
    OccurrenceData_new <- cleaned_filtered  ### if filtered rows are more than 0, then make OccurrenceData var
  } else {
    unique_file[distances <= sea_dist + 16000,] ### if there are 0 rows, save distances in unique_file
    ### that are less or equal to sea_dist + 16000
  }
  print("cleaned_filtered file has been made")
  
  # make new long dataframe and only put in the specific species with only 1 location
  new_long <- filter(long, name == location_name)
  # Filter if there are more than 10 unique locations
  new_long$inrange <- nrow(OccurrenceData_new)
  print(paste0("nrow(OccurrenceData_new) == ", new_long$inrange))
  
  #########################################################
  ### FILTER N CLOSEST COORDINATE CEILING ###
  #########################################################
  if(nrow(OccurrenceData_new) > 10){
    print("start calculating distances...")
    OccurrenceData_new$distance <- pmax(abs(OccurrenceData_new$Longitude - samplelocation$Longitude), 
                                        abs(OccurrenceData_new$Latitude - samplelocation$Latitude))
    print("distances have been calculated, calculating dist and OccurrenceData_new now...")
    # pmax is used to get element-wise maximum, abs() to get absolute differences
    # Calculate the index of the nth sorted distance
    sorted_distances <- sort(OccurrenceData_new$distance)
    n <- length(sorted_distances)  # Replace 'n' with the desired index
    dist <- ceiling(sorted_distances[n])
    print(paste0("dist ====", dist))
    OccurrenceData_new <- OccurrenceData_new[OccurrenceData_new$distance < dist,] 
    # keep rows where distance is less than 'dist'
  }
  ##########################################################
  
  # Save the number of points for which the distance will be calculated
  new_long$pointscalculated <- nrow(OccurrenceData_new)
  # find the shortest route to every point through the sea
  paths <- sapply(1:nrow(OccurrenceData_new), function(i) {
    sampleloc_matrix <- matrix(as.numeric(c(samplelocation$Longitude, samplelocation$Latitude)), ncol = 2)
    occurrencedata_matrix <- matrix(as.numeric(c(OccurrenceData_new[i,]$Longitude, OccurrenceData_new[i,]$Latitude)), ncol = 2)
    
    path <- shortestPath(tr, sampleloc_matrix,
                         occurrencedata_matrix, 
                         output = "SpatialLines")
    return(path)
  })
  print("Shortest routes through sea for every point have been calculated")
  # Find the closest location the point of sampling
  new_long$distance <- min(as.numeric(sapply(paths, function(x) geosphere::lengthLine(x))), na.rm=T)
  return(new_long)
  write.table(paste(new_long,collapse = ","), file = filename, append = TRUE, quote = FALSE, 
              col.names = FALSE, row.names = FALSE)   # file = DistanceOverSea.csv
}, new_long$Specieslist, new_long$name)

##########################################
#### Clean R environment ####
##########################################
rm(list = ls())

#### Additional, optional, cleaning step ####
results <- read.csv("Output/DistanceOverSea.csv")
results <- results[!is.na(results$distance),]
write.csv(results, "Output/DistanceOverSea.csv", quote = F, row.names = F)
