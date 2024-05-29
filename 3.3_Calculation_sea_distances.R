
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
library("maps")
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
# RESHAPE DF WITH LOCATIONS AND SPECIES FROM WIDE TO LONG FORMAT
############################################################################################

long <- pivot_longer(df, !Specieslist)
# from tidyr package, reshape df from wide to long format
long <- long[long$value > 0, ]

############################################################################################
# REVISION: CALCULATE DISTANCES
############################################################################################

# Iterate over species_name and location_name
Calculation_seadistance <- function(species_name, species_location){
  
  #####################################################################################
  # FIRST CHECK IF FILE WITH OCCURRENCE DATA IN OCCURRENCEDATA DIRECTORY EXISTS
  # IF NOT: MAKE ONE AND GET DATA FROM GBIF
  # limit in occ_data() function is changeable for personal preference
  ####################################################################################
  
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  
  # if file exists:        put content into variable res
  # if file doesn't exist: make one, get data, and put into variable res
  occurrence_coord <- paste0("OccurrenceData_test/", species_name, ".csv")
  if (file.exists(occurrence_coord) == TRUE) {
    res <- read.csv(occurrence_coord, header = TRUE)
  } else {
    res <- occ_data(scientificName = species_name, hasCoordinate = TRUE, limit = 10000)
    
    # here, we can split the occ_data function into subsets to get more than 100000 records
    # retrieve data per category
    
    res <- res$data[, c('decimalLongitude', 'decimalLatitude')]
    print(res)
    #rename the column names
    colnames(res) <- c('Longitude', 'Latitude')
    # Remove occurrences where longitude or latitude is NA
    res <- res[!is.na(res$Latitude) & !is.na(res$Longitude),]
    # check if there's no information for a species
    if (nrow(res) == 0) {
      error_message <- paste0("No information found on GBIF for ", species_name, " in ", species_location)
      
      # check if directory with error messages exists, if it doesn't: make one
      if (!dir.exists("test_outputs/errors")){
        dir.create("test_outputs/errors", recursive = TRUE)
        error_file_name <- paste0("test_outputs/errors/error_", species_name, "_", species_location, ".csv")
        # error files written to test_outputs/errors/
        write.csv(error_message, file = error_file_name)
        return(FALSE)
        
        # write error file to the directory
      } else {
        error_file_name <- paste0("test_outputs/errors/error_", species_name, "_", species_location, ".csv")
        # error files written to test_outputs/errors/
        write.csv(error_message, file = error_file_name)
        return(FALSE)
        }
    }
    
    write.csv(res, occurrence_coord)
  
  }
  OccurrenceData <- res

  # getting coordinates from Coordinates dataframe to get longitude and latitude for ARMS location
  ################################################################################################
  
  # use grep to get the row from the Coordinates df where location_name is present
  location_row_index <- grep(species_location, Coordinates$Observatory.ID)
  # save longitude and latitude for the row that you selected with grep
  longitude <- Coordinates[location_row_index, "Longitude"]
  latitude <- Coordinates[location_row_index, "Latitude"]
  # make a dataframe out of the longitude and latitude called samplelocation
  samplelocation <- data.frame(Longitude = longitude, Latitude = latitude)
  
  # Remove duplicate longitudes & latitudes  ==> from original function
  #######################################################################
  unique_file <- OccurrenceData[!duplicated(OccurrenceData[, c("Longitude", "Latitude")]), ]
  
  # check if unique_file has coordinates
  if(nrow(unique_file) < 1) {
    stop("Unique_file has no coordinates (line 139)")
  }
  
  # make matrices out of longitudes and latitudes of samplelocation and unique_file
  # These are used in both flying distances and sea distances calculations
  sampleloc_matrix <- matrix(as.numeric(c(samplelocation$Longitude, samplelocation$Latitude)), ncol = 2)
  uniquefile_matrix <- matrix(as.numeric(c(unique_file$Longitude, unique_file$Latitude)), ncol = 2)
  
  ##########################################################################
  # FLYING DISTANCE CALCULATION
  ##########################################################################
  
  flying_distances <- as.numeric(distm(samplelocation[,c("Longitude", "Latitude")], 
                                unique_file[,c("Longitude", "Latitude")], 
                                fun = distVincentyEllipsoid))
  # 'fun =' this method is used for distances on earth (ellipsoid)
  
  # check if directory and/or csv files already exists with flying distances (files for each species/location)
  flydistance_file <- paste0("test_outputs/fly_distances/", species_name, "_fly_distancesTo_", species_location)
  
  if(!dir.exists("test_outputs/fly_distances")){
    dir.create("test_outputs/fly_distances")
    
    if(!file.exists(flydistance_file)){
      write.table(flying_distances,file = flydistance_file)
    } else {
      warning(paste0(species_name, " and ", species_location, "flying distances file already written."))
    }
    
  } else {
    if(!file.exists(flydistance_file)){
      write.table(flying_distances,file = flydistance_file)
    } else {
      warning(paste0(species_name, " and ", species_location, "flying distances file already written."))
    }
  }
  
  ##########################################################################
  # SEA DISTANCE CALCULATION
  ##########################################################################
  
  # sea_dist calculation and error catching
  sea_distances <- NA
  result <- tryCatch({
    
    ### SEA_DIST ###
    
  ########################################################################## ORIGINAL FUNCTION PART 
    #(transformed to normal running script instead of function)
    ###################
    ### SHORTESTPATH
    ###################
    
    x <- tr
    origin <- sampleloc_matrix
    goal <- uniquefile_matrix
    output <- "SpatialLines"
    
    originCells <- raster::cellFromXY(x, origin)
    goalCells <- raster::cellFromXY(x, goal)
    indexOrigin <- originCells
    indexGoal <- goalCells
    y <- transitionMatrix(x)
    if(isSymmetric(y)) {
      mode <- "undirected"
    }else{
      mode <- "directed"
    }
    
    adjacencyGraph <- igraph::graph_from_adjacency_matrix(y, mode=mode, weighted=TRUE)
    E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
    
    shortestPaths <- shortest_paths(adjacencyGraph,
                                    indexOrigin, indexGoal)$vpath
    if(output=="SpatialLines")
    {
      linesList <- vector(mode="list", length=length(shortestPaths))
      
      for(i in 1:length(shortestPaths))
      {
        sPVector <- shortestPaths[[i]]
        coords <- raster::xyFromCell(x, sPVector)
        linesList[[i]] <- Line(coords)
      }
      
      # Suggested by Sergei Petrov
      LinesObject <- mapply(Lines,
                            slinelist = linesList,
                            ID = as.character(1:length(shortestPaths)),
                            SIMPLIFY = FALSE)
      
      result <- sp::SpatialLines(LinesObject, proj4string = sp::CRS(projection(x)))
    }
    
    ################
    ### LENGTHLINE
    ################
    
    library(sp)
    library(raster)
    
    line <- result
    
    if (inherits(line, 'SpatialPolygons')) {
      requireNamespace('raster')
      line <- raster::geom(methods::as(line, 'SpatialLines'))
    } else if (inherits(line, 'SpatialLines')) {
      requireNamespace('raster')
      line <- raster::geom(line)
    } else {
      line <- cbind(object=1, part=1, cump=1, line[, 1:2])
      colnames(line)[4:5] <- c('x', 'y')
    }
    
    ids <- unique(line[,1])
    len <- rep(0, length(ids))
    
    for (i in 1:length(ids)) {
      d <- line[line[,1] == ids[i], , drop = FALSE]
      
      parts <- unique(d[,2])
      
      for (p in parts) {
        dd <- d[d[,2] == p, ,drop=FALSE]
        
        if(nrow(dd) < 2) {
          warning("Not enough points to calculate distance for part ", p, " of object ", ids[i])
          next
        }
        
        for (j in 1:(nrow(dd)-1)) {
          len[i] <- len[i] + distGeo(dd[j, c('x', 'y'), drop=FALSE], dd[j+1, c('x', 'y'), drop=FALSE])
        }
      }
    }
    sea_distances <- len
  ############################################################################ END FUNCTION PART
    
    }, error = function(e) {
      error_message <- paste0("An error occurred during calculation of sea_dist for ", species_name, " in ", species_location)
      cat(error_message, "\n")  
      
      # Collect the error message
      error_messages <<- append(error_messages, list(error_message))
      
      # Return FALSE to indicate that this iteration should be skipped
      return(FALSE)
      
  }) # trycatch() closed
  
  # check if directory and/or csv files already exists with sea distances(files for each species/location)
  distance_file <- paste0("test_outputs/sea_distances/", species_name, "_distancesTo_", species_location)
  
  if(!dir.exists("test_outputs/sea_distances")){
    dir.create("test_outputs/sea_distances")
    
    if(!file.exists(distance_file)){
      write.table(sea_distances,file = distance_file)
    } else {
      warning(paste0(species_name, " and ", species_location, "sea distances file already written."))
    }
    
  } else {
    if(!file.exists(distance_file)){
      write.table(sea_distances,file = distance_file)
    } else {
      warning(paste0(species_name, " and ", species_location, "sea distances file already written."))
    }
  }
  
  return(sea_distances)
}

error_messages <- c()

result <- Map(Calculation_seadistance, long$Specieslist, long$name)

# write error file
if (length(error_messages) > 0) {
  error_df <- data.frame(error_message = unlist(error_messages))
  write.csv(error_df, "test_outputs/species_error.csv", row.names = FALSE)
}








