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
df <- read.csv("Output/RealFirst10_Species_Location.csv")
# Read coordinates file
Coordinates <- read.csv("Inputs/Coordinates.csv")
# Read the species_location file
Species_Location <- read.csv("Output/Species_Location.csv")

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
  print(species)
  file <- "Output/ShortestPath.csv"
  if(!file.exists("Output/ShortestPath.csv")){
    write.table(paste(c("Speciesname", "Total observations", "Unique locations/observations", 
                        as.character(Coordinates$Observatory.ID), "Errors"),collapse = ","),
                file = file, append = TRUE, quote = FALSE, 
                col.names = FALSE, row.names = FALSE)
  }
  ### CHECK IN FILE
  # Check the occurrence data, If there is an error there it catches it and writes an error in the output file
  contents <- readChar(file, file.info(file)$size)
  
  if (any(grepl(paste(species, collapse = ","), contents))) {
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
      #print(paste("res$data is NULL for species:", species))
      warning("No records found in GBIF for this species")
      # Error handling: write error message to CSV
      write.table(paste(c(species, rep("", 15), "ERROR while getting the occurrence data"),
                        collapse = ","),
                  file = file, append = TRUE, quote = FALSE, 
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
        
    }
  }
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
        geom_hex(aes(x= Longitude, y = Latitude), data = res_new) +
        # coord_cartesian(xlim = c(-30, 60), ylim = c(30,70)) +
        # geom_label_repel(aes(x= Longitude, y = Latitude, label = Observatory.ID), data = Coordinates) +
        geom_point(aes(x= Longitude, y = Latitude), data = Coordinates, col = "red") + 
        ggtitle(species)
      ggsave(plotname)
    }
        
    # next part calculates line distances between samplelocations and closest occurrence data
    row_index <- grep(species, Species_Location$Specieslist) # get the row index for the species in Species_location
    selected_row <- Species_Location[row_index,] # select the row with that row index
    print(selected_row)
    distances <- c() # start a distances list
    for (col in names(selected_row[-1])){  # iterate over columns in this row
      if (selected_row[[col]]  == 0) {  # if value == 0
        distances <- c(distances, "0")  # save zero in the list
      } else {
        Location_coord <- grep(col, Coordinates$Observatory.ID) # get the location from coordinates
        Location <- Coordinates[Location_coord, 2-3] # save Longitude and Latitude in Location variable
        
        # compute distances using distVincentyEllipsoid
        dist_matrix <- distm(res_new[, c("Longitude", "Latitude")], Location, 
                              fun = distVincentyEllipsoid)
        # Find the closest point to every sampling location
        shortest <- round(apply(dist_matrix, 2, min), 0)
        distances <- c(distances, shortest)
      }
    }
    print(distances)
    # Writing results to the ShortestPath.csv file
    write.table(paste(c(species, obs, uobs, distances, ""),collapse = ","), 
                file = file, append = TRUE, quote = FALSE, 
                col.names = FALSE, row.names = FALSE)
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
error_messages <- list()
  # Iterate over the locations and species names in long dataframe
process_species_location <- function(species_name, location_name) {
  # Subset "long" dataframe for the current species and location
  species_data <- long %>% filter(Specieslist == species_name & name == location_name)
  print(species_data)
  # use grep to get the row from the Coordinates df where location_name is present
  location_row_index <- grep(location_name, Coordinates$Observatory.ID)
  # save longitude and latitude for the row that you selected with grep
  longitude <- Coordinates[location_row_index, "Longitude"]
  latitude <- Coordinates[location_row_index, "Latitude"]
  # make a dataframe out of the longitude and latitude and print as a check
  samplelocation <- data.frame(Longitude = longitude, Latitude = latitude)
  #print(paste0("Location => samplelocation$Longitude= ", samplelocation$Longitude))
  #print(paste0("Location => samplelocation$Latitude= ", samplelocation$Latitude))
  
  #####################################################################################
  # CHECK IF FILE WITH OCCURRENCE DATA IN OCCURRENCEDATA DIRECTORY EXISTS
  # IF NOT: MAKE ONE AND GET DATA FROM GBIF
  # limit in occ_data() function is changeable for personal preference
  #####################################################################################
  # if file exists: put content into variable res
  # if file doesn't exist: make one, get data, and put into variable res
  filename_occ <- paste0("OccurrenceData/", species_name, ".csv")
  if (file.exists(filename_occ) == TRUE) {
    #print(filename_occ)  # = check print
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
    write.table(paste(c(colnames(long), "inrange", "pointscalculated","distance"),collapse = ","), 
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
  #print(paste0("distances have been calculated for ", species_name, " and ", location_name, " in line 220"))
  # 'fun =' which method is used = this method is for distances on earth (ellipsoid)
  
  # SOME CHECKS
  if (any(is.na(samplelocation$Longitude)) || any(is.na(samplelocation$Latitude))) {
    stop("NA values found in samplelocation coordinates")
  }
  if (any(is.na(unique_file$Longitude)) || any(is.na(unique_file$Latitude))) {
    stop("NA values found in unique_file coordinates")
  }
  
  ### step2 ### Calculate the length through sea for the closest point
  ### make matrices out of longitudes and latitudes of samplelocation and unique_file
  sampleloc_matrix <- matrix(as.numeric(c(samplelocation$Longitude, samplelocation$Latitude)), ncol = 2)
  uniquefile_matrix <- matrix(as.numeric(c(unique_file$Longitude, unique_file$Latitude)), ncol = 2)
  
  # CHECK
  if (nrow(sampleloc_matrix) != 1 || nrow(uniquefile_matrix) == 0) {
    stop("Incorrect dimensions of the matrices")
  }
  # sea_dist calculation and error catching using function trycatch() up until line 391
  sea_dist <- NA
  result <- tryCatch({
    
    ###################################################################################### SEA_DIST
    #sea_dist <- geosphere::lengthLine(gdistance::shortestPath(tr, sampleloc_matrix, uniquefile_matrix, output = "SpatialLines"))
    
    #### Line above in comments has been replaced by ShortestPath and lengthLine functions original code
    #### Below you find ShortestPath function first and after that lengthLine
    #### This is an important adjustment, because when using the original line, not all organisms
    #### are calculated. There were too many errors. Some debug prints have been put in comments here
    
    ###################
    ### SHORTESTPATH
    ###################
    
    #function(x, origin, goal, output)
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
    #print(shortestPaths)
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
    #print(result)
    
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
    
    # Debug: Check structure of line after conversion
    #print(head(line))
    
    ids <- unique(line[,1])
    
    len <- rep(0, length(ids))
    
    # Debug: Check unique IDs and their count
    #print(ids)
    #print(length(ids))
    
    for (i in 1:length(ids)) {
      d <- line[line[,1] == ids[i], , drop = FALSE]
      
      # debug: check dimensions of d
      #print(paste0("dim(d)", dim(d)))
      
      parts <- unique(d[,2])
      
      # debug: check unique parts
      #print(parts)
      
      for (p in parts) {
        dd <- d[d[,2] == p, ,drop=FALSE]
        
        # debug: check dimensions of dd
        #print(paste0("dim(dd))", dim(dd)))
        
        if(nrow(dd) < 2) {
          warning("Not enough points to calculate distance for part ", p, " of object ", ids[i])
          next
        }
        
        for (j in 1:(nrow(dd)-1)) {
          len[i] <- len[i] + distGeo(dd[j, c('x', 'y'), drop=FALSE], dd[j+1, c('x', 'y'), drop=FALSE])
        }
      }
    }
    sea_dist <- len
    #print(paste0("sea_dist: ", sea_dist))
    #### end of function lengthLine
    
    # Continue processing if no error
    #print(paste0("sea distances have been calculated for ", species_name, " and ", location_name, " in line 230"))
    
    ### Still in trycatch() function: below is the code for the error catching part of trycatch

  }, error = function(e) {
    error_message <- paste0("An error occurred during calculation of sea_dist on line 257 for ", species_name, " in ", location_name)
    cat(error_message, "\n")  
    
    # Collect the error message
    error_messages <<- append(error_messages, list(error_message))
    
    # Return FALSE to indicate that this iteration should be skipped
    return(FALSE)
    
  }) ### trycatch() closed
  
  # Check if sea_dist was calculated successfully
  if (length(sea_dist) == 1 && is.na(sea_dist)){
    return(FALSE) # skip current iteration
  }
  
  ### calculated sea_dist => shortest path between sampleloc and occurrence 
  
  ### step3 ### filter out the datapoints further away than the sea_dist
  filtered <- unique_file[distances <= sea_dist,]
  #print("distances further than sea_dist have been filtered out")
  #print(filtered)
  ### retain distances that are equal to or smaller than sea_dist calculation
  ### NA values coming from bug in shortestPath function
  ### removing NA values
  cleaned_filtered <- na.omit(filtered)

  if(nrow(cleaned_filtered) > 0){
    OccurrenceData_new <- cleaned_filtered  ### if filtered rows are more than 0, then make OccurrenceData var
  } else {
    unique_file[distances <= sea_dist + 16000,] ### if there are 0 rows, save distances in unique_file
    ### that are less or equal to sea_dist + 16000
  }
  #print("cleaned_filtered file has been made")
  
  # save the number of rows in 'inrange'
  species_data$inrange <- nrow(OccurrenceData_new)
  #print(paste0("nrow(OccurrenceData_new) == ", species_data$inrange))
  
  #########################################################
  ### FILTER N CLOSEST COORDINATE CEILING ###
  # Filter if there are more than 10 unique locations
  #########################################################
  if(nrow(OccurrenceData_new) > 10){
    print("start calculating distances...")
    OccurrenceData_new$distance <- pmax(abs(OccurrenceData_new$Longitude - samplelocation$Longitude), 
                                        abs(OccurrenceData_new$Latitude - samplelocation$Latitude))
    #print("distances have been calculated, calculating dist and OccurrenceData_new now...")
    # pmax is used to get element-wise maximum, abs() to get absolute differences
    # Calculate the index of the nth sorted distance
    sorted_distances <- sort(OccurrenceData_new$distance)
    n <- length(sorted_distances)
    dist <- ceiling(sorted_distances[n])
    #print(paste0("dist ====", dist))
    OccurrenceData_new <- OccurrenceData_new[OccurrenceData_new$distance < dist,] 
    # keep rows where distance is less than 'dist'
  }
  ##########################################################
  
  # Save the number of points for which the distance will be calculated
  species_data$pointscalculated <- nrow(OccurrenceData_new)
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
  # Find the closest location to the point of sampling
  species_data$distance <- min(as.numeric(sapply(paths, function(x) geosphere::lengthLine(x))), na.rm=T)
  #return(long)
  write.table(paste(species_data,collapse = ","), file = filename, append = TRUE, quote = FALSE, 
              col.names = FALSE, row.names = FALSE)   # file = DistanceOverSea.csv

  return(TRUE)
}

# Apply the function to each pair of species_name and location_name
result <- Map(process_species_location, long$Specieslist, long$name)

# write error file
if (length(error_messages) > 0) {
  error_df <- data.frame(error_message = unlist(error_messages))
  write.csv(error_df, "Output/removed_species_error_new.csv", row.names = FALSE)
}
##########################################
#### Clean R environment ####
##########################################
rm(list = ls())

#### Additional, optional, cleaning step ####

results <- read.csv("Output/DistanceOverSea_realData.csv")
results <- results[!is.na(results$distance),]
write.csv(results, "Output/DistanceOverSea_realData.csv", quote = F, row.names = F)
