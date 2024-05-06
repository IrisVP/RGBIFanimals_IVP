################################### Load Packages ############################################
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

################################### Functions that work ################################################
create_rastered_world <- function(filename){
  # if the file is already made, it will load the file and return it
  if(file.exists(filename)){
    load(filename)
    return(tr)
  }
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
  save(tr, file = filename)
  return(tr)
}

find_closest_registered_place <- function(species, Coordinates, tr, outputfile, plot=TRUE){
  # Add headers to the output file if it is not made already
  if(!file.exists(outputfile)){
    write.clean.csv(c("Speciesname", "Total observations", "Unique locations", 
                      as.character(Coordinates$Observatory.ID), "Errors"), outputfile)
  }
  # Check the occurrence data, If there is an error there it catches it and writes an error in the output file
  if(check_in_file(species, outputfile)){
    warning(paste(species, "has already been written to the file"))
    return()
  }
  tryCatch(res <- check_occurrence_data(species),    ### line 56
           ### this function has another function with error
           error = function(errormessage){           ### new function in function
             write.clean.csv(c(species, rep("", 15), "ERROR while getting the occurrence data"),
                             outputfile)
           })
  obs <- nrow(res)
  # Plot the distribution
  plotname <- paste0("plots/", species, ".jpeg")
  if(plot & !file.exists(plotname)){
    plot_distribution(Coordinates, res, plotname, species)  ### line 72
  }
  # Remove duplicate coordinates
  res <- res[!duplicated(res),]
  uobs <- nrow(res)  # unique observations
  # Calculate the distances between all the points and all the locations
  distances <- distm(res[, c("Longitude", "Latitude")], Coordinates[, c("Longitude", "Latitude")], 
                     fun = distVincentyEllipsoid)
  # Find the closest point to every sampling location
  shortest <- round(apply(distances, 2, min), 0)
  # Writing results to the file
  write.clean.csv(c(species, obs, uobs, shortest, ""), outputfile) ### outputfile = shortestpath.csv
}

check_in_file <- function(text, file){
  contents <- readChar(file, file.info(file)$size)
  return(length(grep(paste(text, collapse = ","), contents))>0)
  #char <- grep(paste(text, collapse = ",") ,contents)
  #if(length(char) > 0){
  #return(TRUE)  # Retourneer TRUE als er overeenkomsten zijn gevonden
  #} else {
  #return(FALSE)  # Retourneer FALSE als er geen overeenkomsten zijn gevonden
  #}
}

check_occurrence_data <- function(species, filename = TRUE) {
  if (isTRUE(filename)) {
    filename <- paste0("OccurrenceData/", species, ".csv")
  }
  
  if (file.exists(filename)) {
    res <- read.csv(filename, header = TRUE)
  } else {
    try(res <- get_occurrence_data(species))
    # try(res <- get_occurrence_data(check_official_name(species)))
    if (nrow(res) == 0) stop("No information found for this species")
    write.csv(res, filename)
  }
  return(res)
}

get_occurrence_data <- function(species){
  res = occ_data(scientificName = species, hasCoordinate = TRUE, limit = 100000) # limit changeable
  res <- res$data[, c('decimalLongitude', 'decimalLatitude')]
  #rename the column names
  colnames(res) <- c('Longitude', 'Latitude')
  # Remove occurrences where longitude or latitude is NA
  res <- res[!is.na(res$Latitude) & !is.na(res$Longitude),]
  return(res)
}

plot_distribution <- function(Coordinates, res, plotname, title){
  ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group), data = map_data("world")) +
    geom_hex(aes(x= Longitude, y = Latitude), data = res) +
    # coord_cartesian(xlim = c(-30, 60), ylim = c(30,70)) +
    # geom_label_repel(aes(x= Longitude, y = Latitude, label = Observatory.ID), data = Coordinates) +
    geom_point(aes(x= Longitude, y = Latitude), data = Coordinates, col = "red") + 
    ggtitle(title)
  ggsave(plotname)
}

################################ function(s) to work on ##################

find_shortest_route_in_sea <- function(samplelocation, occurrence_data, tr, row, filename){ # filename = Output/DistanceOverSea.csv
  if(!file.exists(filename)){
    write.clean.csv(c(names(row), "inrange", "pointscalculated","distance"), filename) # filename = Output/DistanceOverSea.csv
  }
  ### checks if combination of row1 and 2 is already present, if so, there is a warning
  if(check_in_file(row[1:2], filename)){
    warning(paste(c(as.character(row[1:2]), "has already been written"), collapse=" "))
    return()
  }
  # Remove duplicates
  occurrence_data <- occurrence_data[!duplicated(occurrence_data),]
  
  ## originally filter_on_distance function
  ############### Remove samples taken further away than the closest point ##########
  
  # step1: calculate all distances to samplelocation
  distances <- as.numeric(distm(samplelocation[,c("Longitude", "Latitude")], 
                                occurrence_data[,c("Longitude", "Latitude")], 
                                fun = distVincentyEllipsoid))  
  ### 'fun =' which method is used = this method is for distances on earth (ellipsoid)
  # step2: Calculate the length through sea for the closest point
  ### make matrices out of longitudes and latitudes of samplelocation and unique_file
  sampleloc_matrix <- matrix(as.numeric(c(samplelocation$Longitude, samplelocation$Latitude)), ncol = 2)
  occurrence_matrix <- matrix(as.numeric(c(occurrence_data$Longitude, occurrence_data$Latitude)), ncol = 2)
  
  sea_dist <-  geosphere::lengthLine(gdistance::shortestPath(tr, sampleloc_matrix, #function sp_format!
                                                             occurrence_matrix, 
                                                             output = "SpatialLines"))
  ### calculated sea_dist => shortest path between sampleloc and occurrence 
  # step3: filter out the datapoints further away than the sea_dist
  occurrence_data <- occurrence_data[distances <= sea_dist,]  ### retain distances that are equal to or smaller than sea_dist calculation
  ### why do I get NA values here???: distances has 740 values and sea_dist only 74!
  if(nrow(occurrence_data) <= 0){
    occurrence_data[distances <= sea_dist + 16000,]
  }
  ###################################################################################
  
  # Filter if there are more than 10 unique locations
  df$inrange <- nrow(occurrence_data)
  if(nrow(occurrence_data) > 10){
    #occurrence_data <- filter_n_closest_coordinate_ceiling(10, occurrence_data, samplelocation)
    occurrence_data$distance <- pmax(abs(occurrence_data$Longitude - samplelocation$Longitude), 
                                     abs(occurrence_data$Latitude - samplelocation$Latitude))
    # pmax is used to get element-wise maximum, abs() to get absolute differences
    # Get the n-th smallest distance and round it up
    dist = ceiling(sort(occurrence_data$distance)[n]) #sorted ascending and select smallest distance and round up
    occurrence_data <- occurrence_data[occurrence_data$distance < dist,] # keep rows where distance is less than 'dist'
    
  }
  # Save the number of points for which the distance will be calculated
  df$pointscalculated <- nrow(occurrence_data)
  # find the shortest route to every point through the sea
  paths <- sapply(1:nrow(occurrence_data), function(i) {
    sampleloc_matrix <- matrix(as.numeric(c(samplelocation$Longitude, samplelocation$Latitude)), ncol = 2)
    occurrencedata_matrix <- matrix(as.numeric(c(occurrence_data[i,]$Longitude, occurrence_data[i,]$Latitude)), ncol = 2)
    
    path <- shortestPath(tr, sampleloc_matrix,
                         occurrencedata_matrix, 
                         output = "SpatialLines")
    return(path)
  })
  # Find the closest location the point of sampling
  df$distance <- min(as.numeric(sapply(paths, function(x) geosphere::lengthLine(x))), na.rm=T)
  write.clean.csv(row, filename)
}

write.clean.csv <- function(list, outputfile){
  write.table(paste(c(list),collapse = ","), file = outputfile, append = TRUE, quote = FALSE, 
              col.names = FALSE, row.names = FALSE)
}

check_official_name <- function(species){
  return(worrmsbynames(species)$valid_name)
}


