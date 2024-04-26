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

################################### Functions ################################################

############### part of the main script #################

# Read a species list
df <- read.csv("Output/Species_Location.csv")
# Read coordinates file
Coordinates <- read.csv("Inputs/Coordinates.csv")
# Find the shortest path for each sampling location to the species observation (as the crow flies)
sapply(df$Specieslist, function(species){
  #print(species)
  find_closest_registered_place(species, Coordinates, tr, "Output/ShortestPath.csv")
})
### coordinates are from csv file loaded and species come from df from Species_Location.csv
### tr comes from mapping world => create_rastered_world

############## part of the function script #################

check_in_file <- function(text, file){
  contents <- readChar(file, file.info(file)$size)
  return(length(grep(paste(text, collapse = ","), contents))>0)
}
######################################################
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
######################################################
get_occurrence_data <- function(species){
  res = occ_data(scientificName = species, hasCoordinate = TRUE, limit = 200000)
  res <- res$data[, c('decimalLongitude', 'decimalLatitude')]
  #rename the column names
  colnames(res) <- c('Longitude', 'Latitude')
  # Remove occurrences where longitude or latitude is NA
  res <- res[!is.na(res$Latitude) & !is.na(res$Longitude),]
  return(res)
}


#tried to only run the inside of get_occurrence_data function
sapply(df$Specieslist, function(species){
  res = occ_data(scientificName = species, hasCoordinate = TRUE, limit = 200000)
  res <- res$data[, c('decimalLongitude', 'decimalLatitude')]
  #rename the column names
  colnames(res) <- c('Longitude', 'Latitude')
  # Remove occurrences where longitude or latitude is NA
  res <- res[!is.na(res$Latitude) & !is.na(res$Longitude),]
})

######################################################

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

######################################################
find_closest_registered_place <- function(species, Coordinates, tr, outputfile, plot=TRUE){
  # Add headers to the output file if it is not made already
  if(!file.exists(outputfile)){
    write.clean.csv(c("Speciesname", "Total observations", "Unique locations", 
                      as.character(Coordinates$Observatory.ID), "Errors"), outputfile) #outputfile is shortestPath.csv
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
                             outputfile)   ### this is shortestPath.csv
  })
  obs <- nrow(res)
  # Plot the distribution =>  files made in directory Plots
  plotname <- paste0("plots/", species, ".jpeg")  ### works
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
  write.clean.csv(c(species, obs, uobs, shortest, ""), outputfile) ### shortestPath.csv
}

