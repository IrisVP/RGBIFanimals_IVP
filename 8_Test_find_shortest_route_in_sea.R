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

################ sample location, main script part #########################

# Initialize an empty dataframe to store sample locations
samplelocation <- data.frame(Longitude = numeric(0), Latitude = numeric(0))
for (r in 1:nrow(Coordinates)) {
  # Find the coordinates for the current Observatory ID => extract Longitude and Latitude
  current_coordinates <- Coordinates[r,2-3]
  print(current_coordinates)
  # Append the coordinates to the 'samplelocation' dataframe
  samplelocation <- rbind(samplelocation, current_coordinates)
}
#### part above = SOLVED

  ### CHECK OCCURRENCE DATA ###
for (n in 1:nrow(long)){
  filename <- paste0("OccurrenceData/", long[n, 1], ".csv")
  print(paste("filename (from long) =", filename))
  
    ### GET OCCURRENCE DATA ###
  res = occ_data(scientificName = long[n, 1],
                hasCoordinate = TRUE, #Return only occurrence records with lat/long data (TRUE) or all records (FALSE, default).
                limit = 500) # number of records to return, default = 500
  # Note that there is a hard maximum of 100,000, which is calculated as the limit+start, so start=99,000 and limit=2000 won't work
  print(paste("species name (from long) =", long[n, 1]))
    
  res_new <- res$data[, c('decimalLongitude', 'decimalLatitude')]  # there is no data column made
  #rename the column names      
  colnames(res_new) <- c('Longitude', 'Latitude')
  # Remove occurrences where longitude or latitude is NA
  res_new <- res_new[!is.na(res_new$Latitude) & !is.na(res_new$Longitude),]
  if (nrow(res_new) == 0) stop("No information found for this species")
  write.csv(res_new, filename)
}
#####################################################################################
###find_shortest_route_in_sea(samplelocation, occurrence_data, tr, row, "Output/DistanceOverSea.csv")

for (n in 1:nrow(df)){
  filename <- paste0("OccurrenceData/", df[n, 1], ".csv")
  if(!file.exists(filename)){
    write.clean.csv(c(names(row), "inrange", "pointscalculated","distance"), filename)
  }
  ### checks if combination of row1 and 2 is already present, if so, there is a warning
  if(check_in_file(df[n, 1:2], filename)){   #### check if csv Occurrence data exists
    warning(paste(c(as.character(df[n, 1:2]), "has already been written"), collapse=" "))
  }
  # Remove duplicates
  file <- read.csv(filename, header = TRUE)
  unique_file <- file[!duplicated(file[, c("Longitude", "Latitude")]), ]  # unique longitudes + latitudes
  # OccurrenceData <- OccurrenceData[!duplicated(OccurrenceData),]    ==== original lines
  # Remove samples taken further away than the closest point
  OccurrenceData <- filter_on_distance(tr, samplelocation, unique_file) ### return filtered into OccurrenceData
  # Filter if there are more than 10 unique locations
  df$inrange <- nrow(OccurrenceData)
  if(nrow(OccurrenceData) > 10){
    ###OccurrenceData <- filter_n_closest_coordinate_ceiling(10, OccurrenceData, samplelocation)
    OccurrenceData$distance <- pmax(abs(OccurrenceData$Longitude - samplelocation$Longitude), 
                                     abs(OccurrenceData$Latitude - samplelocation$Latitude))
    # pmax is used to get element-wise maximum, abs() to get absolute differences
    # Get the n-th smallest distance and round it up
    dist = ceiling(sort(OccurrenceData$distance)[n]) #sorted ascending and select smallest distance and round up
    OccurrenceData <- OccurrenceData[OccurrenceData$distance < dist,] # keep rows where distance is less than 'dist'
    return(OccurrenceData)
  }
}


# Save the number of points for which the distance will be calculated
row$pointscalculated <- nrow(occurrence_data)
# find the shortest route to every point through the sea
paths <- sapply(1:nrow(occurrence_data), function(i) {
  path <- shortestPath(tr, sp_format(samplelocation), 
                       sp_format(occurrence_data[i,]), 
                       output = "SpatialLines")
  return(path)
})
# Find the closest location the point of sampling
row$distance <- min(as.numeric(sapply(paths, function(x) geosphere::lengthLine(x))), na.rm=T)
write.clean.csv(row, filename)




