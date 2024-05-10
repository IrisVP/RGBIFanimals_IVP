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
df <- read.csv("Output/Species_Location_half.csv")
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
  ### FIND CLOSEST REGISTERED PLACE
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
                   hasCoordinate = TRUE, limit = 1000) #changed limit, was 200000
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
# MAKE SAMPLELOCATION FILE WITH LATITUDE AND LONGITUDE
############################################################################################

# make a file called samplelocation with only longitude and latitude for each species
samplelocation <- data.frame(Longitude = numeric(0), Latitude = numeric(0))
for (r in 1:nrow(Coordinates)) {
  # Find the coordinates for the current Observatory ID => extract Longitude and Latitude
  current_coordinates <- Coordinates[r,2-3]
  # print(current_coordinates)   = if necessary
  # Append the coordinates to the 'samplelocation' dataframe
  samplelocation <- rbind(samplelocation, current_coordinates)
}  

############################################################################################
# FIND SHORTEST ROUTE IN SEA
############################################################################################

for (n in 1:nrow(df)){
  filename <- paste0("OccurrenceData/", df[n, 1], ".csv")
  if(!file.exists(filename)){
    write.table(paste(c(names(row), "inrange", "pointscalculated","distance"),collapse = ","), 
                file = filename, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
  }
            ### CHECK IN FILE ###
  ### checks if combination of row 1 and 2 in df is already present, if so, there is a warning
  contents <- readChar(filename, file.info(filename)$size)
  length <- length(grep(paste(df[n, 1:2], collapse = ","), contents))>0
  if(length == TRUE){   #### check if csv Occurrence data exists
    warning(paste(c(as.character(df[n, 1:2]), "has already been written"), collapse=" "))
  }
  
  # read file
  file <- read.csv(filename, header = TRUE)
  ### cleaning latitudes and longitudes from file
  file_clean <- filter(file, (file$Longitude >= -180 & file$Longitude <= 180) & (file$Latitude >= -90 & file$Latitude <= 90))
  file_clean <- filter(file, is.numeric(Longitude) & is.numeric(Latitude))
  file_clean <- file %>%
    filter(!grepl("[A-Za-z]", Longitude) & !grepl("[A-Za-z]", Latitude) & !grepl("[A-Za-z]", X))
  #file_clean$Longitude <- as.numeric(file_clean$Longitude)
  #file_clean$Latitude <- as.numeric(file_clean$Latitude)
  str(file_clean)
  print(unique(file_clean$Longitude))
  # Remove duplicates
  unique_file <- file_clean[!duplicated(file_clean[, c("Longitude", "Latitude")]), ]
  
  # Remove samples taken further away than the closest point
            ### FILTER ON DISTANCE ###
  # step1: calculate all distances to samplelocation
  distances <- as.numeric(distm(samplelocation[,c("Longitude", "Latitude")], 
                                unique_file[,c("Longitude", "Latitude")], 
                                fun = distVincentyEllipsoid))  
  ### 'fun =' which method is used = this method is for distances on earth (ellipsoid)
  # step2: Calculate the length through sea for the closest point
  ### make matrices out of longitudes and latitudes of samplelocation and unique_file
  sampleloc_matrix <- matrix(as.numeric(c(samplelocation$Longitude, samplelocation$Latitude)), ncol = 2)
  uniquefile_matrix <- matrix(as.numeric(c(unique_file$Longitude, unique_file$Latitude)), ncol = 2)
  
  sea_dist <-  geosphere::lengthLine(gdistance::shortestPath(tr, sampleloc_matrix, #function sp_format!
                                                             uniquefile_matrix, 
                                                             output = "SpatialLines"))
  ### calculated sea_dist => shortest path between sampleloc and occurrence 
  # step3: filter out the datapoints further away than the sea_dist
  filtered <- unique_file[distances <= sea_dist,]  ### retain distances that are equal to or smaller than sea_dist calculation
  ### why do I get NA values here???: distances has 740 values and sea_dist only 74!
  ### removing NA values
  cleaned_filtered <- na.omit(filtered)
  
  if(nrow(cleaned_filtered) > 0){
    OccurrenceData <- cleaned_filtered  ### if filtered rows are more than 0, then make OccurrenceData var
  } else {
    unique_file[distances <= sea_dist + 16000,] ### if there are 0 rows, save distances in unique_file
    ### that are less or equal to sea_dist + 16000
  }
  
  # Filter if there are more than 10 unique locations
  df$inrange <- nrow(OccurrenceData)
                ### FILTER N CLOSEST COORDINATE CEILING ###
  if(nrow(OccurrenceData) > 10){
    
    if(OccurrenceData$Longitude > samplelocation$Longitude){
      if (OccurrenceData$Latitude > samplelocation$Latitude){
        #print(paste("Occ_data_long - Sample_long", OccurrenceData$Longitude, "-", samplelocation$Longitude))
        #print(paste("Occ_data_lat - sample_lat", OccurrenceData$Latitude, "-", samplelocation$Latitude))
        OccurrenceData$distance <- pmax(abs(OccurrenceData$Longitude - samplelocation$Longitude), 
                                        abs(OccurrenceData$Latitude - samplelocation$Latitude))
        # pmax is used to get element-wise maximum, abs() to get absolute differences
        dist = ceiling(sort(OccurrenceData$distance)[n]) #sorted ascending and select highest distance and round up
        OccurrenceData <- OccurrenceData[OccurrenceData$distance < dist,] # keep rows where distance is less than 'dist'
        ############# DON'T GET THESE LINES #################
      }
    }
  }
  
  # Save the number of points for which the distance will be calculated
  df$pointscalculated <- nrow(OccurrenceData)
  # find the shortest route to every point through the sea
  paths <- sapply(1:nrow(OccurrenceData), function(i) {
    sampleloc_matrix <- matrix(as.numeric(c(samplelocation$Longitude, samplelocation$Latitude)), ncol = 2)
    occurrencedata_matrix <- matrix(as.numeric(c(OccurrenceData[i,]$Longitude, OccurrenceData[i,]$Latitude)), ncol = 2)
    
    path <- shortestPath(tr, sampleloc_matrix,
                         occurrencedata_matrix, 
                         output = "SpatialLines")
    return(path)
  })
  # Find the closest location the point of sampling
  df$distance <- min(as.numeric(sapply(paths, function(x) geosphere::lengthLine(x))), na.rm=T)
  write.table(paste(df,collapse = ","), file = filename, append = TRUE, quote = FALSE, 
              col.names = FALSE, row.names = FALSE)
  
}


#### Clean R environment ####
rm(list = ls())

#### Additional, optional, cleaning step ####
results <- read.csv("Output/DistanceOverSea.csv")
results <- results[!is.na(results$distance),]
write.csv(results, "Output/DistanceOverSea.csv", quote = F, row.names = F)

