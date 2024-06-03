####################################################
### LOAD LIBRARIES

library("tidyr")
library("ggplot2")
library("dplyr")
library("poliscidata")
####################################################

####################################################
### PREPARE LOCATION FILE
# Set working directory to directory where the R-script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# first load in species_location file
species_location <- read.csv("Output/next20_Species_Location.csv")
long <- pivot_longer(species_location, !Specieslist)
# from tidyr package, reshape df from wide to long format
long <- long[long$value > 0, ]




####################################################

# iterate over species_location file again using function and map()
# and read csv file per species

########################## FUNCTION

species_name <- "Aleochara obscurella"
species_location <- "Getxo"

Distribution_seadistance <- function(species_name, species_location){
  # read csv file per species
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  distance_file <- paste0("test_outputs/sea_distances/", species_name, "_distancesTo_", species_location)
  
  distances <- read.table(distance_file, header = TRUE)
  distances$sea_distances <- distances$sea_distances/1000
  #print(summary(distances))
  distances <- subset(distances, sea_distances < 40000)
  #print(summary(distances))
  
  ###############################################################################
  # make histograms of distances per species, with filtering on distance limit 40000
  ###############################################################################
  
  dist_plot <- ggplot(distances, aes(x = sea_distances)) +
    geom_histogram(binwidth = 200, fill = "blue", color = "black", boundary = 0) +
    labs(title = "Histogram of Distances", x = "Distance in km", y = "Frequency of species") +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    ggtitle(paste0("Distribution of ", species_name, " from ", species_location)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20)) +         # Set font size for axis titles
    scale_x_continuous(breaks = seq(0, 40500, by = 1500), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 40500)) # Use coord_cartesian for setting limits
    
  #print(dist_plot)
  # create a directory to save the plots in
  
  if(!dir.exists("test_outputs/sea_distribution_plots")) {
    dir.create("test_outputs/sea_distribution_plots")
    ggsave(filename = paste0("test_outputs/sea_distribution_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = dist_plot, width = 50, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("test_outputs/sea_distribution_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = dist_plot, width = 50, height = 30, units = "cm", dpi = 900)
  }
  
  ##################################################################
  # make histograms for species but with a limit set to 10000km away
  ##################################################################
  lim <- 10000
  distances_subset <- subset(distances, sea_distances < lim)  # take distances less than 10000 km away
  
  # check if there are datapoints in this subset
  if(nrow(distances_subset) == 0) {
    error_message <- paste0("There are no distances shorter than ", lim, " for ", species_name, " from ", species_location)
    error_10000 <- append(error_10000, error_message)
    return(FALSE)
  }
  
  # make histograms of distances per species, filtering on distance limit
  dist_plot_10000 <- ggplot(distances_subset, aes(x = sea_distances)) +
    geom_histogram(binwidth = 200, fill = "blue", color = "black") +
    labs(x = "Distance in km", y = "Frequency of species") +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    ggtitle(paste0("Distribution of ", species_name, " from ", species_location)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20)) +        # Set font size for axis titles
    scale_x_continuous(breaks = seq(0, lim, by = 1500), limits = c(0, lim))

  #print(dist_plot_10000)
  # create a directory to save the plots in
  
  if(!dir.exists(paste0("test_outputs/sea_distribution_limit", lim, "_plots"))) {
    dir.create(paste0("test_outputs/sea_distribution_limit", lim, "_plots"))
    ggsave(filename = paste0("test_outputs/sea_distribution_limit", lim, "_plots/plot_", species_name, "_from_", species_location, lim, ".png"), 
            plot = dist_plot_10000, width = 50, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("test_outputs/sea_distribution_limit", lim, "_plots/plot_", species_name, "_from_", species_location, lim, ".png"), 
          plot = dist_plot_10000, width = 50, height = 30, units = "cm", dpi = 900)
  }
}  
################################# END FUNCTION

# error lists initiations to store error messages in
error_10000 <- c()

# executing iteration for the long file!
plot <- Map(Distribution_seadistance, long$Specieslist, long$name)

# write error files
print(error_1000)

if(!exists("test_outputs/errors_lessthan", lim)) {
  file_error <- paste0("test_outputs/errors_lessthan", lim)
  write.table(error_10000, file = file_error, append = TRUE, quote = FALSE, 
              col.names = FALSE, row.names = FALSE)
}


#############################################################
# Make combined histogram of sea distances and fly distances
#############################################################
species_name <- "Acartia bifilosa"
species_location <- "Koster"

Distribution_combDistance <- function(species_name, species_location){
  sea_distance_file <- paste0("test_outputs/sea_distances/", species_name, "_distancesTo_", species_location)
  fly_distance_file <- paste0("test_outputs/fly_distances/", species_name, "_fly_distancesTo_", species_location)
  
  sea_distances <- read.table(sea_distance_file)
  fly_distances <- read.table(fly_distance_file)
  sea_distances$sea_distances <- sea_distances$sea_distances/1000
  sea_distances <- subset(sea_distances, sea_distances < 40000)
  sea_distances$type <- "sea_distance"
  fly_distances$fly_distances <- fly_distances$fly_distances/1000
  fly_distances <- subset(fly_distances, fly_distances < 40000)
  fly_distances$type <- "fly_distance"
  
  combined_distances <- rbind(sea_distances, fly_distances)
  
  p <- combined_distances %>%
    ggplot( aes(x=x, fill=type)) +
    geom_histogram(binwidth = 200, color="#e9ecef", alpha=0.6, position = 'identity') +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Distance in km", y = "Frequency") +
    ggtitle(paste0("Distribution of ", species_name, " from ", species_location)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
      plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
      axis.text = element_text(size = 16),  # Set font size for axis numbers
      axis.title = element_text(size = 20),
      legend.title = element_text(size = 18, face="bold"),
      legend.text = element_text(size = 16)) +
    scale_x_continuous(breaks = seq(0, 40500, by = 1500), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 40500)) + # Use coord_cartesian for setting limits
    # set legend title and labels
    scale_fill_discrete(
      name = "Distance type",
      breaks = c("sea_distance", "fly_distance"),
      labels = c("Sea distance", "Fly distance")) +
    labs(x = "Distance in km", y = "Frequency")
  
  
  if(!dir.exists(paste0("test_outputs/combined_distribution_plots"))) {
    dir.create(paste0("test_outputs/combined_distribution_plots"))
    ggsave(filename = paste0("test_outputs/combined_distribution_plots/plotcomb_", species_name, "_from_", species_location,".png"), 
           plot = p, width = 60, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("test_outputs/combined_distribution_plots/plotcomb_", species_name, "_from_", species_location,".png"), 
           plot = p, width = 60, height = 30, units = "cm", dpi = 900)
  }
}

# executing iteration for the long file!
plot <- Map(Distribution_combDistance, long$Specieslist, long$name)

#############################################################
# Make histograms of locations
#############################################################
species_name <- "Ancula gibbosa"
# Alcyonium digitatum
species_location <- "Vigo"

Location_histograms <- function(species_name, species_location){
  # read csv file per species
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  distance_file <- paste0("test_outputs/sea_distances/", species_name, "_distancesTo_", species_location)
  
  distance_file <- read.table(distance_file, header = TRUE)
  distance_file$x <- distance_file$x/1000
  
  distances <- subset(distance_file, x < 40000)
  countries <- distances$country
  #############################################################
  
  data(world)
  # Select relevant columns from the world dataset and rename them for clarity
  region_mapping <- world %>%
    filter(country %in% countries) %>%
    dplyr::select(Country = country, Region = regionun)
  
  # Join your data with the region mapping
  # hard-coding some country names from GBIF that are not linked to a region in data(world)
  df_with_regions <- distances %>%
    left_join(region_mapping, by = c("country" = "Country")) %>%
    mutate(Region = if_else(country == "United States of America", "North America", Region)) %>%
    mutate(Region = if_else(country == "Russian Federation", "Russian Federation", Region)) %>%
    mutate(Region = if_else(country == "Tonga", "Australia/New Zealand/Oceania", Region)) %>%
    mutate(Region = if_else(country == "Svalbard and Jan Mayen", "Svalbard", Region)) %>%
    mutate(Region = if_else(country == "Korea, Republic of", "Asia", Region)) %>%
    mutate(Region = if_else(country == "United Kingdom of Great Britain and Northern Ireland", "Europe", Region))
  
  
  ### PLOT
  
  country_plot <- ggplot(df_with_regions, aes(x = x, fill = Region)) +
    geom_histogram(binwidth = 500, boundary = 0, position = "stack") +  # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances/country for ", species_name," in ", species_location),
                        x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +  # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # Increase legend title size
          legend.text = element_text(size = 12),    # Increase legend text size
          legend.key.size = unit(1.5, "lines")) +   # Increase legend key size
    scale_x_continuous(breaks = seq(0, 40500, by = 1500), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 40500)) # Use coord_cartesian for setting limits
  
  print(country_plot)
  # create a directory to save the plots in
  
  if(!dir.exists("test_outputs/sea_distribution_country_plots")) {
    dir.create("test_outputs/sea_distribution_country_plots")
    ggsave(filename = paste0("test_outputs/sea_distribution_country_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = country_plot, width = 60, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("test_outputs/sea_distribution_country_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = country_plot, width = 60, height = 30, units = "cm", dpi = 900)
  }
}

# executing iteration for the long file!
plot <- Map(Location_histograms, long$Specieslist, long$name)



#############################################################
# Make histograms of year categories
#############################################################
#species_name <- "Alcyonium digitatum"
# Alcyonium digitatum
#species_location <- "Koster"

Year_histograms <- function(species_name, species_location){
  # read csv file per species
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  distance_file <- paste0("test_outputs/sea_distances/", species_name, "_distancesTo_", species_location)
  
  distance_file <- read.table(distance_file, header = TRUE)
  distance_file$sea_distances <- distance_file$sea_distances/1000
  
  distances <- subset(distance_file, sea_distances < 40000)
  years <- distances$year
  #############################################################
  assign_category <- function(year) {
    for (category in year_categories) {
      range <- as.numeric(unlist(strsplit(category, "-")))
      if (year >= range[1] & year <= range[2]) {
        return(category)
      }
    }
    return(NA) # If year doesn't fall into any category
  }
  year_categories <- c("1990-1995","1995-2000","2000-2005","2005-2010","2010-2015","2015-2020","2020-2025")
  # Apply function to create new column
  distances$year_category <- sapply(distances$year, assign_category)
  
  #################################
  ### PLOT
  
  year_plot <- ggplot(df_with_regions, aes(x = sea_distances, fill = distances$year_category)) +
    geom_histogram(binwidth = 500, boundary = 0, position = "stack") +  # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances/year for ", species_name," in ", species_location),
         x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") + # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # Increase legend title size
          legend.text = element_text(size = 12),    # Increase legend text size
          legend.key.size = unit(1.5, "lines")) +   # Increase legend key size
    scale_x_continuous(breaks = seq(0, 40500, by = 1500), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 40500)) # Use coord_cartesian for setting limits
  
  print(year_plot)
  # create a directory to save the plots in
  
  if(!dir.exists("test_outputs/sea_distribution_year_plots")) {
    dir.create("test_outputs/sea_distribution_year_plots")
    ggsave(filename = paste0("test_outputs/sea_distribution_year_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = year_plot, width = 60, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("test_outputs/sea_distribution_year_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = year_plot, width = 60, height = 30, units = "cm", dpi = 900)
  }
}

# executing iteration for the long file!
plot <- Map(Location_histograms, long$Specieslist, long$name)

