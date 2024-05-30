####################################################
### LOAD LIBRARIES

library("tidyr")
library("ggplot2")
library("dplyr")
####################################################

####################################################
### PREPARE LOCATION FILE
# Set working directory to directory where the R-script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# first load in species_location file
species_location <- read.csv("Output/RealFirst10_Species_Location.csv")
long <- pivot_longer(species_location, !Specieslist)
# from tidyr package, reshape df from wide to long format
long <- long[long$value > 0, ]
####################################################

# iterate over species_location file again using function and map()
# and read csv file per species

########################## FUNCTION
Distribution_seadistance <- function(species_name, species_location){
  # read csv file per species
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  distance_file <- paste0("test_outputs/sea_distances_more/", species_name, "_distancesTo_", species_location)
  
  distances <- read.table(distance_file, header = TRUE)
  distances$x <- distances$x/1000
  #print(summary(distances))
  distances <- subset(distances, x < 40000)
  #print(summary(distances))
  
  ###############################################################################
  # make histograms of distances per species, with filtering on distance limit 40000
  ###############################################################################
  
  dist_plot <- ggplot(distances, aes(x = x)) +
    geom_histogram(binwidth = 200, fill = "blue", color = "black", boundary = 0) +
    labs(title = "Histogram of Distances", x = "Distance in km", y = "Frequency of species") +
    theme_bw() +
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
  
  if(!dir.exists("test_outputs/sea_distribution_plots_more")) {
    dir.create("test_outputs/sea_distribution_plots_more")
    ggsave(filename = paste0("test_outputs/sea_distribution_plots_more/plot_", species_name, "_from_", species_location, ".png"), 
           plot = dist_plot, width = 50, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("test_outputs/sea_distribution_plots_more/plot_", species_name, "_from_", species_location, ".png"), 
           plot = dist_plot, width = 50, height = 30, units = "cm", dpi = 900)
  }
  
  ##################################################################
  # make histograms for species but with a limit set to 10000km away
  ##################################################################
  lim <- 10000
  distances_subset <- subset(distances, x < lim)  # take distances less than 10000 km away
  
  # check if there are datapoints in this subset
  if(nrow(distances_subset) == 0) {
    error_message <- paste0("There are no distances shorter than ", lim, " for ", species_name, " from ", species_location)
    error_10000 <- append(error_10000, error_message)
    return(FALSE)
  }
  
  # make histograms of distances per species, filtering on distance limit
  dist_plot_10000 <- ggplot(distances_subset, aes(x = x)) +
    geom_histogram(binwidth = 200, fill = "blue", color = "black") +
    labs(x = "Distance in km", y = "Frequency of species") +
    theme_bw() +
    ggtitle(paste0("Distribution of ", species_name, " from ", species_location)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20)) +        # Set font size for axis titles
    scale_x_continuous(breaks = seq(0, lim, by = 1500), limits = c(0, lim))
  
  #print(dist_plot_10000)
  # create a directory to save the plots in
  
  if(!dir.exists(paste0("test_outputs/sea_distribution_limit", lim, "_plots_more"))) {
    dir.create(paste0("test_outputs/sea_distribution_limit", lim, "_plots_more"))
    ggsave(filename = paste0("test_outputs/sea_distribution_limit", lim, "_plots_more/plot_", species_name, "_from_", species_location, lim, ".png"), 
           plot = dist_plot_10000, width = 50, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("test_outputs/sea_distribution_limit", lim, "_plots_more/plot_", species_name, "_from_", species_location, lim, ".png"), 
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

Distribution_combDistance <- function(species_name, species_location){
  sea_distance_file <- paste0("test_outputs/sea_distances_more/", species_name, "_distancesTo_", species_location)
  fly_distance_file <- paste0("test_outputs/fly_distances_more/", species_name, "_fly_distancesTo_", species_location)
  
  sea_distances <- read.table(sea_distance_file)
  fly_distances <- read.table(fly_distance_file)
  sea_distances$x <- sea_distances$x/1000
  sea_distances <- subset(sea_distances, x < 40000)
  sea_distances$type <- "sea_distance"
  fly_distances$x <- fly_distances$x/1000
  fly_distances <- subset(fly_distances, x < 40000)
  fly_distances$type <- "fly_distance"
  
  combined_distances <- rbind(sea_distances, fly_distances)
  
  p <- combined_distances %>%
    ggplot( aes(x=x, fill=type)) +
    geom_histogram(binwidth = 200, color="#e9ecef", alpha=0.6, position = 'identity') +
    theme_bw() +
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
  
  
  if(!dir.exists(paste0("test_outputs/combined_distribution_plots_more"))) {
    dir.create(paste0("test_outputs/combined_distribution_plots_more"))
    ggsave(filename = paste0("test_outputs/combined_distribution_plots_more/plotcomb_", species_name, "_from_", species_location,".png"), 
           plot = p, width = 50, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("test_outputs/combined_distribution_plots_more/plotcomb_", species_name, "_from_", species_location,".png"), 
           plot = p, width = 50, height = 30, units = "cm", dpi = 900)
  }
}

# executing iteration for the long file!
plot <- Map(Distribution_combDistance, long$Specieslist, long$name)
