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
species_name <- "Acartia margalefi"     # organism that has few occurrences??
species_location <- "SwedishWestCoast"

########################## FUNCTION
Distribution_seadistance <- function(species_name, species_location){
  # read csv file per species
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  distance_file <- paste0("test_outputs/sea_distances/", species_name, "_distancesTo_", species_location)
  
  distances <- read.table(distance_file)
  distances$x <- distances$x/1000
  
  ###############################################################################
  # make histograms of distances per species, without filtering on distance limit
  ###############################################################################
  
  dist_plot <- ggplot(distances, aes(x = x)) +
    geom_histogram(binwidth = 200, fill = "blue", color = "black") +
    labs(title = "Histogram of Distances", x = "Distance", y = "Frequency") +
    theme_minimal() +
    ggtitle(paste0("Distribution of ", species_name, " from ", species_location)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 17),           # Set font size for axis numbers
          axis.title = element_text(size = 20),          # Set font size for axis titles
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),) +
    ylab("Frequency of species") +
    xlab("Distance in km")
    
  #print(dist_plot)
  # create a directory to save the plots in
  
  if(!dir.exists("test_outputs/distribution_plots")) {
    dir.create("test_outputs/distribution_plots")
    if (!exists(paste0("test_outputs/distribution_plots/plot_", species_name, "_from_", species_location, ".png"))) {
      ggsave(filename = paste0("test_outputs/distribution_plots/plot_", species_name, "_from_", species_location, ".png"), 
             plot = dist_plot, width = 40, height = 30, units = "cm", dpi = 900)
    } else {
      warning(paste0(species_name, " for ", species_location, " sea distance plot already created"))
    }
  } else {
    if (!exists(paste0("test_outputs/distribution_plots/plot_", species_name, "_from_", species_location, ".png"))) {
      ggsave(filename = paste0("test_outputs/distribution_plots/plot_", species_name, "_from_", species_location, ".png"), 
             plot = dist_plot, width = 40, height = 30, units = "cm", dpi = 900)
    } else {
      warning(paste0(species_name, " for ", species_location, " sea distance plot already created"))
    }
  }
  
  ##################################################################
  # make histograms for species but with a limit set to 10000km away
  ##################################################################
  lim <- 10000
  distances_subset <- as.data.frame(subset(distances, x < lim))  # take distances less than 10000 km away
  
  # check if there are datapoints in this subset
  if(nrow(distances_subset) == 0) {
    error_message <- paste0("There are no distances shorter than ", lim, " for ", species_name, " from ", species_location)
    error_10000 <- append(error_10000, error_message)
    return(FALSE)
  }
  
  # make histograms of distances per species, filtering on distance limit
  dist_plot_10000 <- ggplot(distances_subset, aes(x = x)) +
    geom_histogram(binwidth = 200, fill = "blue", color = "black") +
    labs(x = "Distance", y = "Frequency") +
    theme_minimal() +
    ggtitle(paste0("Distribution of ", species_name, " from ", species_location)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 17),           # Set font size for axis numbers
          axis.title = element_text(size = 20),          # Set font size for axis titles
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),) +
    ylab("Frequency of species") +
    xlab("Distance in km")

  #print(dist_plot_10000)
  # create a directory to save the plots in
  
  if(!dir.exists(paste0("test_outputs/distribution_limit", lim, "_plots"))) {
    dir.create(paste0("test_outputs/distribution_limit", lim, "_plots"))
    if (!exists(paste0("test_outputs/distribution_limit", lim, "_plots/plot_", species_name, "_from_", species_location, lim, ".png"))) {
      ggsave(filename = paste0("test_outputs/distribution_limit", lim, "_plots/plot_", species_name, "_from_", species_location, lim, ".png"), 
              plot = dist_plot_10000, width = 40, height = 30, units = "cm", dpi = 900)
    } else {
      warning(paste0(species_name, " for ", species_location, " sea distance plot with limit already created"))
    }
  } else {
    if (!exists(paste0("test_outputs/distribution_limit", lim, "_plots/plot_", species_name, "_from_", species_location, lim, ".png"))) {
      ggsave(filename = paste0("test_outputs/distribution_limit", lim, "_plots/plot_", species_name, "_from_", species_location, lim, ".png"), 
            plot = dist_plot_10000, width = 40, height = 30, units = "cm", dpi = 900)
    } else {
      warning(paste0(species_name, " for ", species_location, " sea distance plot with limit already created"))
    }
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
species_name <- "Acartia bifilosa"     # organism that has few occurrences??
species_location <- "Roscoff"

Distribution_combDistance <- function(species_name, species_location){
  sea_distance_file <- paste0("test_outputs/sea_distances/", species_name, "_distancesTo_", species_location)
  fly_distance_file <- paste0("test_outputs/fly_distances/", species_name, "_fly_distancesTo_", species_location)
  
  sea_distances <- read.table(sea_distance_file)
  fly_distances <- read.table(fly_distance_file)
  sea_distances$x <- sea_distances$x/1000
  sea_distances$type <- "sea_distance"
  fly_distances$x <- fly_distances$x/1000
  fly_distances$type <- "fly_distance"
  
  combined_distances <- rbind(sea_distances, fly_distances)
  
  p <- combined_distances %>%
    ggplot( aes(x=x, fill=type)) +
    geom_histogram(binwidth = 200, color="#e9ecef", alpha=0.6, position = 'identity') +
    scale_fill_manual(values=c("#69b3a2", "#404080")) +
    theme_bw() +
    labs(x = "Distance", y = "Frequency") +
    ggtitle(paste0("Distribution of ", species_name, " from ", species_location)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
      plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
      axis.text = element_text(size = 12),  # Set font size for axis numbers
      axis.title = element_text(size = 20),
      ) +
    scale_x_continuous(breaks = seq(0, max(combined_distances$x), by = 2000)) # Add extra values to the x-axis
  print(p)
  
  #print(p)
  if(!dir.exists(paste0("test_outputs/combined_distribution_plots"))) {
    dir.create(paste0("test_outputs/combined_distribution_plots"))
    if (!exists(paste0("test_outputs/combined_distribution_plots/plotcomb_", species_name, "_from_", species_location,".png"))) {
      ggsave(filename = paste0("test_outputs/combined_distribution_plots/plotcomb_", species_name, "_from_", species_location,".png"), 
             plot = p, width = 40, height = 30, units = "cm", dpi = 900)
    } else {
      warning(paste0(species_name, " for ", species_location, " combined distance plot already created"))
    }
  } else {
    if (!exists(paste0("test_outputs/combined_distribution_plots/plotcomb_", species_name, "_from_", species_location,".png"))) {
      ggsave(filename = paste0("test_outputs/combined_distribution_plots/plotcomb_", species_name, "_from_", species_location,".png"), 
             plot = p, width = 40, height = 30, units = "cm", dpi = 900)
    } else {
      warning(paste0(species_name, " for ", species_location, " combined distance plot already created"))
    }
  }
}

# executing iteration for the long file!
plot <- Map(Distribution_combDistance, long$Specieslist, long$name)


