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
species_location <- read.csv("Output/1_Species_Location.csv")
long <- pivot_longer(species_location, !Specieslist)
# from tidyr package, reshape df from wide to long format
long <- long[long$value > 0, ]

####################################################

# iterate over species_location file again using function and map()
# and read csv file per species

species_name <- "Acartia bifilosa"
species_location <- "Koster"

Distribution_seadistance <- function(species_name, species_location){
  # read csv file per species
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  distance_file <- paste0("theoretical_data/sea_distances/", species_name, "_distancesTo_", species_location, "_realData.csv")
  
  distances <- read.table(distance_file, header = TRUE)
  # clean dataframe from rows with Inf in them
  distances <- distances[is.finite(distances$x), ]
  distances$x <- distances$x/1000
  distances <- subset(distances, x < 40000)
  
  ###############################################################################
  # make histograms of distances per species, with filtering on distance limit 40000
  ###############################################################################
  
  dist_plot <- ggplot(distances, aes(x = x)) +
    geom_histogram(binwidth = 50, fill = "blue", color = "black", boundary = 0) +
    labs(title = "Histogram of Distances", x = "Distance in km", y = "Frequency of species") +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +
    ggtitle(paste0("Distribution of ", species_name, " from ", species_location)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20)) +         # Set font size for axis titles
    scale_x_continuous(breaks = seq(0, 3000, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 3000)) # Use coord_cartesian for setting limits
    
  print(dist_plot)
  # create a directory to save the plots in
  
  if(!dir.exists("test_outputs/sea_distribution_plots")) {
    dir.create("test_outputs/sea_distribution_plots")
    ggsave(filename = paste0("test_outputs/sea_distribution_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = dist_plot, width = 50, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("test_outputs/sea_distribution_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = dist_plot, width = 50, height = 30, units = "cm", dpi = 900)
  }
}
################################# END FUNCTION

# error lists initiations to store error messages in
error <- c()

# executing iteration for the long file!
plot <- Map(Distribution_seadistance, long$Specieslist, long$name)

# write error files
print(error)

file_error <- paste0("test_outputs/errors_graph_sea_distances.csv",)
write.table(error, file = file_error, append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)



#############################################################
# Make combined histogram of sea distances and fly distances
#############################################################
species_name <- "Acartia bifilosa"
species_location <- "Koster"

Distribution_combDistance <- function(species_name, species_location){
  sea_distance_file <- paste0("theoretical_data/sea_distances/", species_name, "_distancesTo_", species_location, "_realData.csv")
  fly_distance_file <- paste0("theoretical_data/fly_distances/", species_name, "_distancesTo_", species_location, "_realData.csv")
  
  sea_distances <- read.table(sea_distance_file)
  fly_distances <- read.table(fly_distance_file)
  # transform Inf to NA
  sea_distances$x[is.infinite(sea_distances$x)] <- NA
  
  
  sea_distances$x <- sea_distances$x/1000
  sea_distances <- subset(sea_distances, x < 40000)
  sea_distances$type <- "sea_distance"
  fly_distances$x <- fly_distances$x/1000
  fly_distances <- subset(fly_distances, x < 40000)
  fly_distances$type <- "fly_distance"
  
  combined_distances <- rbind(sea_distances, fly_distances)
  
  p <- combined_distances %>%
    ggplot( aes(x=x, fill=type)) +
    geom_histogram(binwidth = 50, color="#e9ecef", alpha=0.6, position = 'identity') +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +
    labs(x = "Distance in km", y = "Frequency") +
    ggtitle(paste0("Distribution of ", species_name, " from ", species_location)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
      plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
      axis.text = element_text(size = 16),  # Set font size for axis numbers
      axis.title = element_text(size = 20),
      legend.title = element_text(size = 18, face="bold"),
      legend.text = element_text(size = 16)) +
    scale_x_continuous(breaks = seq(0, 3000, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 3000)) + # Use coord_cartesian for setting limits
    # set legend title and labels
    scale_fill_discrete(
      name = "Distance type",
      breaks = c("sea_distance", "fly_distance"),
      labels = c("Sea distance", "Fly distance")) +
    labs(x = "Distance in km", y = "Frequency")
  
  print(p)
  
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
species_name <- "Acartia bifilosa"
# Alcyonium digitatum
species_location <- "Koster"

Location_histograms <- function(species_name, species_location){
  # read csv file per species
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  distance_file <- paste0("test_outputs/sea_distances/", species_name, "_distancesTo_", species_location, "_realData.csv")
  
  distance_file <- read.table(distance_file, header = TRUE, sep = ",")
  distance_file$x <- distance_file$x/1000
  
  distances <- subset(distance_file, x < 40000)
  countries <- distances$country
  #############################################################
  
  ### PLOT
  
  country_plot <- ggplot(distances, aes(x = x, fill = country)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack") +  # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances/country for ", species_name," in ", species_location),
                        x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +  # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # Increase legend title size
          legend.text = element_text(size = 12),    # Increase legend text size
          legend.key.size = unit(1.5, "lines")) +   # Increase legend key size
    scale_x_continuous(breaks = seq(0, 3000, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 3000)) # Use coord_cartesian for setting limits
  
  print(country_plot)
  # create a directory to save the plots in
  
  if(!dir.exists("test_outputs/sea_distribution_country_plots")) {
    dir.create("test_outputs/sea_distribution_country_plots")
    ggsave(filename = paste0("theoretical_data/plot_", species_name, "_from_", species_location, ".png"), 
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
species_name <- "Acartia bifilosa"
# Alcyonium digitatum
species_location <- "Roscoff"

Year_histograms <- function(species_name, species_location){
  # read csv file per species
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  distance_file <- paste0("theoretical_data/sea_distances/", species_name, "_distancesTo_", species_location, "_realData.csv")
  
  distance_file <- read.table(distance_file, header = TRUE)
  distance_file$x <- distance_file$x/1000
  
  distances <- subset(distance_file, x < 40000)
  years <- distances$year
  #############################################################
  assign_category <- function(year) {
    if (is.na(year)) {
      return(NA)
    }
    for (category in year_categories) {
      range <- as.numeric(unlist(strsplit(category, "-")))
      if (year >= range[1] & year < range[2]) {
        return(category)
      }
    }
    return(NA) # If year doesn't fall into any category
  }
  year_categories <- c("1965-1985","1985-1990",
                       "1990-1995","1995-2000","2000-2005","2005-2010","2010-2015",
                       "2015-2020","2020-2025")
  # Apply function to create new column
  distances$year_category <- sapply(distances$year, assign_category)
  
  #################################
  ### PLOT
  
  year_plot <- ggplot(distances, aes(x = x, fill = year_category)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack") +  # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances/year for ", species_name," in ", species_location),
         x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    scale_fill_brewer(palette = "YlOrRd", na.value = "black") + # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # Increase legend title size
          legend.text = element_text(size = 12),    # Increase legend text size
          legend.key.size = unit(1.5, "lines")) +   # Increase legend key size
    scale_x_continuous(breaks = seq(0, 3000, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 3000)) # Use coord_cartesian for setting limits
  
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
plot <- Map(Year_histograms, long$Specieslist, long$name)

