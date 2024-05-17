#Load necessary packages
pacman::p_load(dplyr, tidyverse, phyloseqGraphTest, rstatix, ggpubr, vegan, ggplot2)

### phyloseq doesn't exist anymore  => but phyloseqGraphTest does

#Set wd to wherever all files are stored
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Rarefaction curves ####
#Plot rarefaction results
Read_Count_Species_ARMS <- read.csv("Output/Read_Count_Species_ARMS_10.csv")
OTU_Table_Rar <- Read_Count_Species_ARMS
OTU_Table_Rar <- column_to_rownames(OTU_Table_Rar, "Specieslist")
OTU_Table_Rar <- as.matrix(t(OTU_Table_Rar))

Rarefaction_df <- rarecurve(OTU_Table_Rar, step = 20,
                            col = "blue", 
                            cex = 0.6, tidy = T)
ggplot(Rarefaction_df, aes(x = Sample, y = Species, group = Site)) +
  geom_line(color = "blue", alpha = 0.5) + 
  #geom_vline(xintercept = 100, linetype = "dashed", color = "black", alpha = 0.6, size = 0.6) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = seq(0, 80, 20)) +
  scale_x_continuous(expand = c(0, 500)) +
  ylab("Number of species") +
  xlab("Number of sequences")


#### Alien Read Fraction and Alien Species Fraction ####
GBIF_Distance <- read.csv("Output/DistanceOverSea.csv")
#Remove duplicates
GBIF_Distance <- unique(GBIF_Distance)
#Convert distance to km
GBIF_Distance$distance <- GBIF_Distance$distance/1000
#Assign alien status
Alien_Distance_Threshold <- 250
#Most NAs are assigned due to a bug in the ShortestPath function, where the sampling location and GBIF occurrence are too close to be measured.
  #Therefore, change al NAs to the maximum potential error from the ShortestPath and assume these are native species
GBIF_Distance$distance[is.na(GBIF_Distance$distance)] <- 16
GBIF_Distance$Alien_status <- ifelse(GBIF_Distance$distance > Alien_Distance_Threshold, "YES", "NO")
rm(Alien_Distance_Threshold)
#Remove lines where distance is not calculated
GBIF_Distance <- GBIF_Distance[complete.cases(GBIF_Distance), ]
#Create df with only aliens
Only_Aliens <- GBIF_Distance %>% filter(Alien_status == "YES") %>% select(Specieslist, name, value, distance)
Only_Natives <- GBIF_Distance %>% filter(Alien_status == "NO") %>% select(Specieslist, name, value, distance)
rm(GBIF_Distance)
#Create factor to be able to assign read count and metadata to alien status
Alien_Species_Obs <- paste0(Only_Aliens$Specieslist, "_", Only_Aliens$name)

#Spread df
Only_Aliens <- spread(data = Only_Aliens, key = name, value = value)
#Prepare main dataframe with read counts to assign alien status
RC_Species_Obs_Year <- Read_Count_Species_ARMS %>% ungroup() %>% select(Specieslist:ncol(Read_Count_Species_ARMS))
RC_Species_Obs_Year <- RC_Species_Obs_Year %>% gather(key = "No_Fraction", value = "count", 2:ncol(RC_Species_Obs_Year))
MetaData <- read.csv("Inputs/MetaData_Adjusted.csv")
Meta_RC <- MetaData %>% select(-Fraction, -Filename)
Meta_RC <- left_join(RC_Species_Obs_Year, Meta_RC, by = "No_Fraction", relationship = "many-to-many")
### changed upper line to relationship = "many-to-many"
Meta_RC$Species_Obs <- paste0(Meta_RC$Specieslist, "_", Meta_RC$Observatory.ID)
#Filter main dataframe to one with only alien species
Meta_RC_Alien <- Meta_RC %>% filter(Species_Obs %in% Alien_Species_Obs)
#Remove dupes in order to spread
Meta_RC <- Meta_RC[!duplicated(Meta_RC),]
Meta_RC_Alien <- Meta_RC_Alien[!duplicated(Meta_RC_Alien),]
#Calculate Alien Species Fraction
SpeciesCount <- Meta_RC %>% 
  filter(count > 0) %>%
  group_by(No_Fraction) %>%
  dplyr::summarise(Non_Alien = n_distinct(Specieslist))
AlienSpeciesCount <- Meta_RC_Alien %>% 
  filter(count > 0) %>%
  group_by(No_Fraction) %>%
  dplyr::summarise(Alien = n_distinct(Specieslist))
AlienSpeciesFraction <- left_join(SpeciesCount, AlienSpeciesCount, by = "No_Fraction")
AlienSpeciesFraction$AlienSpeciesFraction <- AlienSpeciesFraction$Alien / (AlienSpeciesFraction$Non_Alien + AlienSpeciesFraction$Alien)
rm(SpeciesCount, AlienSpeciesCount)
#Spread dfs
Meta_RC <- Meta_RC %>% select(-Species_Obs) %>% spread(Specieslist, count)
Meta_RC_Alien <- Meta_RC_Alien %>% select(-Species_Obs) %>% spread(Specieslist, count)

#Fill empty cells with 0
Meta_RC[is.na(Meta_RC)] <- 0
Meta_RC_Alien[is.na(Meta_RC_Alien)] <- 0
##################################################################### UNTIL HERE
#Calculate Alien Read Fraction
TotalReadCount <- rowSums(Meta_RC[, 19:ncol(Meta_RC)])
AlienReadCount <- rowSums(Meta_RC_Alien[, 19:ncol(Meta_RC_Alien)])  ### ERROR
AlienReadFraction <- AlienReadCount/TotalReadCount
Meta_RC$AlienReadFraction <- AlienReadFraction
Meta_RC <- left_join(Meta_RC, AlienSpeciesFraction, by = "No_Fraction")
Meta_RC <- Meta_RC %>%
  select(AlienReadFraction, AlienSpeciesFraction, Non_Alien, Alien, everything())
rm(TotalReadCount, AlienReadCount)

#Make dataframe which filters on at least 100 summed reads per ARMS
Meta_RC_100 <- Meta_RC
Meta_RC_100 <- Meta_RC_100[complete.cases(Meta_RC_100), ]
Meta_RC_100 <- filter(Meta_RC_100, rowSums(Meta_RC_100[,19:ncol(Meta_RC_100)]) > 100)

##### NMDS plot BOLDigger ARMS #####
Meta_RC[is.na(Meta_RC)] <- 0
Meta_RC_Alien[is.na(Meta_RC_Alien)] <- 0
#Only keep complete cases
Meta_RC_1000 <- Meta_RC
Meta_RC_1000 <- Meta_RC_1000[complete.cases(Meta_RC_1000), ]
#Filter on RC over 1000
Minimum_RC <- 1000
Meta_RC_1000 <- filter(Meta_RC_1000, rowSums(Meta_RC_1000[,19:ncol(Meta_RC_1000)]) > Minimum_RC)
#Select Species composition and environmental variables
Com <- as.data.frame(Meta_RC_100[, 19:ncol(Meta_RC_100)])
Env <- Meta_RC_100 %>% select(AlienReadFraction, AlienSpeciesFraction, Country, Observatory.ID, Latitude, Longitude, Depth_m, Monitoring_area, Year, Sample_region_country)
Com <- sapply(Com, as.numeric)

#Calculate nmds
nmds <- metaMDS(Com, distance = "bray")
en <- envfit(nmds, Env, permutations = 999, na.rm = T)
#Adjust dataframe for nicer plot
data.scores <- scores(nmds)
data.scores <- as.data.frame(data.scores$sites)
#Select one environmental variable to focus on
data.scores$Sample_region_country = Env$Sample_region_country
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) *3
en_coord_cat = as.data.frame(scores(en, "factors")) 
palette <- pals::cols25(length(unique(data.scores$Sample_region_country)))
#Plot NMDS focussing on the chosen variable
ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = Sample_region_country), size = 5, alpha = 0.7) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Sample region (country)") +
  scale_color_manual(values = palette)


##### STATISTICAL ANALYSES AND VISUALISATION ALIEN READ/SPECIES FRACTION #####
Comparisons_Mon <- list(c("Industrial port", "Marina"), c("Industrial port", "LHI"), c("Industrial port", "MPA"), 
                        c("Marina", "LHI"), c("Marina", "MPA"),
                        c("LHI", "MPA"))
ggline(data = Meta_RC_1000, x = "Monitoring_area", y = "AlienSpeciesFraction", add = c("mean_se", "jitter")) +
  stat_compare_means(ref.group = "MPA")

ggplot(Meta_RC_1000, aes(Monitoring_area, AlienSpeciesFraction, color = Country)) + 
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(ref.group = "Marina") +
  theme_pubr() +
  scale_color_manual(values = palette)

##### Alien species versus native #####
#Measure difference of AlienSpeciesFraction between Monitoring_area
stat.test <- wilcox_test(Meta_RC, AlienSpeciesFraction ~ Monitoring_area) %>%
  add_xy_position()
#plot
ggboxplot(Meta_RC, x = "Monitoring_area", y = "AlienSpeciesFraction", palette = palette, order = c(unique(stat.test$group1), "MPA"), ) +
  stat_pvalue_manual(stat.test, hide.ns = F) +
  scale_y_continuous(expand = c(0, 0.01))

#Stretch df to fit barplot
Barplot_df <- Meta_RC_100 %>% 
  select(No_Fraction, Non_Alien, Alien, Monitoring_area, ) %>% 
  arrange(Monitoring_area, Non_Alien)
#Save PlotOrder
Barplot_df <- Barplot_df %>% arrange(factor(Monitoring_area, levels = c("Industrial port", "Marina", "LHI", "MPA")))
PlotOrderMA <- Barplot_df$Monitoring_area
Barplot_df <- pivot_longer(data = Barplot_df, cols = c(Non_Alien, Alien), names_to = "Species_status")
#Create the plot
ggplot(Barplot_df, aes(x = No_Fraction, y = value, fill = Species_status)) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_x_discrete(labels = PlotOrderMA) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2)
        #, text = element_text(size = 30)
  ) +
  aes(x = forcats::fct_inorder(No_Fraction)) +
  xlab("Monitoring area") +
  ylab("Number of species") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = seq(0, 80, 20)) +
  labs(fill='Species status') +
  scale_fill_manual(labels = c("Alien", "Native"), 
                    values = c("grey60", "grey40"))


##### Phyloseq #####
#OTU Table
library("tidyverse")
Read_Count_Species_Fraction <- read.csv("Output/Read_Count_Species_Fraction.csv")
OTU_Table <- Read_Count_Species_Fraction
OTU_Table <- column_to_rownames(OTU_Table, "Specieslist")
OTU_Table <- as.matrix(t(OTU_Table))

MetaPhylo <- MetaData
rownames(MetaData) <- MetaData$Filename

#Taxonomic info
Data_LOI <- read.csv("Output/Data_LOI.csv")
Tax_Info <- Data_LOI %>% 
  select(Phylum:Specieslist, -Species) %>%
  distinct() 
Tax_Info <- Tax_Info[!duplicated(Tax_Info$Specieslist),]
rownames(Tax_Info) <- Tax_Info$Specieslist
Tax_Info <- as.matrix(Tax_Info)

#Combine into Phyloseq object
physeq <- phyloseq(otu_table(OTU_Table, taxa_are_rows = FALSE), 
                   tax_table(Tax_Info), 
                   sample_data(MetaData))

#Add Plot Order
PlotOrder_phy <- MetaData %>% 
  arrange(Sample_region_country, Observatory.ID, Year) %>%
  select(Filename) %>%
  distinct()
PlotOrder_phy <- as.list(PlotOrder_phy)
sample_data(physeq)$NewID <- factor(sample_names(physeq))
sample_data(physeq)$NewID <- factor(sample_data(physeq)$NewID, levels = PlotOrder_phy$Filename)

#Remove samples with less than 100 reads
physeq <- prune_samples(sample_sums(physeq) > 100, physeq)

#Merge samples based on variable by summing read counts
palette_SRC <- pals::cols25(length(unique(Meta_RC$Sample_region_country)))
#Merge by Sample_region_country...
mergedphyseq_SRC = merge_samples(physeq, "Sample_region_country")
plot_heatmap(mergedphyseq_SRC, taxa.order = "Phylum", max.label = 500, taxa.label = "Phylum")
#... or Monitoring_area
mergedphyseq_MA = merge_samples(physeq, "Monitoring_area")
plot_heatmap(mergedphyseq_MA, taxa.order = "Phylum", max.label = 500)

#### Diversity of native species ####
#Plot diversity
Richness <- estimate_richness(physeq, measures = c("Observed", "Chao1", "Shannon"))
Richness <- left_join(rownames_to_column(Richness), rownames_to_column(MetaData), by = "rowname")
Richness <- Richness %>% arrange(factor(Monitoring_area, levels = c("Industrial port", "Marina", "LHI", "MPA")))
#Measure difference of Observed between Monitoring_area
stat.test <- wilcox_test(Richness, Observed ~ Monitoring_area) %>%
  add_xy_position()
#plot diversity of entire species pool
ggboxplot(Richness, x = "Monitoring_area", y = "Observed", palette = palette) +
  stat_pvalue_manual(stat.test, hide.ns = T, label = "p.adj", size = 10) +
  theme(text = element_text(size = 30))

#plot diversity of just native species
stat.test <- wilcox_test(Meta_RC_100, Non_Alien ~ Monitoring_area) %>%
  add_xy_position()
ggboxplot(Meta_RC_100, x = "Monitoring_area", y = "Non_Alien", palette = palette) +
  stat_pvalue_manual(stat.test, hide.ns = T, label = "p.adj", size = 10)

#### Alien phyloseq prep ####
Aliens_OTU <- Read_Count_Species_Fraction
Aliens_OTU <- Aliens_OTU %>% filter(Specieslist %in% unique(Only_Aliens$Specieslist))
Aliens_OTU <- column_to_rownames(Aliens_OTU, "Specieslist")
Aliens_OTU[Aliens_OTU > 0] <- 1
Aliens_OTU <- as.matrix(Aliens_OTU)

#Add Alien_Type to Taxonomic info
As.the.crow.flies <- read.csv("Output/ShortestPath.csv")
#Extract observation count, to re-add after gather function
Observation_Count <- select(As.the.crow.flies, c("Speciesname", "Total.observations", "Unique.locations"))
As.the.crow.flies <- As.the.crow.flies %>% 
  select(-Total.observations, -Unique.locations, -Errors) %>% 
  gather(key = "Observatory.ID", value = "distance_atcf", 2:13)
As.the.crow.flies <- left_join(As.the.crow.flies, Observation_Count, by = "Speciesname")
#Convert meters to kilometres for easier reading
As.the.crow.flies$distance_atcf <- as.numeric(As.the.crow.flies$distance_atcf)
As.the.crow.flies$distance_atcf <- ceiling(As.the.crow.flies$distance_atcf/1000)
#Find species with no observation on GBIF
Cryptogenic <- As.the.crow.flies[is.na(As.the.crow.flies$Total.observations),]
#Calculate mean distance to closes observation in order to determine Alien_Type
Mean <- As.the.crow.flies %>% na.omit %>% group_by(Speciesname) %>% summarise_at(vars(distance_atcf), list(name = mean))
As.the.crow.flies <- left_join(As.the.crow.flies, Mean, by = "Speciesname")
As.the.crow.flies <- As.the.crow.flies[As.the.crow.flies$Speciesname %in% Only_Aliens$Specieslist,]
As.the.crow.flies$Alien_Type <- ifelse(As.the.crow.flies$Total.observations < 50 | is.na(As.the.crow.flies$Total.observations), As.the.crow.flies$Alien_Type <- "Cryptogenic", 
                                       ifelse(As.the.crow.flies$name > 2500, As.the.crow.flies$Alien_Type <- "Hitchhiker", "Range expander"
                                       )
)
#DF with Alien types
Alien_Type <- As.the.crow.flies %>% select(Speciesname, Alien_Type) %>% unique() %>% arrange(Speciesname)
#Add species with no gbif occurences
Cryptogenic <- as.data.frame(unique(Cryptogenic$Speciesname))
Cryptogenic$Alien_Type <- "Cryptogenic"
colnames(Cryptogenic) <- c("Speciesname", "Alien_Type")

# Plot native biodiversity and alien species fraction
library("SciViews")
palette_4 <- c("#440154", "#12c5ed", "#fde725", "#0fb859")
ggplot(Meta_RC_100, aes(Non_Alien, AlienSpeciesFraction, color = Monitoring_area)) +
  geom_point(size = 8, alpha = 0.7) +
  stat_smooth(method = lm, formula = y ~ ln(x), inherit.aes = F, data = Meta_RC_100, mapping = aes(Non_Alien, AlienSpeciesFraction)) +
  stat_regline_equation(formula = y ~ ln(x), inherit.aes = F, data = Meta_RC_100, 
                        mapping = aes(Non_Alien, AlienSpeciesFraction, label = ..adj.rr.label..), size = 8, label.x = 60, label.y = 0.33) +
  stat_regline_equation(formula = y ~ ln(x), inherit.aes = F, data = Meta_RC_100, 
                        mapping = aes(Non_Alien, AlienSpeciesFraction, label = ..eq.label..), size = 8, label.x = 60, label.y = 0.37) +
  theme_pubclean(base_size = 30) +
  scale_color_manual(values = palette_4) +
  xlab("Species with confirmed presence") +
  ylab("New alien species fraction")

#Add Alien_type to Taxonomix info, prepare for phyloseq
Aliens_Tax <- as.data.frame(Tax_Info)
Aliens_Tax <- Aliens_Tax %>% filter(Specieslist %in% Alien_Type$Speciesname)
Aliens_Tax <- left_join(Aliens_Tax, Alien_Type, by = c("Specieslist" = "Speciesname"))
rownames(Aliens_Tax) <- Aliens_Tax$Specieslist
Alien_Order <- Aliens_Tax %>% arrange(Alien_Type) %>% select(Specieslist)
Aliens_Tax <- as.matrix(Aliens_Tax)

#### Phyloseq of aliens ####
Alien_physeq <- phyloseq(otu_table(Aliens_OTU, taxa_are_rows = TRUE), 
                         tax_table(Aliens_Tax), 
                         sample_data(MetaData))

#Add plotorder of samples
PlotOrder_Al <- MetaData %>% arrange(factor(Monitoring_area, levels = c("Industrial Port", "Marina", "LHI", "MPA"))) %>% select(Filename)
sample_data(Alien_physeq)$NewID <- factor(sample_names(Alien_physeq))
sample_data(Alien_physeq)$NewID <- factor(sample_data(Alien_physeq)$NewID, levels = PlotOrder_Al$Filename)

#Remove samples with less than 100 alien reads
Alien_physeq <- prune_samples(sample_sums(Alien_physeq) > 100, Alien_physeq)

#Plot heatmap of SRC
Alienphyseq_SRC = merge_samples(Alien_physeq, "Sample_region_country")
plot_heatmap(Alienphyseq_SRC, taxa.order = Alien_Order$Specieslist, taxa.label = "Alien_Type", )

#Plot heatmap of MA
Alienphyseq_MA = merge_samples(Alien_physeq, "Monitoring_area")
plot_heatmap(Alienphyseq_MA, taxa.order = Alien_Order$Specieslist, taxa.label = "Alien_Type")

#Merge Alien types
Range_expander <- Alien_Type[Alien_Type$Alien_Type == "Range expander", ]
Hitchhiker <- Alien_Type[Alien_Type$Alien_Type == "Hitchhiker" ,]
Cryptogenic <- Alien_Type[Alien_Type$Alien_Type == "Cryptogenic" ,]

#Merge Alien types and Monitoring areas
x1 = merge_taxa(Alienphyseq_MA, Range_expander$Speciesname, 2)
x2 = merge_taxa(x1, Hitchhiker$Speciesname, 2)
x3 = merge_taxa(x2, Cryptogenic$Speciesname, 2)
#plot heatmap of merged monitoring type and merged alien type
plot_heatmap(x3, taxa.order = c("Acartia hudsonica","Cephalothrix spiralis", "Amblyosyllis madeirensis"), sample.order = c("Industrial port", "Marina", "LHI", "MPA")) +
  scale_y_discrete(labels = c("Range expanders", "Hithchikers", "Cryptogenic"))

#Statistics on Alien type
Alien_Type_Stats <- as.data.frame(otu_table(x3))
colnames(Alien_Type_Stats) <- c("Range expanders", "Cryptogenic", "Hitchhikers")
chisq.test(Alien_Type_Stats)
Alien_Type_Stats <- rownames_to_column(Alien_Type_Stats, "Monitoring_area")
Alien_Type_Stats <- pivot_longer(Alien_Type_Stats, cols = c("Range expanders", "Cryptogenic", "Hitchhikers"))


