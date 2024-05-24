library("ggplot2")

#Input data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
straightLine <- read.csv("Output/ShortestPath.csv")
throughSea <- read.csv("Output/DistanceOverSea.csv")


test <- merge(throughSea, straightLine, by.y = "Speciesname", by.x = "Specieslist")
if(!is.factor(test$name)) {
  test$name <- factor(test$name)
}

print(levels(test$name))
for(item in levels(test$name)){  # levels() is used to retrieve levels from test$name
  test$straight[test$name == item] <- test[test$name == item, item] 
}
# test$straight[test$name == item] => selects straight column values for those rows
# test[test$name == item, item] => extracts the values from the column named 'item' in 'test' for those rows
# complete line: updates 'straight' column for selected rows with corresponding
# values from the column named 'item'

#remove unused sampling locations
test <- test[,!names(test) %in% c(levels(test$name), "Galway", "Gdynia")]

#plot
ggplot(test)+
  geom_point(aes(x = straight , y = distance), col = "green") +
  geom_abline()

ggplot(test) + 
  geom_density(aes(x = log10(straight), fill = is.na(distance)), alpha = 0.5)
