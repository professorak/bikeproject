file_in <- "/Users/ashishkabra/Dropbox/VelibData/ParisData/GooglePlacesData/points_grid_tokens.csv"
places_data <- read.csv(file_in)
for(i in c(1:ncol(places_data))) {
  na_list <- which(is.na(places_data[,i]))
  places_data[na_list,i] <- 0
}

#choosing_tokens-- point_of_interest, museum, food, lodging, grocery_or_supermarket, 
table(places_data$point_of_interest)
table(places_data$museum)
table(places_data$food)
table(places_data$lodging)
table(places_data$grocery_or_supermarket)