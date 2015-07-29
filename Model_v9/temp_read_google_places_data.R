read_googleplaces_data <- function() {
  file_in <- paste0(dropbox_dir,"/../VelibData/ParisData/GooglePlacesData/points_grid_tokens.csv")
  places_data <- read.csv(file_in)
  for(i in c(1:ncol(places_data))) {
    na_list <- which(is.na(places_data[,i]))
    places_data[na_list,i] <- 0
  }
  #assign points to tracts
  wdcCenTrParsed <- readtractboundaryfile() 
  
  #find distance from each points row to each point in census and assign one that 
  #is closest
  places_data$tract <- assign_points_tract(places_data$lat,places_data$lon,wdcCenTrParsed)
  
  return(places_data)
}

