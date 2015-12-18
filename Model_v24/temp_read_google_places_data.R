read_googleplaces_data <- function() {
  file_in <- paste0(dropbox_dir,"/../VelibData/ParisData/GooglePlacesData/points_grid_tokens.csv")
  places_data <- read.csv(file_in)
  for(i in c(1:ncol(places_data))) {
    na_list <- which(is.na(places_data[,i]))
    if(length(na_list)) {
      places_data[na_list,i] <- 0
    }    
  }
  #assign points to tracts
  wdcCenTrParsed <- readtractboundaryfile() 
  
  #find distance from each points row to each point in census and assign one that 
  #is closest
  places_data$tract <- assign_points_tract(places_data$lat,places_data$lon,wdcCenTrParsed)
  
  return(places_data)
}

read_googleplaces_data_dis_points <- function() {
  if(identical(tract,c(1:10))) {
    file_in <- paste0(dropbox_dir,"/../VelibData/ParisData/GooglePlacesData/points_grid_tokens_dis_1to10_dispoints",dis_points,".csv")
  } else if(identical(tract,c(1:20))) {
    file_in <- paste0(dropbox_dir,"/../VelibData/ParisData/GooglePlacesData/points_grid_tokens_dis_1to20_dispoints",dis_points,".csv")
  } else {
    file_in <- paste0(dropbox_dir,"/../VelibData/ParisData/GooglePlacesData/points_grid_tokens_dis_",
      paste0(tract,collapse = "_"),"_dispoints",dis_points,".csv")
  }
  
  
  places_data <- read.csv(file_in)
  for(i in c(1:ncol(places_data))) {
    na_list <- which(is.na(places_data[,i]))
    if(length(na_list)) {
      places_data[na_list,i] <- 0
    }    
  }
  #assign points to tracts
  wdcCenTrParsed <- readtractboundaryfile() 
  
  #find distance from each points row to each point in census and assign one that 
  #is closest
  places_data$tract <- assign_points_tract(places_data$lat,places_data$lon,wdcCenTrParsed)
  
  return(places_data)
}

