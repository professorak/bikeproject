read_googleplaces_data <- function() {
  file_in <- paste0(dropbox_dir,"/../VelibData/ParisData/GooglePlacesData/points_grid_tokens.csv")
  places_data <- read.csv(file_in)
  for(i in c(1:ncol(places_data))) {
    na_list <- which(is.na(places_data[,i]))
    places_data[na_list,i] <- 0
  }
  return(places_data)
}




# #computing coef. of variation for each column
# colmeans <- colMeans(places_data)
# sd <- apply(places_data, 2, FUN=sd)
# coefvar <- sd/colmeans
# coefvar <- coefvar[-c(1:2)]
# 
# coefvar <- coefvar[order(coefvar, decreasing=T)]
# 
# summary(coefvar)


