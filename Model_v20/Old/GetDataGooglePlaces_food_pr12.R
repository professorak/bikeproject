source("GetDataGooglePlaces_pr12.R")
if(is.null(points$food)) {
  stop("food variable doesnt exist in points")
}
points$places_count <- points$food
