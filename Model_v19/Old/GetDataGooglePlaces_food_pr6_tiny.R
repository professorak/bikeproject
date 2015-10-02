source("GetDataGooglePlaces_pr6_tiny.R")

if(is.null(points$food)) {
  stop("food variable doesnt exist in points")
}
points$places_count <- points$food
