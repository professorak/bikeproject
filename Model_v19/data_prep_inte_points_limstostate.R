# this file codes generation of points used for integration
# crude way of generating \del lat, lon for points are certain frequency

generate_integration_points_in <- function(wdcMerged_in=wdcMerged) {
  
  #gen integration boundaries 
  if(!exists("min_lon") | !exists("max_lon") | !exists("min_lat") | !exists("max_lat")) {
    min_lon <- min(wdcMerged_in$lon)
    max_lon <- max(wdcMerged_in$lon)
    min_lat <- min(wdcMerged_in$lat)
    max_lat <- max(wdcMerged_in$lat)    
  }
  
  #distances in kms
  
  R <- 6371 #earths radius
  raddeg <- 180/3.14
  dis_lat <- dis_points/R*raddeg
  dis_lon <- dis_points/R*raddeg/abs(cos(ave(min_lon, max_lon)))
  
  #boundaries of integration region
  r_min_lat <- min_lat - 2*max_walking_dis/dis_points*dis_lat
  r_max_lat <- max_lat + 2*max_walking_dis/dis_points*dis_lat
  r_min_lon <- min_lon - 2*max_walking_dis/dis_points*dis_lon
  r_max_lon <- max_lon + 2*max_walking_dis/dis_points*dis_lon
  
  lat_seq <- seq(r_min_lat, r_max_lat, by=dis_lat)
  lon_seq <- seq(r_min_lon, r_max_lon, by=dis_lon)
  
  points <- merge(lat_seq,lon_seq)
  colnames(points) <- c("lat","lon")
  #points$density <- 0.1
  points$type <- 1

  #keep only points which are at distance less than max distance from nearest station
  wdcTemp <- unique(wdcMerged_in[,c("lat","lon")])
  lat2 <- wdcTemp[,"lat"]
  lon2 <- wdcTemp[,"lon"]
  
  keep_index_func <- function(t) {
    dis_v <- latlondistance(t[1], t[2], lat2, lon2)    
    if(min(dis_v)<max_points_distance ) return(1)
    return(0)
  }
  points$keep_index <- apply(points[,c("lat","lon")],1,FUN=keep_index_func)
  points <- subset(points, keep_index==1)
  points$keep_index <- NULL
  
  #assign points to tracts
  wdcCenTrParsed <- readtractboundaryfile()
  wdcCenTrDensity <- readtractdensityfile()
  
  #find distance from each points row to each point in census and assign one that 
  #is closest
  points$tract <- assign_points_tract(points$lat,points$lon,wdcCenTrParsed)
  
  #merge points with density data
  points <- merge(points,wdcCenTrDensity,by="tract")
  
  #scale density according to the area around point. 
  #it is a rectangle of dis_points*dis_points*10^6 as distance is in kms and density in /m2
  points$density <- points$density *dis_points*dis_points
  points$density <- points$density /6 #to account for number of time windows in a day
    
  unlink(wdcCenTrParsed)
  unlink(wdcCenTrDensity)
  return(points)
}

generate_integration_points <- function(wdcMerged_in=wdcMerged) {

  points <- generate_integration_points_in(wdcMerged_in)
  points_metro <- generate_metro_density_points()
  points_metro$type <- 2
  points_bus <- generate_bus_density_points()
  points_bus$type <- 3
  points_tourist <- generate_tourist_density_points()
  points_tourist$type <- 4
  
  #keep only points which are at distance less than max distance from nearest station
  wdcTemp <- unique(wdcMerged_in[,c("lat","lon")])
  lat2 <- wdcTemp[,"lat"]
  lon2 <- wdcTemp[,"lon"]
  
  keep_index_func <- function(t) {
    dis_v <- latlondistance(t[1], t[2], lat2, lon2)    
    if(min(dis_v)<max_points_distance ) return(1)
    return(0)
  }

  points_metro$keep_index <- apply(points_metro[,c("lat","lon")],1,FUN=keep_index_func)
  points_metro <- subset(points_metro, keep_index==1)
  points_metro$keep_index <- NULL
  points_bus$keep_index <- apply(points_bus[,c("lat","lon")],1,FUN=keep_index_func)
  points_bus <- subset(points_bus, keep_index==1)
  points_bus$keep_index <- NULL
  points_tourist$keep_index <- apply(points_tourist[,c("lat","lon")],1,FUN=keep_index_func)
  points_tourist <- subset(points_tourist, keep_index==1)
  points_tourist$keep_index <- NULL
  
  #assign points to tracts
  wdcCenTrParsed <- readtractboundaryfile()
  
  points_metro$tract <- assign_points_tract(points_metro$lat,points_metro$lon,wdcCenTrParsed)
  points_bus$tract <- assign_points_tract(points_bus$lat,points_bus$lon,wdcCenTrParsed)
  points_tourist$tract <- assign_points_tract(points_tourist$lat,points_tourist$lon,wdcCenTrParsed)
  
  points <- rbind(points, points_metro, points_bus, points_tourist)

  unlink(wdcCenTrParsed)
  return(points)
}

generate_metro_density_points <- function() {
  metrotrafficloc_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/ParisMetroData/metro_stations_data_curated.csv") 
  metro_traffic_locations <- read.csv(metrotrafficloc_file_name)
#   #keep only metro stations. There are both Metro and RER stations in the file
#   metro_traffic_locations <- subset(metro_traffic_locations, station_type=="Metro")
  metro_traffic_locations <- metro_traffic_locations[,c("lat","lon","traffic")]
  colnames(metro_traffic_locations) <- c("lat","lon","density")
#  unlink(wdcCenTrParsed)
  return(metro_traffic_locations)
}

generate_bus_density_points <- function() {
  busstops_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/RATPOpenData/Accessibility RATP bus stops/bus_stops_data_dis1to10.csv") 
  busstops_locations <- read.csv(busstops_file_name)
  
  busstops_locations <- busstops_locations[,c("Latitude","Longitude")]
  busstops_locations$density <- 1
  colnames(busstops_locations) <- c("lat","lon","density")
  return(busstops_locations)
}

generate_tourist_density_points <- function() {
  tourist_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/ParisTourismeData/tourist_curated_list.csv") 
  tourist_locations <- read.csv(tourist_file_name)
  
  tourist_locations <- tourist_locations[,c("lat","lon","traffic")]
  colnames(tourist_locations) <- c("lat","lon","density")
  return(tourist_locations)
}


generate_v0  <- function(no_vals = 40) {
  set.seed(2)    
  return(rnorm(no_vals)  )
}

