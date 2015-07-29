# this file codes generation of points used for integration
# crude way of generating \del lat, lon for points are certain frequency
generate_integration_points <- function(wdcMerged_in=wdcMerged,  dis_points=0.050) {

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
  points_metro <- generate_metro_density_points()
  points_metro$type <- 2
  points_bus <- generate_bus_density_points()
  points_bus$type <- 3
  #keep only points which are at distance less than max distance from nearest station
  wdcTemp <- unique(wdcMerged_in[,c("lat","lon")])
  lat2 <- wdcTemp[,"lat"]
  lon2 <- wdcTemp[,"lon"]
  
  keep_index_func <- function(t) {
    dis_v <- latlondistance(t[1], t[2], lat2, lon2)    
    if(min(dis_v)<max_walking_dis ) return(1)
    return(0)
  }
  points$keep_index <- apply(points[,c("lat","lon")],1,FUN=keep_index_func)
  points <- subset(points, keep_index==1)
  points$keep_index <- NULL
  points_metro$keep_index <- apply(points_metro[,c("lat","lon")],1,FUN=keep_index_func)
  points_metro <- subset(points_metro, keep_index==1)
  points_metro$keep_index <- NULL
  points_bus$keep_index <- apply(points_bus[,c("lat","lon")],1,FUN=keep_index_func)
  points_bus <- subset(points_bus, keep_index==1)
  points_bus$keep_index <- NULL
  
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

  points_metro$tract <- assign_points_tract(points_metro$lat,points_metro$lon,wdcCenTrParsed)
  points_bus$tract <- assign_points_tract(points_bus$lat,points_bus$lon,wdcCenTrParsed)
  
  points <- rbind(points, points_metro, points_bus)

  unlink(wdcCenTrParsed)
  unlink(wdcCenTrDensity)
  return(points)
}

generate_metro_density_points <- function() {
  metrotrafficloc_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/ParisTransportData/metro_traffic_locations.csv") 
  metro_traffic_locations <- read.csv(metrotrafficloc_file_name)
  
#   #keep only points which are at distance less than max distance from nearest station
#   wdcTemp <- unique(wdcMerged_in[,c("lat","lon")])
#   lat2 <- wdcTemp[,"lat"]
#   lon2 <- wdcTemp[,"lon"]
#   
#   keep_index_func <- function(t) {
#     dis_v <- latlondistance(t[1], t[2], lat2, lon2)    
#     if(min(dis_v)<max_walking_dis ) return(1)
#     return(0)
#   }
#   metro_traffic_locations$keep_index <- apply(metro_traffic_locations[,c("lat","lon")],1,FUN=keep_index_func)
#   metro_traffic_locations <- subset(metro_traffic_locations, keep_index==1)
#   metro_traffic_locations$keep_index <- NULL
#   
#   #assign to tracts
#   wdcCenTrParsed <- readtractboundaryfile()
#   
#   #find distance from each points row to each point in census and assign one that 
#   #is closest
#   metro_traffic_locations$tract <- 
#     assign_points_tract(metro_traffic_locations$lat,metro_traffic_locations$lon,
#                         wdcCenTrParsed)
#   
  metro_traffic_locations <- metro_traffic_locations[,c("lat","lon","Traffic")]
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

# generate_integration_points_2 <- function(wdcMerged, dis_points=0.050) {
#   
#   #gen integration boundaries 
#   min_lon <- min(wdcMerged$lon)
#   max_lon <- max(wdcMerged$lon)
#   min_lat <- min(wdcMerged$lat)
#   max_lat <- max(wdcMerged$lat)
#   
#   #distances in kms
#   
#   R <- 6371 #earths radius
#   raddeg <- 180/3.14
#   dis_lat <- dis_points/R*raddeg
#   dis_lon <- dis_points/R*raddeg/abs(cos(ave(min_lon, max_lon)))
#   
#   #boundaries of integration region
#   r_min_lat <- min_lat - 2*max_walking_dis/dis_points*dis_lat
#   r_max_lat <- max_lat + 2*max_walking_dis/dis_points*dis_lat
#   r_min_lon <- min_lon - 2*max_walking_dis/dis_points*dis_lon
#   r_max_lon <- max_lon + 2*max_walking_dis/dis_points*dis_lon
#   
#   lat_seq <- seq(r_min_lat, r_max_lat, by=dis_lat)
#   lon_seq <- seq(r_min_lon, r_max_lon, by=dis_lon)
#   
#   points <- merge(lat_seq,lon_seq)
#   colnames(points) <- c("lat","lon")
#   
#   #keep only points which are at distance less than max distance from nearest station
#   wdcTemp <- unique(wdcMerged[,c("lat","lon")])
#   lat2 <- wdcTemp[,"lat"]
#   lon2 <- wdcTemp[,"lon"]
#   
#   keep_index_func <- function(t) {
#     dis_v <- latlondistance(t[1], t[2], lat2, lon2)    
#     if(min(dis_v)<max_walking_dis ) return(1)
#     return(0)
#   }
#   points$keep_index <- apply(points[,c("lat","lon")],1,FUN=keep_index_func)
#   points <- subset(points, keep_index==1)
#   points$keep_index <- NULL
#   
#   #assign points to tracts
#   wdcCenTrParsed <- readtractboundaryfile()
#   wdcCenTrDensity <- readtractdensityfile()
#   
#   #find distance from each points row to each point in census and assign one that 
#   #is closest
#   points$tract <- assign_points_tract(points$lat,points$lon,wdcCenTrParsed)
#   
#   #merge points with density data
#   points <- merge(points,wdcCenTrDensity,by="tract")
#   
#   #scale density according to the area around point. 
#   #it is a rectangle of dis_points*dis_points*10^6 as distance is in kms and density in /m2
#   points$density <- points$density *dis_points*dis_points*10^6
#   points$density <- points$density /60 #to get better fixed effects in numerical computation range.
#   
#   #for each point keep on lat lon density
#   #points$tract <- NULL
#   points$Population <- NULL
#   points$ALAND10 <- NULL
#   unlink(wdcCenTrParsed)
#   unlink(wdcCenTrDensity)
#   return(points)
# }


generate_v0  <- function(no_vals = 40) {
  set.seed(2)    
  return(rnorm(no_vals)  )
}

# generate_integration_points_stdata <- function(station_data,  dis_points=0.050) {
#   
#   #gen integration boundaries 
#   min_lon <- min(station_data$lon)
#   max_lon <- max(station_data$lon)
#   min_lat <- min(station_data$lat)
#   max_lat <- max(station_data$lat)
#   
#   #distances in kms
#   
#   R <- 6371 #earths radius
#   raddeg <- 180/3.14
#   dis_lat <- dis_points/R*raddeg
#   dis_lon <- dis_points/R*raddeg/abs(cos(ave(min_lon, max_lon)))
#   
#   #boundaries of integration region
#   r_min_lat <- min_lat - 2*max_walking_dis/dis_points*dis_lat
#   r_max_lat <- max_lat + 2*max_walking_dis/dis_points*dis_lat
#   r_min_lon <- min_lon - 2*max_walking_dis/dis_points*dis_lon
#   r_max_lon <- max_lon + 2*max_walking_dis/dis_points*dis_lon
#   
#   lat_seq <- seq(r_min_lat, r_max_lat, by=dis_lat)
#   lon_seq <- seq(r_min_lon, r_max_lon, by=dis_lon)
#   
#   points <- merge(lat_seq,lon_seq)
#   colnames(points) <- c("lat","lon")
#   
#   #keep only points which are at distance less than max distance from nearest station
#   wdcTemp <- unique(station_data[,c("lat","lon")])
#   lat2 <- wdcTemp[,"lat"]
#   lon2 <- wdcTemp[,"lon"]
#   
#   keep_index_func <- function(t) {
#     dis_v <- latlondistance(t[1], t[2], lat2, lon2)    
#     if(min(dis_v)<max_walking_dis ) return(1)
#     return(0)
#   }
#   points$keep_index <- apply(points[,c("lat","lon")],1,FUN=keep_index_func)
#   points <- subset(points, keep_index==1)
#   points$keep_index <- NULL
#   
#   #assign points to tracts
#   wdcCenTrParsed <- readtractboundaryfile()
#   wdcCenTrDensity <- readtractdensityfile()
#   
#   #find distance from each points row to each point in census and assign one that 
#   #is closest
#   points$tract <- assign_points_tract(points$lat,points$lon,wdcCenTrParsed)
#   
#   #merge points with density data
#   points <- merge(points,wdcCenTrDensity,by="tract")
#   
#   #scale density according to the area around point. 
#   #it is a rectangle of dis_points*dis_points*10^6 as distance is in kms and density in /m2
#   points$density <- points$density *dis_points*dis_points*10^6
#   points$density <- points$density /60 #to get better fixed effects in numerical computation range.
#   
#   #for each point keep on lat lon density
#   #points$tract <- NULL
#   points$Population <- NULL
#   points$ALAND10 <- NULL
#   unlink(wdcCenTrParsed)
#   unlink(wdcCenTrDensity)
#   return(points)
# }


# generate_integration_points_cluster <- function(  dis_points=0.050) {
#   
#   #gen integration boundaries 
#   min_lon <- min(wdcMerged$lon)
#   max_lon <- max(wdcMerged$lon)
#   min_lat <- min(wdcMerged$lat)
#   max_lat <- max(wdcMerged$lat)
#   
#   #distances in kms
#   
#   R <- 6371 #earths radius
#   raddeg <- 180/3.14
#   dis_lat <- dis_points/R*raddeg
#   dis_lon <- dis_points/R*raddeg/abs(cos(ave(min_lon, max_lon)))
#   
#   #boundaries of integration region
#   r_min_lat <- min_lat - 2*max_walking_dis/dis_points*dis_lat
#   r_max_lat <- max_lat + 2*max_walking_dis/dis_points*dis_lat
#   r_min_lon <- min_lon - 2*max_walking_dis/dis_points*dis_lon
#   r_max_lon <- max_lon + 2*max_walking_dis/dis_points*dis_lon
#   
#   lat_seq <- seq(r_min_lat, r_max_lat, by=dis_lat)
#   lon_seq <- seq(r_min_lon, r_max_lon, by=dis_lon)
#   
#   points <- merge(lat_seq,lon_seq)
#   colnames(points) <- c("lat","lon")
#   
#   #keep only points which are at distance less than max distance from nearest station
#   wdcTemp <- unique(wdcMerged[,c("lat","lon")])
#   lat2 <- wdcTemp[,"lat"]
#   lon2 <- wdcTemp[,"lon"]
#   
#   keep_index_func <- function(t) {
#     dis_v <- latlondistance(t[1], t[2], lat2, lon2)    
#     if(min(dis_v)<max_walking_dis ) return(1)
#     return(0)
#   }
#   points$keep_index <- apply(points[,c("lat","lon")],1,FUN=keep_index_func)
#   points <- subset(points, keep_index==1)
#   points$keep_index <- NULL
#   
#   #assign points to tracts
# #   wdcCenTrParsed <- readtractboundaryfile()
# #   wdcCenTrDensity <- readtractdensityfile()
#   
#   #find distance from each points row to each point in census and assign one that 
#   #is closest
#   load(paste0(csv_dir,"/ffdb/Paris/st_cluster"))    
#   points$tract <- assign_points_cluster_tract(points$lat,points$lon,st_cluster)
#   
#   #merge points with density data
#   #points <- merge(points,wdcCenTrDensity,by="tract")
#   
#   #scale density according to the area around point. 
#   #it is a rectangle of dis_points*dis_points*10^6 as distance is in kms and density in /m2
# #   points$density <- points$density *dis_points*dis_points*10^6
# #   points$density <- points$density /60 #to get better fixed effects in numerical computation range.
#   points$density <- 1
#   #for each point keep on lat lon density
#   #points$tract <- NULL
#   points$Population <- NULL
#   points$ALAND10 <- NULL
# #   unlink(wdcCenTrParsed)
# #   unlink(wdcCenTrDensity)
#   return(points)
# }

