tract <- c(1:20)

st_list <<- select_stations_in_tract(c(tract))
length(st_list)
#temporary - load st_6_7_8 file, should be able to generate this realtime, 
save_dir <- paste0(csv_dir,"/ffdb/Paris/st_5_8")
load(save_dir)
st_list <<- intersect(st_5_8,st_list)

#get station info: (lat,lon) from saved files
load(file="wdcMerged_lin_weather_dis1to20_aggmonths_averaged_pr6.RData")
st_data <- wdcMerged[which(!duplicated(wdcMerged$station_id)),c("lat","lon","station_id")]

source("tract_assignment_functions.R")

#set parameters
discrete_point_size <- 0.025 #distance between consecutive points of ray and also the boundary of tracts
intersection_gap <- 0.025 # 20mts
skip_upon_intersection <- 0.200 #200mts
wdcCenTrParsed <- readtractboundaryfile_finegrain(discrete_point_size)

tracts_st_data <- list()
for(i in 1:nrow(st_data)) {
  tracts_st_data <- c(tracts_st_data,
                      list(classify_to_tract(c(st_data$lat[i],st_data$lon[i]),wdcCenTrParsed,discrete_point_size,
                                        intersection_gap, skip_upon_intersection)))
}
