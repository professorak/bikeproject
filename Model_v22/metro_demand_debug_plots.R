#ran t_rb30.1_p6_ls_stepf_full_debug.R before, plotting here
source("GetDataGooglePlaces_20dis_aggmonth_averaged_pr6.R")
# source("GetDataGooglePlaces_20dis_aggmonth_bytw_averaged_pr6.R")
# tw_in <- 1
# wdcMerged <- subset(wdcMerged, tw==tw_in)
points_metro <- subset(points, type==2)
wdcCenTrParsed <- readtractboundaryfile_presentation()
#
#Plot metro stations and velib stations with bubble size as their demand (both normalized so highest is 1)
#use different colors
#Then another plot with demand coming from metro stations
#and another with demand from intercept
rad_vec <- c(rep(0.1,nrow(wdcCenTrParsed)), sqrt(wdcMerged$out_dem_mean/max(wdcMerged$out_dem_mean)), 
             sqrt(points_metro$weight/max(points_metro$weight)))
rad_vec <- c(rep(0.1,nrow(wdcCenTrParsed)), sqrt(demand_division_df[,"density metro"]/max(demand_division_df[,"density metro"])), 
             sqrt(points_metro$weight/max(points_metro$weight)))
rad_vec <- c(rep(0.1,nrow(wdcCenTrParsed)), sqrt(demand_division_df[,"density intercept"]/max(demand_division_df[,"density metro"])), 
             sqrt(points_metro$weight/max(points_metro$weight)))


lon_vec <- c(wdcCenTrParsed$lon, wdcMerged$lon, points_metro$lon)
lat_vec <- c(wdcCenTrParsed$lat, wdcMerged$lat, points_metro$lat)

fg_vec <- c(rep("grey",nrow(wdcCenTrParsed)),rep("black",nrow(wdcMerged)),rep("red",nrow(points_metro)))
bg_vec <- fg_vec

symbols(lon_vec, lat_vec, circles=rad_vec,
        inches=0.06, fg=fg_vec, bg=bg_vec,
        cex.axis=0.7,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
#textxy(wdcMerged$lon,wdcMerged$lat,labs=c(1:nrow(wdcMerged))) 


summary(lm(wdcMerged$out_dem_mean ~ (wdcMerged$metro_den_on_1>0)))
cor(wdcMerged$out_dem_mean,wdcMerged$metro_den_on_1)
cor(wdcMerged$out_dem_mean,wdcMerged$metro_den_on_1>0)
cor(wdcMerged$out_dem_mean,(wdcMerged$metro_den_on_1+wdcMerged$metro_den_on_2)>0)
cor(wdcMerged$out_dem_mean,(wdcMerged$metro_den_on_1+wdcMerged$metro_den_on_2+wdcMerged$metro_den_on_3)>0)



#stid_list1 <- c(9004,9005,18041,9018,18042,18043,9027,18004,18114,18003,18001,18113,18111,18017,18016,18015,18002,18006,18005)
stid_list1 <- c(15118 ,15002, 9106, 9018, 8017, 4018, 10016, 10007, 8039, 8013, 8030, 8010)

View(cbind(wdcMerged$out_dem_mean[which(wdcMerged$station_id %in% stid_list1)],
wdcMerged$station_id[which(wdcMerged$station_id %in% stid_list1)]))


dem_nearby_metro_stations <- cbind(wdcMerged$out_dem_mean[which(wdcMerged$metro_den_on_1>0)],
                                   wdcMerged$metro_den_on_1[which(wdcMerged$metro_den_on_1>0)],
                                   wdcMerged$station_id[which(wdcMerged$metro_den_on_1>0)])

dem_nearby_metro_stations <- dem_nearby_metro_stations[order(dem_nearby_metro_stations[,2]),]

wdcMerged$metro_den_on_1_2 <- wdcMerged$metro_den_on_1 + wdcMerged$metro_den_on_2

dem_nearby_metro_stations <- cbind(wdcMerged$out_dem_mean[which(wdcMerged$metro_den_on_1_2>0)],
                                   wdcMerged$metro_den_on_1_2[which(wdcMerged$metro_den_on_1_2>0)],
                                   wdcMerged$station_id[which(wdcMerged$metro_den_on_1_2>0)])

dem_nearby_metro_stations <- dem_nearby_metro_stations[order(dem_nearby_metro_stations[,2]),]







