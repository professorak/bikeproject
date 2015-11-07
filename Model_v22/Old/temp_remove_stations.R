wdcCenTrParsed <- readtractboundaryfile_presentation()

rad_vec <- c(rep(0.1,nrow(wdcCenTrParsed)), 
             (wdcMerged$out_dem_mean/max(wdcMerged$out_dem_mean)))

lon_vec <- c(wdcCenTrParsed$lon, wdcMerged$lon)
lat_vec <- c(wdcCenTrParsed$lat, wdcMerged$lat)

fg_vec <- c(rep("grey",nrow(wdcCenTrParsed)),rep("black",nrow(wdcMerged)))

symbols(lon_vec, lat_vec, circles=rad_vec,
        inches=0.06, fg=fg_vec, bg=bg_vec,
        cex.axis=0.7,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")

outliers_idx <- which(wdcMerged$tract==12 | wdcMerged$tract==16)

textxy(wdcMerged$lon[outliers_idx],wdcMerged$lat[outliers_idx],
       labs=outliers_idx, cex=0.6) 

outliers_idx <- c(723,724,725,959,726,727,442,455,460,456,457,458,459,960,961)
fg_vec[outliers_idx+nrow(wdcCenTrParsed)] <- "red"


save_dir <- paste0(csv_dir,"/ffdb/Paris/paris_stations.RData")
load(save_dir)
paris_stations <- paris_stations[-outliers_idx]
save(paris_stations, file=save_dir)


list <- which(wdcMerged$station_id %in% paris_stations)
