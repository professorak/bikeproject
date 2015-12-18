require("calibrate")
outliers_idx <- which(xi^2>0.3)


#plots
deltain_stin_cf <- deltain_stin - xi + mean(xi)
lambda_cf <- eval_lambda_full(deltain_stin_cf, theta1, wdcMerged, points)

points_metro <- subset(points, type==2)
wdcCenTrParsed <- readtractboundaryfile_presentation()

rad_vec <- c(rep(0.1,nrow(wdcCenTrParsed)), 
             (lambda_cf/max(lambda_cf)), 
             ((points_metro$weight/max(points_metro$weight))^2)*0.3)
rad_vec <- c(rep(0.1,nrow(wdcCenTrParsed)), 
             (wdcMerged$out_dem_mean/max(wdcMerged$out_dem_mean)), 
             ((points_metro$weight/max(points_metro$weight))^2)*0.3)


lon_vec <- c(wdcCenTrParsed$lon, wdcMerged$lon, points_metro$lon)
lat_vec <- c(wdcCenTrParsed$lat, wdcMerged$lat, points_metro$lat)

fg_vec <- c(rep("grey",nrow(wdcCenTrParsed)),rep("black",nrow(wdcMerged)),rep("red",nrow(points_metro)))
fg_vec[nrow(wdcCenTrParsed)+outliers_idx] <- "blue"
bg_vec <- fg_vec

symbols(lon_vec, lat_vec, circles=rad_vec,
        inches=0.06, fg=fg_vec, bg=bg_vec,
        cex.axis=0.7,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")

textxy(wdcMerged$lon[outliers_idx],wdcMerged$lat[outliers_idx],
       labs=outliers_idx, cex=0.6) 



##################################################################
##################################################################
#major change between theta1 and theta1_case2
outliers_idx_1 <- which((xi_case2^2-xi^2)>0.05)
outliers_idx_2 <- which((xi_case2^2-xi^2) < -0.05)

deltain_stin_cf <- deltain_stin - xi + mean(xi)
lambda_cf <- eval_lambda_full(deltain_stin_cf, theta1, wdcMerged, points)

deltain_stin_cf_case2 <- deltain_stin_case2 - xi_case2 + mean(xi_case2)
lambda_cf_case2 <- eval_lambda_full(deltain_stin_cf_case2, theta1_case2, wdcMerged, points)

points_metro <- subset(points, type==2)
wdcCenTrParsed <- readtractboundaryfile_presentation()

rad_vec <- c(rep(0.1,nrow(wdcCenTrParsed)), 
             (lambda_cf/max(lambda_cf)), 
             ((points_metro$weight/max(points_metro$weight))^2)*0.3)
rad_vec <- c(rep(0.1,nrow(wdcCenTrParsed)), 
             (lambda_cf_case2/max(lambda_cf_case2)), 
             ((points_metro$weight/max(points_metro$weight))^2)*0.3)


lon_vec <- c(wdcCenTrParsed$lon, wdcMerged$lon, points_metro$lon)
lat_vec <- c(wdcCenTrParsed$lat, wdcMerged$lat, points_metro$lat)

fg_vec <- c(rep("grey",nrow(wdcCenTrParsed)),rep("black",nrow(wdcMerged)),rep("red",nrow(points_metro)))
fg_vec[nrow(wdcCenTrParsed)+outliers_idx_1] <- "blue"
fg_vec[nrow(wdcCenTrParsed)+outliers_idx_2] <- "green"
bg_vec <- fg_vec

symbols(lon_vec, lat_vec, circles=rad_vec,
        inches=0.06, fg=fg_vec, bg=bg_vec,
        cex.axis=0.7,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")

textxy(wdcMerged$lon[c(outliers_idx_1,outliers_idx_2)],
       wdcMerged$lat[c(outliers_idx_1,outliers_idx_2)],
       labs=c(outliers_idx_1,outliers_idx_2), cex=0.6) 



