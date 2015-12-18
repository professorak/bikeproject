dim(wdcMerged)
#stations that are too close to each other
dis_next <- c()
for(i in c(1:nrow(wdcMerged))) {
  dis_next <- c(dis_next, min(latlondistance(wdcMerged$lat[i],wdcMerged$lon[i],
                                   wdcMerged$lat[-i],wdcMerged$lon[-i])))
  
}

dis_next[order(dis_next)[1:10]]
length(which(dis_next<=0.100))
length(which(dis_next<=0.050))
length(which(dis_next<=0.020))

points_metro <- subset(points, type==2)
wdcCenTrParsed <- readtractboundaryfile_presentation()

st_idx <- which(dis_next<=0.050)

rad_vec <- c(rep(0.1,nrow(wdcCenTrParsed)), 
             (wdcMerged$out_dem_mean[st_idx]/max(wdcMerged$out_dem_mean)), 
             ((points_metro$weight/max(points_metro$weight))^2)*0.3)


lon_vec <- c(wdcCenTrParsed$lon, wdcMerged$lon[st_idx], points_metro$lon)
lat_vec <- c(wdcCenTrParsed$lat, wdcMerged$lat[st_idx], points_metro$lat)

fg_vec <- c(rep("grey",nrow(wdcCenTrParsed)),rep("black",length(st_idx)),rep("red",nrow(points_metro)))
#fg_vec[nrow(wdcCenTrParsed)+outliers_idx] <- "blue"
bg_vec <- fg_vec

symbols(lon_vec, lat_vec, circles=rad_vec,
        inches=0.06, fg=fg_vec, bg=bg_vec,
        cex.axis=0.7,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")

textxy(wdcMerged$lon[st_idx],wdcMerged$lat[st_idx],
       labs=st_idx, cex=0.6) 

