require("calibrate")

#plot demand for stations and metro moment variables to decude functional form.
wdcMerged$metro_den_on_1_2 <- wdcMerged$metro_den_on_1 + wdcMerged$metro_den_on_2
wdcMerged$out_dem_mean

plot(wdcMerged$out_dem_mean, type="l")
plot(wdcMerged$metro_den_on_1_2, type="l")
cor(wdcMerged$metro_den_on_1_2, wdcMerged$out_dem_mean)

ord <- order(wdcMerged$tract, wdcMerged$metro_den_on_1_2)
plot(wdcMerged$out_dem_mean[ord], type="l")
plot(wdcMerged$metro_den_on_1_2[ord], type="l")
plot(log(wdcMerged$metro_den_on_1_2+1e-6)[ord], type="l")

cor(log(wdcMerged$metro_den_on_1_2+1e-6), wdcMerged$out_dem_mean)
cor((wdcMerged$metro_den_on_1_2>0), wdcMerged$out_dem_mean)


#looking at xi and demand values of stations nearby gare du nord
#changing to gard de l'est (48.8765863  2.3591431)
#idx <- which(points$type==2 & points$weight==1) #gare du nord
#idx <- which(points$lat==48.8809807 & points$lon==2.3553404) #gare du nord
idx <- which(points$lat==48.8765863 & points$lon==2.3591431) #gare de l'est
idx_locst <- splitchar(points$local_stations[idx])
xi[idx_locst]
wdcMerged$out_dem_sum[idx_locst]
latlondistance(wdcMerged$lat[idx_locst], wdcMerged$lon[idx_locst],
               points$lat[idx],points$lon[idx])

symbols(c(wdcMerged$lon[idx_locst],points$lon[idx]), c(wdcMerged$lat[idx_locst],points$lat[idx]),
        circles=c(sqrt(wdcMerged$out_dem_sum[idx_locst]),1), inches=0.5)
textxy(c(wdcMerged$lon[idx_locst],points$lon[idx]), c(wdcMerged$lat[idx_locst],points$lat[idx]),
       labs=c(c(1:length(idx_locst)),"S")) 

# #possible exits for gare du nord, based on google maps data
# nord_exits_df <- as.data.frame(rbind(c(48.880680, 2.356703),
#                       c(48.879918, 2.356363),
#                       c(48.880941, 2.354059)))
# colnames(nord_exits_df)  <- c("lat","lon")
# 
# symbols(c(wdcMerged$lon[idx_locst],points$lon[idx],nord_exits_df$lon), 
#         c(wdcMerged$lat[idx_locst],points$lat[idx],nord_exits_df$lat),
#         circles=c(sqrt(wdcMerged$out_dem_sum[idx_locst]),1,rep(1,nrow(nord_exits_df))), inches=0.5,
#         xlab="lon",ylab="lat")
# textxy(c(wdcMerged$lon[idx_locst],points$lon[idx],nord_exits_df$lon), 
#        c(wdcMerged$lat[idx_locst],points$lat[idx],nord_exits_df$lat),
#        labs=c(c(1:length(idx_locst)),"S",rep("E",nrow(nord_exits_df)))) 



