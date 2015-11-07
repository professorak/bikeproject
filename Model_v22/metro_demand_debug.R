#note down closest station to metro, its demand, its distance
#note down 2 closest station to metro, its demand, its distance to 2nd station
#metro station
points_metro <- subset(points, type==2)
points_metro <- points_metro[order(points_metro$weight, decreasing = T),]
metro_demand_df <- c()
for(i in c(1:nrow(points_metro))) {
  #also not lat,lon, tract
  dis_v <- latlondistance(points_metro$lat[i],points_metro$lon[i],
                          wdcMerged$lat, wdcMerged$lon)
  top_3_st <- order(dis_v)[1:3]
  metro_demand_df <- rbind(metro_demand_df,
  c(points_metro[i,c("tract","lat","lon","weight")],
  top_3_st[1],wdcMerged$out_dem_mean[top_3_st[1]],dis_v[top_3_st[1]],
  top_3_st[2],wdcMerged$out_dem_mean[top_3_st[2]],dis_v[top_3_st[2]],
  top_3_st[3],wdcMerged$out_dem_mean[top_3_st[3]],dis_v[top_3_st[3]]))  
}
colnames(metro_demand_df) <- c("tract","lat","lon","weight",
                               "st_1","dem_1","dis_1",
                               "st_2","dem_2","dis_2",
                               "st_3","dem_3","dis_3") 
metro_demand_df <- as.data.frame(metro_demand_df)

cor(as.numeric(metro_demand_df$weight), as.numeric(metro_demand_df$dem_1))
cor(as.numeric(metro_demand_df$weight[which(metro_demand_df$tract<=10)])[-c(1:10)], 
    as.numeric(metro_demand_df$dem_2[which(metro_demand_df$tract<=10)])[-c(1:10)])

plot(as.numeric(metro_demand_df$dem_2[which(metro_demand_df$tract<=10)]))
par(new=T)
plot(as.numeric(metro_demand_df$weight[which(metro_demand_df$tract<=10)]))


