#get total demand at stations in nearby 100,200,..,600mts at a metro station

###### Load saved file if exists
file <- "total_station_demand_linearreg"
filepath <- paste0(csv_dir,"/Paris/",file,".RData")
if(!file.exists(filepath)) {
  wdcMerged <- get_total_station_demand()  
  save(wdcMerged,file=filepath)
}
###
load(filepath)  

wdcMerged$metro_nearby <- (wdcMerged$metro_den_on_1>0)
wdcMerged$metro_nearby_1_2 <- ((wdcMerged$metro_den_on_1+wdcMerged$metro_den_on_2)>0)

points <- generate_integration_points()

#keep only metro stations
points <- subset(points, type==2)
points <- subset(points, tract %in% c(1:10))

#generate total demand in 100mts etc. and number of stations
points_nearby_dem <- points
radius_vec <- c(1:6)/10 #100mts to 600mts
for(j in 1:length(radius_vec)) {
  t_j <- generate_points_demand(radius_vec[j])
  colnames(t_j) <- paste(c("demand","no_st"),rep(radius_vec[j],2),sep="_")
  points_nearby_dem <- cbind(points_nearby_dem,t_j)
}

#order points_nearby_dem by tract
points_nearby_dem <- points_nearby_dem[order(points_nearby_dem$tract),]

# plot.new()
# par(mar=c(5,4,4,5)+.1)
# plot(points_nearby_dem$density, type="l",col="blue",ylab="metro demand")
# par(new=TRUE)
# yrange=range(c(points_nearby_dem$demand_0.1,points_nearby_dem$demand_0.2,
#                points_nearby_dem$demand_0.3,points_nearby_dem$demand_0.4,
#                points_nearby_dem$demand_0.5,points_nearby_dem$demand_0.6))
# plot(points_nearby_dem$demand_0.1,type="l",col="red",xaxt="n",yaxt="n",xlab="",ylab="",
#     ylim=yrange)
# lines(points_nearby_dem$demand_0.2,type="l",col="orange",xaxt="n",yaxt="n",xlab="",ylab="",
#       ylim=yrange)
# lines(points_nearby_dem$demand_0.3,type="l",col="green",xaxt="n",yaxt="n",xlab="",ylab="",
#       ylim=yrange)
# lines(points_nearby_dem$demand_0.4,type="l",col="darkorchid",xaxt="n",yaxt="n",xlab="",ylab="",
#       ylim=yrange)
# lines(points_nearby_dem$demand_0.5,type="l",col="deeppink",xaxt="n",yaxt="n",xlab="",ylab="",
#       ylim=yrange)
# lines(points_nearby_dem$demand_0.6,type="l",col="firebrick",xaxt="n",yaxt="n",xlab="",ylab="",
#       ylim=yrange)
# axis(4)
# mtext("total station demand",side=4,line=3)
# legend("topleft",col=c("blue","red"),lty=1,legend=c("metro demand","station demand"), cex=0.5)


#plot in order of tracts
plot.new()
par(mar=c(5,4,4,5)+.1)
plot(points_nearby_dem$density, type="l",col="blue",ylab="metro demand")
par(new=TRUE)
# yrange=range(c(points_nearby_dem$demand_0.1,points_nearby_dem$demand_0.2,
#                points_nearby_dem$demand_0.3,points_nearby_dem$demand_0.4,
#                points_nearby_dem$demand_0.5,points_nearby_dem$demand_0.6))
# plot(points_nearby_dem$demand_0.1,type="l",col="red",xaxt="n",yaxt="n",xlab="",ylab="",
#      ylim=yrange)
# lines(points_nearby_dem$demand_0.2,type="l",col="orange",xaxt="n",yaxt="n",xlab="",ylab="",
#       ylim=yrange)
#plot(points_nearby_dem$demand_0.3,type="l",col="green",xaxt="n",yaxt="n",xlab="",ylab="")
# lines(points_nearby_dem$demand_0.4,type="l",col="darkorchid",xaxt="n",yaxt="n",xlab="",ylab="",
#       ylim=yrange)
# lines(points_nearby_dem$demand_0.5,type="l",col="deeppink",xaxt="n",yaxt="n",xlab="",ylab="",
#       ylim=yrange)
plot(points_nearby_dem$demand_0.6,type="l",col="firebrick",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("total station demand",side=4,line=3)
legend("topleft",col=c("blue","red"),lty=1,legend=c("metro","station"), cex=0.3)


#plot in order of metro density
order_metro_den <- order(points_nearby_dem$density)
plot.new()
par(mar=c(5,4,4,5)+.1)
plot(points_nearby_dem$density[order_metro_den], type="l",col="blue",ylab="metro demand")
par(new=TRUE)
# yrange=range(c(points_nearby_dem$demand_0.1,points_nearby_dem$demand_0.2,
#                points_nearby_dem$demand_0.3,points_nearby_dem$demand_0.4,
#                points_nearby_dem$demand_0.5,points_nearby_dem$demand_0.6))
# plot(points_nearby_dem$demand_0.1[order_metro_den],type="l",col="red",xaxt="n",yaxt="n",xlab="",ylab="",
#      ylim=yrange)
# lines(points_nearby_dem$demand_0.2[order_metro_den],type="l",col="orange",xaxt="n",yaxt="n",xlab="",ylab="",
#       ylim=yrange)
plot(points_nearby_dem$demand_0.3[order_metro_den],type="l",col="green",xaxt="n",yaxt="n",xlab="",ylab="")
# lines(points_nearby_dem$demand_0.4[order_metro_den],type="l",col="darkorchid",xaxt="n",yaxt="n",xlab="",ylab="",
#       ylim=yrange)
# lines(points_nearby_dem$demand_0.5[order_metro_den],type="l",col="deeppink",xaxt="n",yaxt="n",xlab="",ylab="",
#       ylim=yrange)
# plot(points_nearby_dem$demand_0.6[order_metro_den],type="l",col="firebrick",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("total station demand",side=4,line=3)
legend("topleft",col=c("blue","red"),lty=1,legend=c("metro","station"), cex=0.3)




#correlation
cor(points_nearby_dem$density,points_nearby_dem$demand_0.3)
cor(log(points_nearby_dem$density),points_nearby_dem$demand_0.3)
cor(sqrt(points_nearby_dem$density),points_nearby_dem$demand_0.3)

cor(points_nearby_dem$density[order_metro_den][1:80],points_nearby_dem$demand_0.3[order_metro_den][1:80])
cor(points_nearby_dem$density[order_metro_den][1:80],points_nearby_dem$demand_0.6[order_metro_den][1:80])

cor(points_nearby_dem$density[order_metro_den][1:80],
    (points_nearby_dem$demand_0.6-points_nearby_dem$demand_0.3)[order_metro_den][1:80])



#plot catchment area
colors <- rainbow(20)
colors_count <- 1
order_outdem <- order(wdcMerged$tract,wdcMerged$out_dem_mean/wdcMerged$catchment_area_step2)
yrange <- c(-2,10)
ycount <- 9
leg_text <- c()

plot.new()
par(mar=c(5,4,4,5)+.1)
plot(normalize_vec(wdcMerged$out_dem_mean[order_outdem])+ycount, type="l",col=colors[colors_count],ylab="station demand",ylim=yrange)
ycount <- ycount - 1
colors_count <- colors_count + 1
leg_text <- c(leg_text, "out_dem_mean")

par(new=TRUE)
plot(normalize_vec(wdcMerged$catchment_area_step2[order_outdem])+ycount,
     type="l",col=colors[colors_count],xaxt="n",yaxt="n",xlab="",ylab="",ylim=yrange)
ycount <- ycount - 1
colors_count <- colors_count + 1
leg_text <- c(leg_text, "catchment_area_step2")

par(new=TRUE)
plot(normalize_vec((wdcMerged$out_dem_mean/wdcMerged$catchment_area_step2)[order_outdem])+ycount,
     type="l",col=colors[colors_count],xaxt="n",yaxt="n",xlab="",ylab="",ylim=yrange)
ycount <- ycount - 1
colors_count <- colors_count + 1
leg_text <- c(leg_text, "dem/catchment_area_step2")

par(new=TRUE)
plot(normalize_vec(wdcMerged$metro_den_on_1[order_outdem])+ycount
                   ,type="l",col=colors[colors_count],xaxt="n",yaxt="n",xlab="",ylab="",ylim=yrange)
ycount <- ycount - 1
colors_count <- colors_count + 1
leg_text <- c(leg_text, "metro_den_on_1")

par(new=TRUE)
plot(normalize_vec(wdcMerged$census_density[order_outdem])+ycount,type="l",col=colors[colors_count],xaxt="n",yaxt="n",xlab="",ylab="",ylim=yrange)
ycount <- ycount - 1
colors_count <- colors_count + 1
leg_text <- c(leg_text, "census_density")

par(new=TRUE)
plot(normalize_vec(wdcMerged$googleplaces_den_1[order_outdem])+ycount,type="l",col=colors[colors_count],xaxt="n",yaxt="n",xlab="",ylab="",ylim=yrange)
ycount <- ycount - 1
colors_count <- colors_count + 1
leg_text <- c(leg_text, "googleplaces_den_1")

par(new=TRUE)
plot(normalize_vec(wdcMerged$googleplaces_food_den_1[order_outdem])+ycount,
     type="l",col=colors[colors_count],xaxt="n",yaxt="n",xlab="",ylab="",ylim=yrange)
ycount <- ycount - 1
colors_count <- colors_count + 1
leg_text <- c(leg_text, "googleplaces_food_den_1")

par(new=TRUE)
plot(normalize_vec(wdcMerged$googleplaces_grocery_den_1[order_outdem])+ycount,
     type="l",col=colors[colors_count],xaxt="n",yaxt="n",xlab="",ylab="",ylim=yrange)
ycount <- ycount - 1
colors_count <- colors_count + 1
leg_text <- c(leg_text, "googleplaces_grocery_den_1")

par(new=TRUE)
plot(normalize_vec(wdcMerged$googleplaces_government_den_1[order_outdem])+ycount,
     type="l",col=colors[colors_count],xaxt="n",yaxt="n",xlab="",ylab="",ylim=yrange)
ycount <- ycount - 1
colors_count <- colors_count + 1
leg_text <- c(leg_text, "googleplaces_government_den_1")

legend("bottomright",
       col=colors[1:length(leg_text)],lty=1,legend=leg_text, cex=0.5)


cor(wdcMerged$out_dem_mean,wdcMerged$catchment_area_step2)
summary(lm(wdcMerged$out_dem_mean ~ wdcMerged$catchment_area_step2))
summary(lm(wdcMerged$out_dem_mean ~ factor(wdcMerged$tract)))
summary(lm(wdcMerged$out_dem_mean ~ factor(wdcMerged$tract) + wdcMerged$catchment_area_step2))

#outliers
wdcMerged[which(wdcMerged$out_dem_mean>=0.20),
          c("lat","lon","station_id","station_id_index","tract","out_dem_mean","catchment_area_step2","metro_den_on_1","metro_den_on_2",
            "googleplaces_den_1","googleplaces_food_den_1","googleplaces_grocery_den_1")]

wdcMerged[which(wdcMerged$station_id==10029),
          c("lat","lon","station_id","station_id_index","tract","out_dem_mean","catchment_area_step2","metro_den_on_1","metro_den_on_2",
            "googleplaces_den_1","googleplaces_food_den_1","googleplaces_grocery_den_1")]



#location of gare du nord station
points_metro <- points[which(points$type==2),]
  #48.87997,2.354703, also has the highest demand.

cor(wdcMerged$out_dem_mean,wdcMerged$metro_den_on_2)
cor(wdcMerged$out_dem_mean,log(wdcMerged$metro_den_on_2+1))

wdcMerged$metro_den_on_1_2 <- wdcMerged$metro_den_on_1*4.83542598489239e-05 + wdcMerged$metro_den_on_2*0.000170904722330068
cor(wdcMerged$out_dem_mean,wdcMerged$metro_den_on_1_2)
cor(wdcMerged$out_dem_mean,wdcMerged$metro_den_on_1_2!=0)
summary(lm(wdcMerged$out_dem_mean ~ wdcMerged$metro_den_on_1_2))
summary(lm(wdcMerged$out_dem_mean ~ wdcMerged$tract + wdcMerged$metro_den_on_1_2))
summary(lm(wdcMerged$out_dem_mean ~ wdcMerged$tract + wdcMerged$metro_den_on_1))
summary(lm(wdcMerged$out_dem_mean ~ wdcMerged$tract + wdcMerged$metro_den_on_1+ wdcMerged$metro_den_on_2))
summary(lm(wdcMerged$out_dem_mean ~ wdcMerged$tract + wdcMerged$metro_den_on_1+ wdcMerged$metro_den_on_1_2))
summary(lm(wdcMerged$out_dem_mean ~ wdcMerged$tract + wdcMerged$metro_den_on_1+ wdcMerged$metro_den_on_1_2 + wdcMerged$metro_nearby + 
             wdcMerged$metro_nearby_1_2))
summary(lm(wdcMerged$out_dem_mean ~ wdcMerged$census_density + wdcMerged$metro_den_on_1+ wdcMerged$metro_den_on_1_2 + wdcMerged$metro_nearby + 
             wdcMerged$metro_nearby_1_2))
summary(lm(wdcMerged$out_dem_mean ~ wdcMerged$census_density + log(wdcMerged$metro_den_on_1*10000+1)+ log(wdcMerged$metro_den_on_1_2*100000+1) + wdcMerged$metro_nearby + 
             wdcMerged$metro_nearby_1_2))
summary(lm(wdcMerged$out_dem_mean ~ wdcMerged$census_density +  wdcMerged$metro_nearby + 
             wdcMerged$metro_nearby_1_2))

#plot outdemmean agains metro_on_1+metro_on_2 variable to deduce shape of metro variable.
order_metro_den <- order(wdcMerged$tract,wdcMerged$out_dem_mean)
plot.new()
par(mar=c(5,4,4,5)+.1)
plot(wdcMerged$out_dem_mean[order_metro_den], type="l",col="blue",ylab="station demand")
par(new=TRUE)
plot((wdcMerged$metro_den_on_1_2)[order_metro_den],type="l",col="red",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("nearby metro demand",side=4,line=3)
legend("topleft",col=c("blue","red"),lty=1,legend=c("station","metro"), cex=1.0)




