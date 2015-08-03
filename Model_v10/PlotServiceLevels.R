source("data_estimation_2.6_weather_saved.R")
points_save <- points
wdcMerged_store <- wdcMerged
user_serv_lvl_store <- user_serv_lvl
user_serv_lvl <- user_serv_lvl[order(user_serv_lvl$tw_group, user_serv_lvl$station_id_index),]

user_serv_lvl_6 <- user_serv_lvl[which(user_serv_lvl$tw_group %like% "6_4"),]
user_serv_lvl_7 <- user_serv_lvl[which(user_serv_lvl$tw_group %like% "7_4"),]
user_serv_lvl_8 <- user_serv_lvl[which(user_serv_lvl$tw_group %like% "8_4"),]

ord <- order(user_serv_lvl_6$serv_lvl)

plot(user_serv_lvl_6$serv_lvl[ord], cex=0.5, col=rainbow(3)[1], xlab="Not station id",
     ylab="Avg Service level time-window 4-8pm")
lines(user_serv_lvl_7$serv_lvl[ord], cex=0.5, col=rainbow(3)[2], type="p")
lines(user_serv_lvl_8$serv_lvl[ord], cex=0.5, col=rainbow(3)[3], type="p")
legend("topright",c("Jun","Jul","Aug"), col= rainbow(3), lty=1, cex=0.5)
