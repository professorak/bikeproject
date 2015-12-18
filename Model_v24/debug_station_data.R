stid <- 10151
load(file=paste0("station_data_",stid,".RData"))

summary(station_data$out_dem)
summary(station_data$out_dem[which(station_data$stocked_out==F)])

#see how does demand look like by taking various stockout definition (>0), (>5)

#plot #bikes against out_dem.
plot(station_data$out_dem[1:1000])
par(new=T)
plot(station_data$bikes[1:1000])

#average demand in different time windows using (>0) stockout definition
station_data_stocked <- subset(station_data, stocked_out==F)
station_data_twdem <- binned_sum(station_data_stocked$out_dem, bin=station_data_stocked$tw+1)
station_data_twdem <- station_data_twdem[,"sum"]/station_data_twdem[,"count"]
station_data_twdem_df <- data.frame(tw=c(0:23), dem=station_data_twdem)
station_data_twdem_df$tw_group <- as.integer(station_data_twdem_df$tw/4)
by(station_data_twdem_df$dem, station_data_twdem_df$tw_group, FUN=mean)


length(which(station_data$bikes==0))/nrow(station_data)
length(which(station_data$bikes<=3))/nrow(station_data)
length(which(station_data$bikes<=5))/nrow(station_data)

sum(station_data$out_dem[which(station_data$bikes>0)])
sum(station_data$out_dem[which(station_data$bikes>5)])

sum(station_data$out_dem[which(station_data$bikes>0)])/length(which(station_data$bikes>0))
sum(station_data$out_dem[which(station_data$bikes>5)])/length(which(station_data$bikes>5))


load(file="station_data_stid_list2.RData")

stid <- 10016

station_data_stid <- subset(station_data, station_id==stid)
mean(station_data_stid$out_dem)
length(which(station_data_stid$bikes>0))/nrow(station_data_stid)
length(which(station_data_stid$bikes>3))/nrow(station_data_stid)
length(which(station_data_stid$bikes>5))/nrow(station_data_stid)
summary(station_data_stid$out_dem)
sum(station_data_stid$out_dem[which(station_data_stid$bikes>0)])/length(which(station_data_stid$bikes>0))
sum(station_data_stid$out_dem[which(station_data_stid$bikes>5)])/length(which(station_data_stid$bikes>5))
write.csv(station_data_stid,file="station_data_stid.csv")

#average demand in different time windows using (>0) stockout definition
station_data_stocked <- subset(station_data_stid, stocked_out==F)
station_data_twdem <- binned_sum(station_data_stocked$out_dem, bin=station_data_stocked$tw+1)
station_data_twdem <- station_data_twdem[,"sum"]/station_data_twdem[,"count"]
station_data_twdem_df <- data.frame(tw=c(0:23), dem=station_data_twdem)
station_data_twdem_df$tw_group <- as.integer(station_data_twdem_df$tw/4)
by(station_data_twdem_df$dem, station_data_twdem_df$tw_group, FUN=mean)
#computing service level
station_data_twservlvl <- binned_sum(station_data_stid$stocked_out==F, bin=station_data_stid$tw+1)
station_data_twservlvl <- station_data_twservlvl[,"sum"]/station_data_twservlvl[,"count"]
station_data_twservlvl_df <- data.frame(tw=c(0:23), serv_lvl=station_data_twservlvl)
station_data_twservlvl_df$tw_group <- as.integer(station_data_twservlvl_df$tw/4)
station_data_twdem_df$serv_lvl <- station_data_twservlvl_df$serv_lvl
#plot #bikes against out_dem.
plot(station_data_stid$out_dem,type="l")
par(new=T)
plot(station_data_stid$bikes,type="l")

a <- by(station_data_stid$out_dem, 
        station_data_stid[,c("tw", "day")], FUN=sum)

write.csv(station_data_stid, file="station_data_stid.csv")


#measure demand conditional in time windows 6,7,8.. on there being more than 10 bikes at 6AM.
#for each day, get number of bikes at 6
station_data_stid_6am <- subset(station_data_stid, n_time_int %% 86400 == 14400)
station_data_stid_6am <- subset(station_data_stid_6am, bikes>=10)

station_data_stid_highbikes <- subset(station_data_stid, day %in% station_data_stid_6am$day)

#average demand in different time windows using (>0) stockout definition
station_data_stocked <- subset(station_data_stid_highbikes, stocked_out==F)
station_data_twdem <- binned_sum(station_data_stocked$out_dem, bin=station_data_stocked$tw+1)
station_data_twdem <- station_data_twdem[,"sum"]/station_data_twdem[,"count"]
station_data_twdem_df <- data.frame(tw=c(0:23), dem=station_data_twdem)
station_data_twdem_df$tw_group <- as.integer(station_data_twdem_df$tw/4)
by(station_data_twdem_df$dem, station_data_twdem_df$tw_group, FUN=mean)
#computing service level
station_data_twservlvl <- binned_sum(station_data_stid_highbikes$stocked_out==F, bin=station_data_stid_highbikes$tw+1)
station_data_twservlvl <- station_data_twservlvl[,"sum"]/station_data_twservlvl[,"count"]
station_data_twservlvl_df <- data.frame(tw=c(0:23), serv_lvl=station_data_twservlvl)
station_data_twservlvl_df$tw_group <- as.integer(station_data_twservlvl_df$tw/4)
station_data_twdem_df$serv_lvl <- station_data_twservlvl_df$serv_lvl


#create a table of demand for above stations with 
# 1. average demand based on 0 stockout
# 2. average demand based on 5 stockout
# 3. seperate average demand in 1 hr tw based on 5 stockout and average of these


station_data

# 0. average demand based on 0 stockout and rellocation at 4
dem_prev <- by(station_data$out_dem[which(station_data$bikes>0 & station_data$out_dem<=4)],
            station_data$station_id[which(station_data$bikes>0 & station_data$out_dem<=4)],
            FUN=mean)
dem_prev <- data.frame(station_id=as.numeric(names(dem_prev)), demand=as.numeric(dem_prev))

# 1. average demand based on 0 stockout
dem_0 <- by(station_data$out_dem[which(station_data$bikes>0)],
            station_data$station_id[which(station_data$bikes>0)],FUN=mean)
dem_0 <- data.frame(station_id=as.numeric(names(dem_0)), demand=as.numeric(dem_0))
# 2. average demand based on 5 stockout
dem_5 <- by(station_data$out_dem[which(station_data$bikes>5)],
            station_data$station_id[which(station_data$bikes>5)],FUN=mean)
dem_5 <- data.frame(station_id=as.numeric(names(dem_5)), demand=as.numeric(dem_5))
# 3. seperate average demand in 1 hr tw based on 5 stockout and average of these
dem_tw_5 <- ave(station_data$out_dem[which(station_data$bikes>5)],
                station_data[which(station_data$bikes>5),c("station_id","tw")],
                FUN=mean)
dem_tw_5 <- data.frame(station_id=station_data$station_id[which(station_data$bikes>5)],
                       tw=station_data$tw[which(station_data$bikes>5)],
                       demand=dem_tw_5)
dem_tw_5 <- dem_tw_5[which(!duplicated(dem_tw_5[,c("station_id","tw")])),]
dem_tw_5 <- dem_tw_5[order(dem_tw_5$station_id, dem_tw_5$tw),]
#check it is full
if(nrow(dem_tw_5) != length(unique(dem_tw_5$station_id)) * length(unique(dem_tw_5$tw))) stop("")
dem_tw_5 <- as.data.frame(matrix(dem_tw_5$demand, length(unique(dem_tw_5$station_id)),length(unique(dem_tw_5$tw)), byrow = T))
colnames(dem_tw_5) <- paste("tw_",c(0:23), split="")
dem_tw_5$demand_avg_tw <- rowMeans(dem_tw_5)
# 3. seperate average demand in 4 hr tw based on 5 stockout and average of these
station_data$tw_group <- as.integer(station_data$tw/4)
dem_twgroup_5 <- ave(station_data$out_dem[which(station_data$bikes>5)],
                     station_data[which(station_data$bikes>5),c("station_id","tw_group")],
                     FUN=mean)
dem_twgroup_5 <- data.frame(station_id=station_data$station_id[which(station_data$bikes>5)],
                            tw_group=station_data$tw_group[which(station_data$bikes>5)],
                            demand=dem_twgroup_5)
dem_twgroup_5 <- dem_twgroup_5[which(!duplicated(dem_twgroup_5[,c("station_id","tw_group")])),]
dem_twgroup_5 <- dem_twgroup_5[order(dem_twgroup_5$station_id, dem_twgroup_5$tw_group),]
#check it is full
if(nrow(dem_twgroup_5) != length(unique(dem_twgroup_5$station_id)) * length(unique(dem_twgroup_5$tw_group))) stop("")
dem_twgroup_5 <- as.data.frame(matrix(dem_twgroup_5$demand, length(unique(dem_twgroup_5$station_id)),length(unique(dem_twgroup_5$tw_group)), byrow = T))
colnames(dem_twgroup_5) <- paste("twgroup_",c(0:5), split="")
dem_twgroup_5$demand_avg_tw <- rowMeans(dem_twgroup_5)



dem_0 <- dem_0[order(dem_0$station_id),]
dem_5 <- dem_5[order(dem_5$station_id),]
dem_all <- cbind(dem_0$station_id,dem_tw_5, dem_0=dem_0$demand, dem_5=dem_5$demand,
                 dem_prev=dem_prev$demand,
                 dem_twgroup_5_avg = dem_twgroup_5$demand_avg_tw)




