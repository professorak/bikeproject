source("constants_mw6.R")
source("data_estimation_2.6_weather_functions.R")
#v3.0 two stockstate observations per stations, if results are consistent, can experiment with other fixed effects etc.
point_range <- 6

max_walking_dis <- 2.0
max_points_distance <- 0.6 #600mts is the maximum a point can be far way from closest station.
dis_points <<- 0.050
no_bootstrap_runs <- 1
reg_runs_save_file <- "reg_runs_discoef_2.89_8states_mw6_bootstrap.txt"
#num_top_states_select <- 8
stocked_out_threshold <<- 5

#length of stockout state for each station. Will still have other stations like before to use the same framework like before
hyperlocal_range <<- 6

source("eval_func_4_GMM_new_temp_2.89.R")
source("data_prep_functions_ffdf_limstostate.R")
source("data_prep_functions_ffdf_limstostate_hyperlocalstostate_withreallocthresh.R")
source("data_prep_inte_points_limstostate.R")
source("eval_func_2_cpp_cntrt_map_2.86_2.R")
source("eval_func_3_cpp_new_2.86_2.R")
source("temp_read_google_places_data.R")

#intergration points generation data
min_lon <<- 2.249428
max_lon <<- 2.46408
min_lat <<- 48.81626
max_lat <<- 48.89964

#max_walking_dis = 0.1
tract <- c(1:20)
  
st_list <<- select_stations_in_tract(c(tract))

#temporary - load st_6_7_8 file, should be able to generate this realtime, 
save_dir <- paste0(csv_dir,"/ffdb/Paris/st_5_8")
load(save_dir)
st_list <<- intersect(st_5_8,st_list)

#stations which are in interior of paris
save_dir <- paste0(csv_dir,"/ffdb/Paris/paris_stations.RData")
load(save_dir)
st_list <<- intersect(paris_stations,st_list)

filelist_all <- c("Paris/ind_paris_5.csv","Paris/ind_paris_6.csv",
                  "Paris/ind_paris_7.csv","Paris/ind_paris_8.csv")
filelist_all_mons <- c(5:8)
#   filelist_all <- c("Paris/ind_paris_7.csv")
#   filelist_all_mons <- c(7)
filelist_all_weather <- c("weatherdata_2013-5.csv","weatherdata_2013-6.csv",
                          "weatherdata_2013-7.csv","weatherdata_2013-8.csv")

#filelist_servlvl <- c("Paris/ind_paris_4.csv","Paris/ind_paris_5.csv")
wdcMerged <- c()
for( i in 1:length(filelist_all_mons)) {
  filelist <- filelist_all[i]
  mon <- filelist_all_mons[i]
  save_file <- paste0(csv_dir,"/ffdb/wdcMerged_bootstrap_hyperlocalstostate_withreallocthresh_stkoutthresh_den",
                      point_range,dis_points,"/mw6_tract_",paste0(tract,collapse="_"),"_mon_",mon)  
  print(save_file)
  #load files
  wdcMerged_bootstrap <- reg_run_bootstrap_simple(NULL,save_file)
  dim(wdcMerged_bootstrap)
  wdcMerged <- rbind(wdcMerged,wdcMerged_bootstrap)
  rm(wdcMerged_bootstrap)
}  
wdcMerged$n_time_int <- NULL #its not correct and sometimes used incorrectly.

# round wdcMerged$n_time_int to closest half hour interval corresponding to weather times
wdcMerged$time_halfhour_int <- wdcMerged$tw + round(wdcMerged$seconds/1800)/2
idx <- ffwhich(wdcMerged,time_halfhour_int==24)
wdcMerged$time_halfhour_int[idx] <- as.ff(rep(23.5, length(idx)))
wdcMerged$day_time_halfhour_int_fac <- as.ff(as.factor(paste(wdcMerged$day[],
        wdcMerged$time_halfhour_int[])))
#
#combine with weather data
# weather_data <- get_weather_data(filelist_all_weather)
# nrow_wdcMerged <- nrow(wdcMerged)
# system.time({
#   wdcMerged <- merge(wdcMerged,as.ffdf(weather_data), by=c("day","time_halfhour_int"), all.x=T)    
# })
# if(nrow(wdcMerged) < nrow_wdcMerged) {stop("ROWS missing in weather data")}
#deweatherize data
weather_scale_df <- get_weather_scale(wdcMerged)
#merge with wdcMerged
nrow_wdcMerged <- nrow(wdcMerged)
system.time({
  wdcMerged <- merge(wdcMerged,as.ffdf(weather_scale_df), by=c("day_time_halfhour_int_fac"))    
})
if(nrow(wdcMerged) < nrow_wdcMerged) {stop("ROWS missing in weather scale data")}

#scale demands
wdcMerged$out_dem <- wdcMerged$out_dem/wdcMerged$outdem_weather_scale
wdcMerged$in_dem <- wdcMerged$in_dem/wdcMerged$indem_weather_scale
#remove weather_state
wdcMerged$weather_state <- NULL
wdcMerged$weather_state <- as.ff(as.factor(rep("0_0_0_0",nrow(wdcMerged))))
#deweatherize data complete


wdcMerged_save <- wdcMerged
#wdcMerged_save <- subset(wdcMerged_save, tw>=6 & tw <10)
#wdcMerged_save <- droplevels(wdcMerged_save)  

#v0_vec <- generate_v0() #necessary to generate outside as seed is otherwise set inside, bad
v0_vec <- c(0)
set.seed(2)
i=1

print(i)   
wdcMerged <- wdcMerged_save
if(i==1) {
  a = as.ff(1:nrow(wdcMerged))
} else {
  a = as.ff(sample.int(nrow(wdcMerged), size = nrow(wdcMerged), replace = TRUE))
  wdcMerged = wdcMerged[a,]  
  wdcMerged <- droplevels(wdcMerged)
}
#construct alternative measures to stocked out, <=5 bikes, <=5% bikes
wdcMerged$station_size <- wdcMerged$bikes + wdcMerged$spaces_available
wdcMerged$less_5_bikes <- (wdcMerged$bikes<=5)
wdcMerged$less_3_bikes <- (wdcMerged$bikes<=3)
wdcMerged$less_7_bikes <- (wdcMerged$bikes<=7)
wdcMerged$less_5perc_bikes <- wdcMerged$bikes/wdcMerged$station_size
wdcMerged$less_5perc_bikes <- (wdcMerged$less_5perc_bikes<=0.05)

wdcMerged <- subset(wdcMerged, month>=5 & month<=8)
wdcMerged <- transform(wdcMerged, month=1) # to combine all months
wdcMerged$week <- wdcMerged$month

#assigning tw. tw 0:3 are night time windows when metro is off. 
list_time_halfhour_int <- list(c(0.5, 1.0,  1.5),  c(2.0,  2.5),  c(3.0,  3.5),
  c(4.0,  4.5, 5.0), c(5.5, 6.0,  6.5),  c(7.0,  7.5),  c(8.0,  8.5),  c(9.0,  9.5), c(10.0, 10.5),
  c(11.0, 11.5), c(12.0, 12.5), c(13.0, 13.5), c(14.0, 14.5), c(15.0, 15.5),
  c(16.0, 16.50), c(17.0, 17.5), c(18.0, 18.5), c(19.0, 19.5), c(20.0, 20.5), c(21.0, 21.5), c(22.0, 22.5), c(23.0, 23.5, 0.0))
map_tw <- c(0:(length(list_time_halfhour_int)-1))
for(j in c(1:length(list_time_halfhour_int))) {
  idx <- ffwhich(wdcMerged, time_halfhour_int %in% list_time_halfhour_int[[j]])
  wdcMerged$tw[idx] <- as.ff(rep(map_tw[j], length(idx)))
} 

wdcMerged$tw_fac <-  as.ff(as.character(wdcMerged$tw))
wdcMerged_save_i <- wdcMerged
#compute demand at the level of st,week,tw,sto_state

#wdcMerged <- subset(wdcMerged, stocked_out==FALSE)
#since aggdata aggregates over stid_twg_locstate_fac, 
#including weather in stid_twg_locstate_fac
wdcMerged$stid_twg_locstate_fac <-  as.ff(as.factor(paste(wdcMerged$station_id_index[],
                                                          wdcMerged$week[],wdcMerged$tw[],wdcMerged$stocked_out[],
                                                          wdcMerged$weather_state[]))) 
wdcMerged$stid_twg_locstate_fac <- droplevels(wdcMerged$stid_twg_locstate_fac)

wdcMerged <- aggdata(wdcMerged)
wdcMerged <- droplevels(wdcMerged)
wdcMerged$sto_state_local <- tapply(wdcMerged$sto_state_local,c(1:nrow(wdcMerged)),FUN=convert_sto_state_to_zero)
wdcMerged$sto_state_local_fac <- as.factor(wdcMerged$sto_state_local)

wdcMerged <- droplevels(wdcMerged)

wdcMerged$out_dem_mean <- wdcMerged$out_dem_sum/wdcMerged$obs_weight

wdcCenTrParsed <- readtractboundaryfile()
station_data <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
station_data$tract  <- assign_points_tract(station_data$lat,station_data$lon,wdcCenTrParsed)
station_data <- station_data[,c("station_id_index","tract")]
wdcMerged <- merge(wdcMerged,station_data,by="station_id_index")

wdcMerged$tract_tw <- as.factor(paste0(wdcMerged$tract[],"_",wdcMerged$tw))
wdcMerged$station_id_index_fac <- as.factor(wdcMerged$station_id_index)
wdcMerged$week_fac <- as.factor(wdcMerged$week)
wdcMerged$tract_tw_fac <- as.factor(wdcMerged$tract_tw)

#wdcMerged <- subset(wdcMerged, tw>=1 & tw<5)
#wdcMerged <- droplevels(wdcMerged)
wdcMerged_dem <- wdcMerged

#construct the service level at station level
wdcMerged <- wdcMerged_save_i
#wdcMerged <- subset(wdcMerged, tw>=1 & tw<5)    
#wdcMerged <- droplevels(wdcMerged)
#generate group by factor
wdcMerged$groupbyfactorservlvl <- as.ff(as.factor(paste(wdcMerged$station_id_index[],wdcMerged$tw[]))) 

agg_serv_lvl = binned_sum(wdcMerged$stocked_out, bin=wdcMerged$groupbyfactorservlvl)
agg_serv_lvl = as.data.frame(agg_serv_lvl)
agg_serv_lvl$groupbyfactorservlvl = row.names(agg_serv_lvl)
agg_serv_lvl$serv_lvl <- 1-(agg_serv_lvl$sum/agg_serv_lvl$count)  
length(which(agg_serv_lvl$serv_lvl==0))
agg <- agg_serv_lvl
agg_serv_lvl <- NULL

agg_serv_lvl_less_5_bikes = binned_sum(wdcMerged$less_5_bikes, bin=wdcMerged$groupbyfactorservlvl)
agg_serv_lvl_less_5_bikes = as.data.frame(agg_serv_lvl_less_5_bikes)
agg_serv_lvl_less_5_bikes$groupbyfactorservlvl = row.names(agg_serv_lvl_less_5_bikes)
agg_serv_lvl_less_5_bikes$serv_lvl <- 1-(agg_serv_lvl_less_5_bikes$sum/agg_serv_lvl_less_5_bikes$count)

agg_serv_lvl_less_5_bikes <- agg_serv_lvl_less_5_bikes[order(agg_serv_lvl_less_5_bikes$groupbyfactorservlvl),]

if(!identical(agg$groupbyfactorservlvl,agg_serv_lvl_less_5_bikes$groupbyfactorservlvl))     stop("agg_serv_lvl_less_5_bikes not in order") 

agg$serv_lvl_less_5_bikes <- agg_serv_lvl_less_5_bikes$serv_lvl
#     agg$serv_lvl_less_5_bikes_lag <- agg_serv_lvl_less_5_bikes$serv_lvl_lag
agg_serv_lvl_less_5_bikes <- NULL

agg_in_dem = binned_sum(wdcMerged$in_dem, bin=wdcMerged$groupbyfactorservlvl)
agg_in_dem = as.data.frame(agg_in_dem)
agg_in_dem$groupbyfactorservlvl = row.names(agg_in_dem)
agg_in_dem$in_dem_rate <- agg_in_dem$sum/agg_in_dem$count

agg_in_dem <- agg_in_dem[order(agg_in_dem$groupbyfactorservlvl),]

if(!identical(agg$groupbyfactorservlvl,agg_in_dem$groupbyfactorservlvl))     stop("agg_in_dem not in order") 

agg$in_dem_rate <- agg_in_dem$in_dem_rate
#     agg$serv_lvl_less_5_bikes_lag <- agg_serv_lvl_less_5_bikes$serv_lvl_lag
agg_in_dem <- NULL

agg_out_dem = binned_sum(wdcMerged$out_dem, bin=wdcMerged$groupbyfactorservlvl)
agg_out_dem = as.data.frame(agg_out_dem)
agg_out_dem$groupbyfactorservlvl = row.names(agg_out_dem)
agg_out_dem$out_dem_rate <- agg_out_dem$sum/agg_out_dem$count

agg_out_dem <- agg_out_dem[order(agg_out_dem$groupbyfactorservlvl),]

if(!identical(agg$groupbyfactorservlvl,agg_out_dem$groupbyfactorservlvl))     stop("agg_out_dem not in order") 

agg$out_dem_rate <- agg_out_dem$out_dem_rate
#     agg$serv_lvl_less_5_bikes_lag <- agg_serv_lvl_less_5_bikes$serv_lvl_lag
agg_out_dem <- NULL

agg$groupbyfactorservlvl <- as.factor(agg$groupbyfactorservlvl)
agg <- as.ffdf(agg[,c("groupbyfactorservlvl","serv_lvl","serv_lvl_less_5_bikes",
  "in_dem_rate","out_dem_rate")])
wdcMerged$dup <- duplicated(wdcMerged$groupbyfactorservlvl)
wdcMerged <- subset(wdcMerged,dup==FALSE)
wdcMerged <- merge(wdcMerged, agg, by="groupbyfactorservlvl")

agg <- NULL
wdcMerged <- as.data.frame(wdcMerged)

station_serv_lvl_effects <- wdcMerged  

wdcMerged <- wdcMerged_dem
#     wdcMerged <- subset(wdcMerged, week!=min_week)
wdcMerged$week_fac <- droplevels(wdcMerged$week_fac)    

#  	wdcMerged$week_tw <- paste0(wdcMerged$week,"_",wdcMerged$tw)
wdcMerged$groupbyfactor <- paste0(wdcMerged$week,"_",wdcMerged$tw,"_",
                                  wdcMerged$station_id_index)
wdcMerged <- wdcMerged[order(wdcMerged$groupbyfactor,wdcMerged$stocked_out,-wdcMerged$obs_weight),]
wdcMerged$index <- ave(wdcMerged$station_id_index, wdcMerged$groupbyfactor, FUN=function(x)c(1:length(x)))
#drop stocked out observation unless they are the only observation in certain groupbyfactor, i.e. at index 1
wdcMerged <- subset(wdcMerged, index==1 | stocked_out==F)
#wdcMerged <- subset(wdcMerged, index<=num_top_states_select)
points <- generate_integration_points()

#for each points, find out stations that are in locality and store them as a string
station_data <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
station_data <- station_data[order(station_data$station_id_index),]
points$local_stations <- rep(NA,nrow(points))
for(i in 1:nrow(points)) {
  lat1 = points$lat[i]
  lon1 = points$lon[i]
  dis_v <- latlondistance(lat1,lon1,station_data$lat,station_data$lon)  
  order_dis_v <- order(dis_v)
  
  points$local_stations[i] =paste(station_data$station_id_index[sort(order_dis_v[which(dis_v[order_dis_v[1:point_range]] <= max_walking_dis)])]
                                  , collapse="_")
}  


wdcMerged$tw_group <- paste0(wdcMerged$week,"_",wdcMerged$tw)
#  wdcMerged$tw_group <- paste0(wdcMerged$tw)
wdcMerged$tw_group_fac <- as.factor(wdcMerged$tw_group)
wdcMerged$out_dem_mean <- wdcMerged$out_dem_sum/wdcMerged$obs_weight
wdcMerged$st_tw <- as.factor(paste0(wdcMerged$station_id_index,"_",wdcMerged$tw))
#wdcMerged$st_tw_index <- as.numeric(wdcMerged$st_tw)
wdcMerged <- wdcMerged[order(wdcMerged$tw_group,wdcMerged$station_id_index, wdcMerged$sto_state_local),]
######################################################################

#market_share = 0.10
# tract_in <- tract
# points <- subset(points, tract %in% tract_in)
#     for(tract_in in tract) {
#       wdcMerged_tract <- wdcMerged[which(wdcMerged$tract == tract_in),]
#       #net_demand = max(as.numeric(by(wdcMerged_tract$out_dem_mean, wdcMerged_tract$tw_group_fac, FUN=sum)))
#       wdcMerged_tract <- subset(wdcMerged_tract,stocked_out==F)
#       #at each stationXtw_group_fac level create a out_dem_mean
#       wdcMerged_tract$st_tw_group_fac <- as.factor(paste(wdcMerged_tract$station_id_index, wdcMerged_tract$tw_group_fac))
#       wdcMerged_tract$out_dem_mean_dem <- ave(wdcMerged_tract$out_dem_sum, wdcMerged_tract$st_tw_group_fac, FUN=sum)
#       wdcMerged_tract$out_dem_mean_obs <- ave(wdcMerged_tract$obs_weight, wdcMerged_tract$st_tw_group_fac, FUN=sum)
#       wdcMerged_tract$out_dem_mean <- wdcMerged_tract$out_dem_mean_dem/wdcMerged_tract$out_dem_mean_obs
#       wdcMerged_tract$dup <- duplicated(wdcMerged_tract$st_tw_group_fac)
#       wdcMerged_tract <- subset(wdcMerged_tract, dup==F)
#       dem_sum <- as.numeric(by(wdcMerged_tract$out_dem_mean, wdcMerged_tract$tw_group_fac, FUN=sum))
#       net_demand <- max(dem_sum)      
#       idx <- which(points$tract %in% tract_in & points$type==1)
#       points$density[idx] = net_demand/market_share/length(idx)
#     }
points$weight <- points$density
points$density <- NULL #changing density to weight to track where density is being used
points$weight[which(points$type==2)] <- points$weight[which(points$type==2)]/max(points$weight[which(points$type==2)])
points$weight[which(points$type==4)] <- points$weight[which(points$type==4)]/max(points$weight[which(points$type==4)])
points <- points[order(points$type, points$lat, points$lon),]

#google places data
places_data <- read_googleplaces_data_dis_points()
# places_data <- subset(places_data, tract %in% tract_in)
# #places_cols_select <- c("lat","lon","cafe", "grocery_or_supermarket", "local_government_office")
# places_cols_select <- c("lat","lon","places_count","cafe", "grocery_or_supermarket","local_government_office","food")
# places_data <- places_data[,places_cols_select]
places_data <- places_data[order(places_data$lat, places_data$lon),]
if(!identical(round(points$lat[c(1:nrow(places_data))],4), round(places_data$lat,4))) stop("points lat dont match")
if(!identical(round(points$lon[c(1:nrow(places_data))],4), round(places_data$lon,4))) stop("points lon dont match")
places_data_temp <- places_data
places_data_temp$lat <- NULL
places_data_temp$lon <- NULL
places_colnames <- colnames(places_data_temp)
places_data_temp_full <- as.data.frame(matrix(0,nrow=nrow(points),ncol=ncol(places_data_temp)))
colnames(places_data_temp_full) <- places_colnames
places_data_temp_full[c(1:nrow(places_data_temp)),] <- places_data_temp
points <- cbind(points, places_data_temp_full)

# and convert to coefficient based 
station_data_tract <- unique(wdcMerged[,c("tract","station_id_index")])  

current_serv_lvl <- data.frame(station_id_index= station_serv_lvl_effects$station_id_index)
#current_serv_lvl$tw_group <- paste0(station_serv_lvl_effects$tw)
current_serv_lvl$tw <- station_serv_lvl_effects$tw
current_serv_lvl$serv_lvl <- station_serv_lvl_effects$serv_lvl_less_5_bikes
current_serv_lvl$st_tw <- as.factor(paste0(current_serv_lvl$station_id_index,"_",current_serv_lvl$tw))
#current_serv_lvl$st_tw_index <- as.numeric(current_serv_lvl$st_tw)
current_serv_lvl$instr_serv_lvl <- station_serv_lvl_effects$in_dem_rate
current_serv_lvl$in_dem_rate <- station_serv_lvl_effects$in_dem_rate
current_serv_lvl$out_dem_rate <- station_serv_lvl_effects$out_dem_rate
current_serv_lvl <- merge(current_serv_lvl,station_data_tract, by="station_id_index")
current_serv_lvl$tract_tw <- as.factor(paste0(current_serv_lvl$tract,"_",current_serv_lvl$tw))
#current_serv_lvl <- current_serv_lvl[order(current_serv_lvl$st_tw_index),]

#lagged in_dem_rate and out_dem_rate
current_serv_lvl <- current_serv_lvl[order(current_serv_lvl$station_id_index,
                                           current_serv_lvl$tw),]
current_serv_lvl$in_dem_rate_lagtw <- ave(current_serv_lvl$in_dem_rate,
    current_serv_lvl$station_id_index, FUN=function(x) {return(c(x[-1],x[1]))})
current_serv_lvl$out_dem_rate_lagtw <- ave(current_serv_lvl$out_dem_rate,
    current_serv_lvl$station_id_index, FUN=function(x) {return(c(x[-1],x[1]))})
current_serv_lvl$diff_dem_rate_lagtw <- current_serv_lvl$in_dem_rate_lagtw - 
  current_serv_lvl$out_dem_rate_lagtw  
current_serv_lvl$diffdummy_dem_rate_lagtw <- ((current_serv_lvl$in_dem_rate_lagtw - 
  current_serv_lvl$out_dem_rate_lagtw)>0)

st_tw_table <- current_serv_lvl[,c("station_id_index","tw")]
st_tw_table <- st_tw_table[order(st_tw_table$tw, st_tw_table$station_id_index),]
st_tw_table$st_tw_index <- c(1:nrow(st_tw_table))
wdcMerged <- merge(wdcMerged, st_tw_table, by=c("station_id_index","tw"))
current_serv_lvl <- merge(current_serv_lvl, st_tw_table, by=c("station_id_index","tw"))

x0_start <<- NULL
prev_theta <<- NULL

#add density attributes for writing moment conditions
wdcMerged <- get_local_attributes_st_state_alt_tw(wdcMerged, points,c(0:3)) #assining 0030-0530 as metro off hours 


#if(!identical(unique(wdcMerged$tw_group),unique(current_serv_lvl$tw_group))) stop("error")
#source("reg_results_GMM_discoef_2.89.R")
wdcMerged$out_dem_sum[which(wdcMerged$out_dem_sum<0.0001 & 
                              wdcMerged$stocked_out==FALSE)] <- 0.0001

wdcMerged <- wdcMerged[order(wdcMerged$tw_group,wdcMerged$station_id_index,wdcMerged$sto_state_local),]
current_serv_lvl <- current_serv_lvl[order(current_serv_lvl$tw,current_serv_lvl$station_id_index),]
user_serv_lvl <- current_serv_lvl

#add serv_lvl to wdcMerged
if(!identical(c(1:nrow(current_serv_lvl)),current_serv_lvl$st_tw_index)) stop("error Line 320")
wdcMerged$serv_lvl <- current_serv_lvl$serv_lvl[wdcMerged$st_tw_index]
wdcMerged <- droplevels(wdcMerged)
gc()



###save commands
# save(wdcMerged, file="wdcMerged_lin_dis1to20_aggmonths_bytw_1hr_hyperlocalstate_withreallocthresh_preweatherreg_stkoutthresh_five_pr6.RData")
# save(current_serv_lvl, file="current_serv_lvl_lin_dis1to20_aggmonths_bytw_1hr_hyperlocalstate_withreallocthresh_preweatherreg_stkoutthresh_five_pr6.RData")
# save(user_serv_lvl, file="user_serv_lvl_lin_dis1to20_aggmonths_bytw_1hr_hyperlocalstate_withreallocthresh_preweatherreg_stkoutthresh_five_pr6.RData")
# save(points, file="points_lin_dis1to20_aggmonths_bytw_1hr_hyperlocalstate_withreallocthresh_preweatherreg_stkoutthresh_five_pr6.RData")

# save(wdcMerged, file="wdcMerged_lin_dis1to20_aggmonths_bytw_1hr_hyperlocalstate_withreallocthresh_preweatherreg_stkoutthresh_five_pr3.RData")
# save(current_serv_lvl, file="current_serv_lvl_lin_dis1to20_aggmonths_bytw_1hr_hyperlocalstate_withreallocthresh_preweatherreg_stkoutthresh_five_pr3.RData")
# save(user_serv_lvl, file="user_serv_lvl_lin_dis1to20_aggmonths_bytw_1hr_hyperlocalstate_withreallocthresh_preweatherreg_stkoutthresh_five_pr3.RData")
# save(points, file="points_lin_dis1to20_aggmonths_bytw_1hr_hyperlocalstate_withreallocthresh_preweatherreg_stkoutthresh_five_pr3.RData")

