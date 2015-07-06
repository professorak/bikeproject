source("constants_mw6.R")
#v3.0 two stockstate observations per stations, if results are consistent, can experiment with other fixed effects etc.
point_range <- 3
no_bootstrap_runs <- 1
reg_runs_save_file <- "reg_runs_discoef_2.89_8states_mw6_bootstrap.txt"


source("eval_func_4_GMM_new_temp_2.89.R")
source("data_prep_functions_ffdf_limstostate.R")
source("data_prep_inte_points_limstostate.R")
source("eval_func_2_cpp_cntrt_map_2.86_2.R")
source("eval_func_3_cpp_new_2.86_2.R")

#max_walking_dis = 0.1
tract_list <- list(1:10)
for(tract in tract_list) {

  st_list = select_stations_in_tract(c(tract))
  
  #temporary - load st_6_7_8 file, should be able to generate this realtime, 
  save_dir <- paste0(csv_dir,"/ffdb/Paris/st_5_8")
  load(save_dir)
  st_list <- intersect(st_5_8,st_list)
  
  filelist_all <- c("Paris/ind_paris_5.csv","Paris/ind_paris_6.csv",
                    "Paris/ind_paris_7.csv","Paris/ind_paris_8.csv")
  filelist_all_mons <- c(5:8)
#   filelist_all <- c("Paris/ind_paris_7.csv")
#   filelist_all_mons <- c(7)
  #filelist_servlvl <- c("Paris/ind_paris_4.csv","Paris/ind_paris_5.csv")
  wdcMerged <- c()
  for( i in 1:length(filelist_all_mons)) {
    filelist <- filelist_all[i]
    mon <- filelist_all_mons[i]
    save_file <- paste0(csv_dir,"/ffdb/wdcMerged_bootstrap_limstostate_5_8_den",
          point_range,"/mw6_tract_",paste0(tract,collapse="_"),"_mon_",mon)  
    print(save_file)
    #load files
    wdcMerged_bootstrap <- reg_run_bootstrap_simple(NULL,save_file)
    dim(wdcMerged_bootstrap)
    wdcMerged <- rbind(wdcMerged,wdcMerged_bootstrap)
    rm(wdcMerged_bootstrap)
  }  

  wdcMerged_save <- wdcMerged
  #wdcMerged_save <- subset(wdcMerged_save, tw>=6 & tw <10)
  wdcMerged_save <- droplevels(wdcMerged_save)
  
  
  #v0_vec <- generate_v0() #necessary to generate outside as seed is otherwise set inside, bad
  v0_vec <- c(0)
  set.seed(2)
  i=1
  for(i in 1:no_bootstrap_runs) {
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
    wdcMerged$week <- wdcMerged$month
    wdcMerged$tw <- floor(wdcMerged$tw/4)
    wdcMerged$tw_fac <-  as.ff(as.character(wdcMerged$tw))
    wdcMerged_save_i <- wdcMerged
    #compute demand at the level of st,week,tw,sto_state
    
    #wdcMerged <- subset(wdcMerged, stocked_out==FALSE)
    #since aggdata aggregates over stid_twg_locstate_fac, 
    wdcMerged$stid_twg_locstate_fac <-  as.ff(as.factor(paste(wdcMerged$station_id_index[],
                    wdcMerged$week[],wdcMerged$tw[],wdcMerged$sto_state_local[]))) 
    wdcMerged$stid_twg_locstate_fac <- droplevels(wdcMerged$stid_twg_locstate_fac)
    
    wdcMerged <- aggdata(wdcMerged)
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
    wdcMerged$groupbyfactorservlvl <- as.ff(as.factor(paste(wdcMerged$station_id_index[],wdcMerged$tw[],
                                                            wdcMerged$week[]))) 
    
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
    
    groupbyfactorservlvl_data <- wdcMerged[,c("groupbyfactorservlvl",
                                              "station_id_index","tw","week")]
    groupbyfactorservlvl_data$dup <- duplicated(groupbyfactorservlvl_data$groupbyfactorservlvl) 
    groupbyfactorservlvl_data <- subset(groupbyfactorservlvl_data, dup==F)
    agg_serv_lvl_less_5_bikes <- merge(agg_serv_lvl_less_5_bikes, groupbyfactorservlvl_data,
                                       by="groupbyfactorservlvl")

    agg_serv_lvl_less_5_bikes <- agg_serv_lvl_less_5_bikes[
      order(agg_serv_lvl_less_5_bikes$station_id_index, agg_serv_lvl_less_5_bikes$tw,
            agg_serv_lvl_less_5_bikes$week),]
    agg_serv_lvl_less_5_bikes$serv_lvl_lag <- c(NA,agg_serv_lvl_less_5_bikes$serv_lvl[-nrow(agg_serv_lvl_less_5_bikes)])
    min_week <- min(agg_serv_lvl_less_5_bikes$week)
    agg_serv_lvl_less_5_bikes$serv_lvl_lag[which(agg_serv_lvl_less_5_bikes$week==min_week)] = NA
    agg_serv_lvl_less_5_bikes <- agg_serv_lvl_less_5_bikes[order(agg_serv_lvl_less_5_bikes$groupbyfactorservlvl),]
    
    if(!identical(agg$groupbyfactorservlvl,agg_serv_lvl_less_5_bikes$groupbyfactorservlvl))     stop("agg_serv_lvl_less_5_bikes not in order") 
      
    agg$serv_lvl_less_5_bikes <- agg_serv_lvl_less_5_bikes$serv_lvl
    agg$serv_lvl_less_5_bikes_lag <- agg_serv_lvl_less_5_bikes$serv_lvl_lag
    agg_serv_lvl_less_5_bikes <- NULL
    
    agg$groupbyfactorservlvl <- as.factor(agg$groupbyfactorservlvl)
    agg <- as.ffdf(agg[,c("groupbyfactorservlvl","serv_lvl","serv_lvl_less_5_bikes"
                          ,"serv_lvl_less_5_bikes_lag")])
    wdcMerged$dup <- duplicated(wdcMerged$groupbyfactorservlvl)
    wdcMerged <- subset(wdcMerged,dup==FALSE)
    wdcMerged <- merge(wdcMerged, agg, by="groupbyfactorservlvl")
    
    agg <- NULL
    wdcMerged <- as.data.frame(wdcMerged)
    wdcMerged <- subset(wdcMerged, week!=min_week)
    
  	station_serv_lvl_effects <- wdcMerged	
	
    wdcMerged <- wdcMerged_dem
    wdcMerged <- subset(wdcMerged, week!=min_week)
    wdcMerged$week_fac <- droplevels(wdcMerged$week_fac)    
    
  	wdcMerged$week_tw <- paste0(wdcMerged$week,"_",wdcMerged$tw)
  	wdcMerged$groupbyfactor <- paste0(wdcMerged$week,"_",wdcMerged$tw,"_",
                                      wdcMerged$station_id_index)
    wdcMerged <- wdcMerged[order(wdcMerged$groupbyfactor,wdcMerged$stocked_out,-wdcMerged$obs_weight),]
    wdcMerged$index <- ave(wdcMerged$station_id_index, wdcMerged$groupbyfactor, FUN=function(x)c(1:length(x)))
    wdcMerged <- subset(wdcMerged, index<=8)
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
  wdcMerged$tw_group_fac <- as.factor(wdcMerged$tw_group)
  wdcMerged$out_dem_mean <- wdcMerged$out_dem_sum/wdcMerged$obs_weight
	wdcMerged <- wdcMerged[order(wdcMerged$week_tw,wdcMerged$station_id_index),]
    ######################################################################
    
    #market_share = 0.10
    tract_in <- tract
    points <- subset(points, tract %in% tract_in)
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
  
    # and convert to coefficient based 
  
    current_serv_lvl <- data.frame(station_id_index= station_serv_lvl_effects$station_id_index)
    current_serv_lvl$tw_group <- paste0(station_serv_lvl_effects$week,"_",station_serv_lvl_effects$tw)
    current_serv_lvl$serv_lvl <- station_serv_lvl_effects$serv_lvl_less_5_bikes_lag
    current_serv_lvl$stid_twg_fac <- as.factor(paste(current_serv_lvl$station_id_index, current_serv_lvl$tw_group))
    current_serv_lvl$instr_serv_lvl <- current_serv_lvl$serv_lvl
	  current_serv_lvl <- current_serv_lvl[order(current_serv_lvl$tw_group,current_serv_lvl$station_id_index),]
#     stop("failing due to less than no_st observations in most of the tw_group,
#          i.e. some of the stations are always stocked in particular tw_group.") 
    sttw_no <- ave(wdcMerged$station_id_index, wdcMerged$tw_group, 
                   FUN=function(x) length(unique(x)))
  
    x0_start <<- NULL
    prev_theta <<- NULL
    if(!identical(unique(wdcMerged$tw_group),unique(current_serv_lvl$tw_group))) stop("error")
    #source("reg_results_GMM_discoef_2.89.R")
    wdcMerged$out_dem_sum[which(wdcMerged$out_dem_sum<=0.0001 & 
                                  wdcMerged$stocked_out==FALSE)] <- 0.0001
    
    wdcMerged <- wdcMerged[order(wdcMerged$tw_group,wdcMerged$station_id_index,wdcMerged$sto_state_local),]
    current_serv_lvl <- current_serv_lvl[order(current_serv_lvl$tw_group,current_serv_lvl$station_id_index),]
    user_serv_lvl <- current_serv_lvl
    gc()
  }
}


