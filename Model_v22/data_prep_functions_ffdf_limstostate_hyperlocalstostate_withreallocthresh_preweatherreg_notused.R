source("data_prep_functions_ffdf_limstostate_withreallocthresh_preweatherreg.R")

reg_ready_data <- function(wdcMergedin=NULL,filelist) {
  #filter only weekdays, generate lat,lon data and stockout state
  
  wdcMerged <- readdata(wdcMergedin,filelist)
  writefile <- paste0(csv_dir,"/filelogs/observationnumberlog.txt")
  write(paste("after readdata",nrow(wdcMerged)),file=writefile, append=T)
  write(paste("after readdata trips",sum(wdcMerged$out_dem)),file=writefile, append=T)
  
  #keeping only weekdays  
  wdcMerged <- subset(wdcMerged, wday >=1 & wday <=5)
  writefile <- paste0(csv_dir,"/filelogs/observationnumberlog.txt")
  write(paste("after only weekdays",nrow(wdcMerged)),file=writefile, append=T)
  write(paste("after only weekdays trips",sum(wdcMerged$out_dem)),file=writefile, append=T)
  
  #compute the local stations to a station
  station_data <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
  #order station_data by "station_id_index" and check that is it same as row numbers
  station_data <- station_data[order(station_data$station_id_index),]
  if(!identical(as.integer(station_data$station_id_index), c(1:nrow(station_data)))) {
    stop("ERROR 5489: station_data not in station_id_index order")
  }  
  points <- generate_integration_points(station_data)  
  station_data$local_points <- rep("",nrow(station_data))
  points$local_stations <- rep(NA,nrow(points))
  for(i in 1:nrow(points)) {
    lat1 = points$lat[i]
    lon1 = points$lon[i]
    dis_v <- latlondistance(lat1,lon1,station_data$lat,station_data$lon)  
    order_dis_v <- order(dis_v)
    
    points$local_stations[i] =paste(station_data$station_id_index[sort(order_dis_v[which(dis_v[order_dis_v[1:point_range]] <= max_walking_dis)])]
                                    , collapse="_")
    for(stid in sort(order_dis_v[which(dis_v[order_dis_v[1:point_range]] <= max_walking_dis)])) {
      #add local points to stid
      if(station_data$local_points[stid]!="") {
        station_data$local_points[stid] <- paste(c(station_data$local_points[stid],i),collapse="_")
      } else {
        station_data$local_points[stid] <- i
      }    
    }
  }  
  
  #for each station, go to each of the local points, and merge all their local stations.
  station_data$local_stations <- NA
  #station_data$no_local_stations <- NA
  for(stid in 1:nrow(station_data)) {
    points_list <- as.integer(strsplit(as.character(station_data$local_points[stid]),split="_")[[1]])
    st_list <- c()
    for(point_id in points_list) {
      st_list <- c(st_list,as.integer(strsplit(as.character(points$local_stations[point_id]),split="_")[[1]]))
    }
    st_list <- sort(unique(st_list))
    station_data$local_stations[stid] <- paste(st_list, collapse="_")
    #station_data$no_local_stations[stid] <- length(st_list)
  }
  
  station_data$hyperlocal_stations <- NA
  station_data$nonhyperlocal_stations <- NA
  #generate hyperlocal_stations and nonhyperlocal_stations in station_data from among local_stations  
  for(stid in 1:nrow(station_data)) {
    #split the local stations, find the closest #hyperlocal_range of stations
    st_list <- as.integer(strsplit(as.character(station_data$local_stations[stid]),split="_")[[1]])
    dis_v <- latlondistance(station_data$lat[stid], station_data$lon[stid], 
                  station_data$lat[st_list], station_data$lon[st_list])
    
    station_data$hyperlocal_stations[stid] <- paste(st_list[order(dis_v)][c(1:min(hyperlocal_range,length(dis_v)))], 
                                                    collapse="_")
    station_data$nonhyperlocal_stations[stid] <- paste(st_list[order(dis_v)][-c(1:min(hyperlocal_range,length(dis_v)))], 
                                                    collapse="_")
    #reassignign local_stations so that the initial part is hyperlocal_stations and rest is nonhyperlocal_stations
    station_data$local_stations[stid] <- paste(c(station_data$hyperlocal_stations[stid], station_data$nonhyperlocal_stations[stid]),
                                               collapse="_")
  }
  
  station_data <- station_data[,c("station_id_index","local_stations")]
  #station_data <- as.ffdf(station_data)
  station_data <- ffdf(station_id_index=as.ff(station_data$station_id_index),
    local_stations=as.ff(as.factor(station_data$local_stations)))

  wdcMerged <- merge(wdcMerged,station_data,by="station_id_index")
  
  
  #generate so_state for each station
  #sort by n_time_int, station_id_index
  wdcMerged <- wdcMerged[fforder(wdcMerged$n_time_int, wdcMerged$station_id_index),]
  wdcMerged$groupbyfactor <- as.character(wdcMerged$n_time_int)
  
  agg <- ffdfdply(wdcMerged, split=wdcMerged$groupbyfactor, FUN=function(x){
    x <- as.data.table(x)
    result <- x[,list(station_id_index,
                sto_state_local = gen_sto_state_hyperlocal_char(stocked_out,local_stations,hyperlocal_range)),
                by = list(groupbyfactor)]
    result <- as.data.frame(result)
    result
  })
  colnames(agg) <- gsub("groupbyfactor", "n_time_int", colnames(agg)) #replace groupbyfactor by n_time_int
  
  wdcMerged$groupbyfactor <- NULL
  
  wdcMerged <- merge(wdcMerged, agg, by=c("n_time_int","station_id_index"))

  
#   #renumber sto_state
#   sto_state_vec <- unique(agg$sto_state)
#   sto_state_vec <- sto_state_vec[fforder(sto_state_vec)]
#   sto_state_vec <- ffdf(sto_state_index=as.ff(c(1:length(sto_state_vec))),
#                           sto_state=sto_state_vec)
#   agg <- merge(agg,sto_state_vec,by="sto_state",all.x=TRUE) #it doesnt seem to support
#     #all.y=TRUE, says only innner joins allowed.
#   wdcMerged = merge(wdcMerged,agg,by="n_time_int",all.x=TRUE)  
  
  #spell length
  #state_vec <- gen_local_state_spell_length(wdcMerged)  
  #wdcMerged <- merge(wdcMerged, state_vec, by=c("n_time_int","station_id_index"),all.x=TRUE) 
  wdcMerged <- wdcMerged[fforder(wdcMerged$station_id,wdcMerged$n_time_int),]
  return(wdcMerged)  
}


