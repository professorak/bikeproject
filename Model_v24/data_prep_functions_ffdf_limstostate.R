#select a sample of stations 
require(ffbase)
require(data.table)
#select a set of stations around a given lat and lon
select_near_stations <- function(med_lat=48.8582, med_lon=2.2945) {
  #48.8582, 2.2945 is default eifel tower location  
  wdcLatLon <- readcsvlatlon()
  dis_v = latlondistance(med_lat, med_lon, wdcLatLon$lat, wdcLatLon$lon)
  st_list = wdcLatLon$station_id[which(dis_v<1.3)]
  return(st_list)
}

#select a set of stations around a given lat and lon
select_stations_in_tract <- function(tract_list=c(7,8)) {
  wdcLatLon <- readcsvlatlon()
  wdcCenTrParsed <- readtractboundaryfile()
  tract <- assign_points_tract(wdcLatLon$lat,wdcLatLon$lon,wdcCenTrParsed)
  st_list = wdcLatLon$station_id[which(tract %in% tract_list)]
  return(st_list)
}


readcsvlatlon <- function() {  
  filename <- file.path(csv_dir,st_id_file, fsep = .Platform$file.sep)
  wdcLatLonRaw <- read.csv(file=filename,head=TRUE,sep=",")	
  wdcLatLon <- wdcLatLonRaw[,c("id","lat","lon")]
  colnames(wdcLatLon) <- c("station_id","lat","lon")
  print(paste0("read ",filename))
  return (wdcLatLon)
}


readcsvdemand <- function(wdcRaw=NULL,filelist) {
  wdcRawOrg_all <- c()
  if(is.null(wdcRaw)) {
    for(file in filelist) {
      save_dir <- paste0(csv_dir,"/ffdb/",file)
      ######
      if(!file.exists(save_dir)) {
        filename <- file.path(csv_dir,file, fsep = .Platform$file.sep)
        setClass("myDate")
        setAs("character","myDate", function(from) as.POSIXct(from, format='%Y-%m-%d %H:%M:%S',tz="Europe/London") )
        
        wdcRaw <- read.csv.ffdf(file=filename,head=FALSE,sep=",",
                                nrows=-1,colClasses=c(V1="integer",V2="myDate",
                                                      V3="integer",V4="integer",V5="integer"))  
        colnames(wdcRaw) <- c("station_id" ,"time" ,"bikes" ,"spaces_available" ,"unbalanced")
        save_dir <- paste0(csv_dir,"/ffdb/",file)
        save.ffdf(wdcRaw,dir=save_dir)
      }
      ###
      load.ffdf(save_dir)
      wdcRawOrg_city <- subset(wdcRaw, station_id %in% st_list & bikes!=-1)
      #wdcRawOrg_city <- as.data.frame(wdcRawOrg_city)
      wdcRawOrg_all <- rbind(wdcRawOrg_all,wdcRawOrg_city)
    }
  } else {
    print("in else")
    wdcRawOrg_city <- subset(wdcRaw, station_id %in% st_list & bikes!=-1)
    #wdcRawOrg_city <- as.data.frame(wdcRawOrg_city)
    wdcRawOrg_all <- rbind(wdcRawOrg_all,wdcRawOrg_city)
  }
  wdcRawOrg_all$n_time_int <- as.ff(as.numeric(wdcRawOrg_all$time[]))
  #remove duplicates from combining datasets
  wdcRawOrg_all$dup <- duplicated.ffdf(as.ffdf(wdcRawOrg_all[,c("station_id","n_time","bikes")]))
  wdcRawOrg <- subset.ffdf(wdcRawOrg_all,dup == "FALSE")
  wdcRawOrg$dup <- NULL
  wdcRawOrg_all <- NULL
  
  #doing the following way as POSIXlt cannot be as ff vector
  wdcRawOrg$tw <- as.ff(as.POSIXlt(wdcRawOrg$time[],"%Y-%m-%d %H:%M:%S",
                                   tz="Europe/Paris")$hour)
  wdcRawOrg$seconds <- as.ff(as.numeric(wdcRawOrg$time[])) %% 3600
  wdcRawOrg$month <- as.ff(as.POSIXlt(wdcRawOrg$time[],"%Y-%m-%d %H:%M:%S",
                                      tz="Europe/Paris")$mon + 1)
  wdcRawOrg$day <- as.ff(as.POSIXlt(wdcRawOrg$time[],"%Y-%m-%d %H:%M:%S",
                                    tz="Europe/Paris")$yday)
  wdcRawOrg$wday <- as.ff(as.POSIXlt(wdcRawOrg$time[],"%Y-%m-%d %H:%M:%S",
                                     tz="Europe/Paris")$wday)
  #by station id, generate time diff and change bikes
  #order by station, time 
  wdcRawOrg$groupbyfactor <- as.character(wdcRawOrg$station_id)
  wdcRawOrg <- wdcRawOrg[ffdforder(wdcRawOrg[c("station_id", "n_time_int")]),]  
  
  diff_ntimeint <- ffdfdply(wdcRawOrg[c("station_id", "n_time_int","bikes","groupbyfactor")], 
                            wdcRawOrg$groupbyfactor, #BATCHBYTES = 167772160, #10 times default = 0.15GB
                            FUN=function(x){       
                              x$groupbyfactor <- droplevels(x$groupbyfactor)
                              x <- split(x, x$groupbyfactor)
                              x <- lapply(x, FUN=function(onlyonedeal){
                                ret <- c(diff(onlyonedeal$n_time_int),NA)
                                ret2 <- c(diff(onlyonedeal$bikes),NA)
                                ret3 <- c(nrow(onlyonedeal):1)
                                data.frame(time_diff= ret,change_bikes=ret2,index=ret3,
                                           station_id=onlyonedeal$station_id, n_time_int=onlyonedeal$n_time_int)
                              })
                              x <- do.call(rbind, x)      
                              
                              return(x)
                            }
  )
  #check if wdcRawOrg and diff_ntimeint are in same order so they could be row combined
  #diff_ntimeint <- diff_ntimeint[ffdforder(diff_ntimeint[c("station_id", "n_time_int")]),]  
  #should sort diff_ntimeint on n_time_int as the assumption is n_time_int sorting will be maintained
  #and logic of computing diff would be intact
  #my guess is the order of station_id is not maintained because internally above function might be 
  #doing some parallel processing for each split and doesnt keep the order
  diff_ntimeint <- diff_ntimeint[ffdforder(diff_ntimeint[c("station_id")]),]  
  if(!identical(diff_ntimeint$station_id[],
                wdcRawOrg$station_id[])) stop("not in order")
  if(!identical(diff_ntimeint$n_time_int[],
                wdcRawOrg$n_time_int[])) stop("not in order")

  wdcRawOrg$time_diff <- diff_ntimeint$time_diff
  wdcRawOrg$change_bikes <- diff_ntimeint$change_bikes
  wdcRawOrg$index <- diff_ntimeint$index
  wdcRawOrg$change_bikes_avg <- wdcRawOrg$change_bikes/wdcRawOrg$time_diff*120
  
  #drop the last row for every station_id	
  wdcRawOrg <- subset(wdcRawOrg, index!=1)
  wdcRawOrg$index <- NULL
  
  #stocked out flag
  wdcRawOrg <- transform(wdcRawOrg,stocked_out=bikes==0)
  wdcRawOrg <- transform(wdcRawOrg,stocked_full=spaces_available==0)
  
  #calculate outgoing demand
  wdcRawOrg <- transform(wdcRawOrg,
                out_dem = floor(-change_bikes_avg + abs(change_bikes_avg))/2,
                in_dem = floor(change_bikes_avg + abs(change_bikes_avg))/2 )
  
  return(wdcRawOrg)
}

readdata <- function(wdcRaw=NULL,filelist) {

  wdcRawOrg <- readcsvdemand(wdcRaw,filelist) 
  wdcRawOrg <- subset(wdcRawOrg, month==mon)
  writefile <- paste0(csv_dir,"/filelogs/observationnumberlog.txt")
  write("",file=writefile, append=T)
  write(paste("filelist: ",filelist),file=writefile, append=T)
  write(paste("after readcsvdemand",nrow(wdcRawOrg)),file=writefile, append=T)
  write(paste("after readcsvdemand trips",sum(wdcRawOrg$out_dem)),file=writefile, append=T)

  #create station id and day indexes to be used by fixed effects
  station_id_list <- unique(wdcRawOrg$station_id) 
  wdcstidindex <- data.frame(cbind(station_id_list[],1:length(station_id_list)))
  colnames(wdcstidindex) <- c("station_id","station_id_index")
  wdcRawOrg <- merge(wdcRawOrg,as.ffdf(wdcstidindex), by="station_id")
  
  #fill up the missing n_time_int with stockouts for missing entries 
  #generate list of n_time_ints
  #merge with wdcRawOrg and create complete join.
  #replace NA with stockouts
  wdcRawOrg$dup <- duplicated((wdcRawOrg$n_time_int))
  n_time_list <- subset.ffdf(wdcRawOrg,dup == "FALSE")
  n_time_list <- n_time_list[c("time","tw","seconds","month","day","wday","time_diff","n_time_int")]
  wdcRawOrg <- wdcRawOrg[c("station_id","bikes","spaces_available","unbalanced",
                "change_bikes","change_bikes_avg","stocked_out",
                "stocked_full","out_dem","in_dem","n_time_int","station_id_index")]
  
  expand_idx <- ffrep.int(1:nrow(n_time_list),length(st_list))
  n_time_list_ffdf <- n_time_list[expand_idx,]
  n_time_list_ffdf$station_id_index <- as.ff(1:nrow(n_time_list_ffdf))
  n_time_list_ffdf <- 
    transform(n_time_list_ffdf, station_id_index=ceiling(station_id_index/nrow(n_time_list)))
  
  wdcRawOrg <- merge(n_time_list_ffdf,wdcRawOrg,by=c("station_id_index","n_time_int"),all.x=TRUE)
  #fill up stocked_out with TRUE and out_dem with 0 for all NA entries
  na_list <- ffwhich(wdcRawOrg,is.na(wdcRawOrg$stocked_out)==TRUE) 
  wdcRawOrg[na_list,"stocked_out"] <- as.ff(!as.logical(1:length(na_list)))
  wdcRawOrg[na_list,"stocked_full"] <- as.ff(!as.logical(1:length(na_list)))
  wdcRawOrg[na_list,"out_dem"] <- ffrep.int(0,length(na_list))
  
  st_index_pairs <- unique(wdcRawOrg[c("station_id","station_id_index")])
  st_index_pairs <- st_index_pairs[ffwhich(st_index_pairs,
                                           !is.na(st_index_pairs$station_id)),]
  
  wdcLatLon <- readcsvlatlon()
  wdcLatLon <- merge(wdcLatLon,st_index_pairs,by="station_id")
  #merge demand and latlon file
  wdcRawOrg$station_id <- NULL
  wdcMerged <- merge(wdcRawOrg, as.ffdf(wdcLatLon), by ="station_id_index", all.x=TRUE)  
  #check for any NA's - indication of unmatched station ids in merge
  if(length(ffwhich(wdcMerged,is.na(wdcMerged$lat))) >0 ) stop("unmatched station ids")

  writefile <- paste0(csv_dir,"/filelogs/observationnumberlog.txt")
  write(paste("after expanding",nrow(wdcMerged)),file=writefile, append=T)
  write(paste("after expanding trips",sum(wdcMerged$out_dem)),file=writefile, append=T)

  #report number of observations for each station as some of the stations might have 
  #a lot of missing data and a lot of observations are getting eliminated in next 
  #steps due to this
  #a=by(wdcMerged$station_id,wdcMerged$station_id,FUN=length)
  #print("no of observations of stations:")
  #print(as.numeric(a))
  ##drop observation which dont have for entire range of stations for given time  
  no_st <- length(unique(wdcMerged$station_id))

  wdcMerged$groupbyfactor <- as.character(wdcMerged$n_time_int)
  agg <- binned_sum(x=wdcMerged$station_id, bin=wdcMerged$groupbyfactor)
  aggg <- as.data.frame(agg,n_time_int=as.numeric(row.names(agg))) 
  aggg$groupbyfactor <- as.character(row.names(agg))
  d <- aggg[which(aggg$count!=no_st),"groupbyfactor"]
  wdcMerged = subset(wdcMerged,!(groupbyfactor %in% d))
  
  
  #order wdcMerged by days, st, n_time. 
  wdcMerged <- wdcMerged[ffdforder(wdcMerged[c("day","station_id","n_time_int")]),]
  
  #consistency checks
  #out demand at stocked out points 
  if(sum(wdcMerged$out_dem[ffwhich(wdcMerged,wdcMerged$stocked_out==TRUE)]) != 0) stop("out dem non-zero at stockouts")
  #if(sum(wdcMerged$in_dem[ffwhich(wdcMerged,wdcMerged$stocked_full==TRUE)]) != 0) stop("in dem non-zero at stockfulls")
  print("not doing stocked full consistency check for now")
  #there are some instances where the consistency checks of incoming demand fails as 
  #false stockfulls appear at many places where spaces available suddenly shrink to 0 from 1 or 2
  
  
  #there seems to be inconsistencty sometimes in bikes+spaces_available count and not systematic
  
  return (wdcMerged )  
}

gen_local_state_spell_length <- function(wdcMerged) { 
#not ffdf optimized
  state_vec_all <- c()    
  for(st_id in 1:max(wdcMerged$station_id_index)) {
    state_vec <- wdcMerged[ffwhich(wdcMerged,station_id_index==st_id),c("n_time_int","sto_state_local")]
    state_vec <- state_vec[fforder(state_vec$n_time_int),]
    state_vec$continuity_index <- ffrep.int(NA, times=nrow(state_vec))
    state_vec$continuity_index[1] = 1
    state_vec <- as.data.frame(state_vec[])
    state_vec$sto_state_local_prev <- state_vec$sto_state_local
    state_vec$sto_state_local_prev[2:nrow(state_vec)] <- 
        state_vec$sto_state_local[1:(nrow(state_vec)-1)]
    state_vec$continuity_index <- state_vec$sto_state_local!=state_vec$sto_state_local_prev 
    state_vec$continuity_index <- as.integer(state_vec$continuity_index)
    state_vec$continuity_index[1] = 1
    state_vec$continuity_index <- cumsum(state_vec$continuity_index)
    
    state_vec$spell_length <- ave(state_vec$continuity_index, state_vec$continuity_index, FUN=length)
    state_vec <- state_vec[,c("n_time_int","spell_length")]
    state_vec$station_id_index <- st_id
    state_vec_all <- ffdfappend(state_vec_all,as.ffdf(state_vec),adjustvmode=F)
  }  
  return(state_vec_all)
}


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
  station_data$no_local_stations <- NA
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
  
  station_data <- station_data[,c("station_id_index","local_stations")]
  station_data <- ffdf(station_id_index=as.ff(station_data$station_id_index),
                       local_stations=as.ff(as.factor(station_data$local_stations)))
  wdcMerged <- merge(wdcMerged,station_data,by="station_id_index")
  
  
  #generate so_state for each station
  #sort by n_time_int, station_id_index
  wdcMerged <- wdcMerged[fforder(wdcMerged$n_time_int, wdcMerged$station_id_index),]
  wdcMerged$groupbyfactor <- as.character(wdcMerged$n_time_int)
  
  agg <- ffdfdply(wdcMerged, split=wdcMerged$groupbyfactor, FUN=function(x){
    x <- as.data.table(x)
    result <- x[,list(station_id_index,sto_state_local = gen_sto_state_local_char(stocked_out,local_stations)), by = list(groupbyfactor)]
    result <- as.data.frame(result)
    result
  })

  wdcMerged$groupbyfactor <- NULL
  colnames(agg) <- c("n_time_int","station_id_index","sto_state_local")
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

gen_serv_lvls <- function(wdcMerged,tw_low=8,tw_high=19) {
  #service levels are generated for all states and not just observations
  #which are used for rest of computations, since they tend to 1 when filtered 
  #with those observations. 
  #only filteration i will do is for time windows.
  list_all <- ffwhich(wdcMerged, tw >=tw_low & tw <=tw_high)
  
  wdcMerged_listall <- wdcMerged[list_all,]
  
  wdcMerged_listall$groupbyfactor <- as.ff(as.factor(paste(wdcMerged_listall$station_id_index[], 
                          wdcMerged_listall$tw_group[]))) 
  
  #wdcMerged_listall$groupbyfactor <- as.character(wdcMerged_listall$station_id_index)
  agg <- binned_sum(x=wdcMerged_listall$stocked_out, bin=wdcMerged_listall$groupbyfactor, 
                    nbins = length(levels(wdcMerged_listall$groupbyfactor)))
  #extract station_id and tw_group from row names of agg
  st_twg = lapply(rownames(agg), FUN= function(x) strsplit(as.character(x),split=" ")[[1]])
  st_twg = (do.call(rbind,st_twg))
  
  serv_lvl_arr <- data.frame(station_id_index=as.numeric(st_twg[,1]), tw_group=as.numeric(st_twg[,2]),
                             serv_lvl=agg[,"sum"]/agg[,"count"],row.names = NULL)
  serv_lvl_arr$serv_lvl <- 1-serv_lvl_arr$serv_lvl
  serv_lvl_arr <- serv_lvl_arr[order(serv_lvl_arr$station_id_index, serv_lvl_arr$tw_group),]
  return(serv_lvl_arr)
}

aggdata <- function(wdcMerged) {  
  #this function aggregates the data by the state of the system in terms of stock outs
  #wdcMerged <- wdcMerged[ fforder(wdcMerged$n_time_int, wdcMerged$station_id),]  

  #aggregating across days for now so now use of day_index,
  #but setting it to 1 for rest of the code which is designed to have different day 
  #fixed effect - delta. this could be like a month or tw later where different fixed
  #effects are there.
  #wdcMerged$day_index <- ffrep.int(1,nrow(wdcMerged))
  #average demand at a station for a state  
  agg = binned_sum(wdcMerged$out_dem,wdcMerged$stid_twg_locstate_fac)  
  agg_in_dem = binned_sum(wdcMerged$in_dem,wdcMerged$stid_twg_locstate_fac)
  if(!identical(rownames(agg),rownames(agg_in_dem))) stop("in aggdata, aggregated data not in same order")
  wdcMerged_sum <- as.ffdf(data.frame(stid_twg_locstate_fac=rownames(agg),
      out_dem_sum=agg[,"sum"],in_dem_sum=agg_in_dem[,"sum"],obs_weight=agg[,"count"]))
  
  #keep only one obs per station, sto_state and merge with wdcMerged_sum using stid_sto_state_fac  
  wdcMerged$dup <- duplicated(wdcMerged$stid_twg_locstate_fac)
  wdcMerged <- subset.ffdf(wdcMerged,dup == "FALSE")
  wdcMerged$dup <- NULL
  if(nrow(wdcMerged)!=nrow(wdcMerged_sum)) stop("no obs dont match in aggdata")
  wdcMerged =merge(wdcMerged,wdcMerged_sum,by="stid_twg_locstate_fac")
  'wdcMerged <- wdcMerged[,c("station_id", "stocked_out", 
      "lat", "lon", "station_id_index", "n_time_int", "stid_twg_locstate_fac",
      "out_dem_sum","out_dem_sum_scaled","obs_weight","serv_lvl")]
  '
  wdcMerged <- wdcMerged[fforder(wdcMerged$station_id_index,wdcMerged$stid_twg_locstate_fac),]

  return(as.data.frame(wdcMerged))
}


readtractboundaryfile_raw <- function() {
  #filename <- file.path(csv_dir,"tract_st_mod_parsed.csv", fsep = .Platform$file.sep)
  dir <- paste0(dropbox_dir,"/../VelibData/ParisData/ParisShapeFiles") 
  filename <- file.path(dir,"tract_st_mod_parsed.csv", fsep = .Platform$file.sep)
  wdcCenTrParsed <- read.csv(file=filename,head=FALSE,sep=",")
  colnames(wdcCenTrParsed) <- c("tract","lat","lon")
  return(wdcCenTrParsed)
}

readtractdensityfile <- function() {
  density_filename <- paste0(dropbox_dir,"/../VelibData/ParisData/CensusData/Insee/DistrictDensityData.csv") 
  wdcCenTrDensity <- read.csv(file=density_filename,head=TRUE,sep=",")
  colnames(wdcCenTrDensity) <- c("tract","density")
  return(wdcCenTrDensity)
}


assign_points_tract <- function(lat1,lon1,wdcCenTrParsed) {
  lat2 <- wdcCenTrParsed[,"lat"]
  lon2 <- wdcCenTrParsed[,"lon"]
  tract <- rep(NA,length(lat1))
  for(i in 1:length(lat1)) {
    dis_v <- latlondistance(lat1[i], lon1[i], lat2, lon2)
    tract[i] <- wdcCenTrParsed[which.min(dis_v),"tract"]  
  }
  return(tract)
}

assign_stid_cluster_tract <- function(station_id_list, st_cluster) {
  station_id_list_df <- data.frame(station_id=station_id_list)
  tract_df <- merge(station_id_list_df, st_cluster, by="station_id")
  tract_df <- tract_df[order(tract_df$station_id_index),]
  return(tract_df$cluster_id)
}

assign_points_cluster_tract <- function(lat1, lon1, st_cluster) {
  lat2 <- st_cluster[,"lat"]
  lon2 <- st_cluster[,"lon"]
  tract <- rep(NA,length(lat1))
  for(i in 1:length(lat1)) {
    dis_v <- latlondistance(lat1[i], lon1[i], lat2, lon2)
    tract[i] <- st_cluster[which.min(dis_v),"cluster_id"]  
  }
  return(tract)
}


second_stage_reg_data <- function(thetaout,serv_lvl_arr) {
  serv_lvl_arr <- serv_lvl_arr[order(serv_lvl_arr$tw_group,serv_lvl_arr$station_id_index),]
  #generate instruments for serv_lvl
  #use average service level of stations between 400 mts and 800 mts
  latlontemp <- unique(wdcMerged[,c("station_id_index","lat","lon")])
  #merge latlontemp with serv_lvl_arr
  serv_lvl_arr_merged <- merge(serv_lvl_arr, latlontemp, by ="station_id_index"
                               , all.x=TRUE, all.y=FALSE)
  
  instr_serv_lvl <- rep(NA,nrow(serv_lvl_arr_merged))
  instr_st_count <- rep(NA,nrow(serv_lvl_arr_merged))
  instr_dis_low = 0.4
  instr_dis_high = 0.8
  for(i in 1:nrow(serv_lvl_arr_merged)) {
    lat1 = serv_lvl_arr_merged[i,"lat"]
    lon1 = serv_lvl_arr_merged[i,"lon"]
    lat2 = serv_lvl_arr_merged$lat
    lon2 = serv_lvl_arr_merged$lon
    tw_groupin = serv_lvl_arr_merged$tw_group[i]
    dis_v <- latlondistance(lat1,lon1,lat2,lon2)
    instr_serv_lvl[i] = mean(serv_lvl_arr_merged$serv_lvl[which(dis_v>=instr_dis_low &
                            dis_v<=instr_dis_high & serv_lvl_arr_merged$tw_group==tw_groupin)])
    instr_st_count[i] = length(which(dis_v>=instr_dis_low &
                            dis_v<=instr_dis_high & serv_lvl_arr_merged$tw_group==tw_groupin))
  }
  
  instr_serv_lvl[is.nan(instr_serv_lvl)]=0
  serv_lvl_arr_merged$instr_serv_lvl <- instr_serv_lvl
  serv_lvl_arr_merged <- serv_lvl_arr_merged[order(serv_lvl_arr_merged$tw_group,serv_lvl_arr_merged$station_id_index),]
  cor(serv_lvl_arr_merged$instr_serv_lvl, serv_lvl_arr_merged$serv_lvl)
  return(serv_lvl_arr_merged)
}

predic_instr_serv_lvl <- function() {
  #generate instruments for serv_lvl
  #use average service level of stations between 400 mts and 800 mts
  latlontemp <- unique(wdcMerged[,c("station_id_index","lat","lon")])
  latlontemp <- latlontemp[order(latlontemp$station_id_index),]
  no_st <- nrow(latlontemp)
  tw_group_list <- unique(wdcMerged$tw_group)
    
  instr_serv_lvl <- rep(NA,no_st)
  instr_st_count <- rep(NA,no_st)
  instr_dis_low = 0.4
  instr_dis_high = 0.8
  instr_serv_lvl <- c()
  for(tw_in in tw_group_list) {
    for(i in 1:no_st) {
      lat1 = latlontemp[i,"lat"]
      lon1 = latlontemp[i,"lon"]
      lat2 = latlontemp$lat
      lon2 = latlontemp$lon
      dis_v <- latlondistance(lat1,lon1,lat2,lon2)
      instr_serv_lvl_i = mean(serv_lvl_arr$serv_lvl[which(dis_v>=instr_dis_low &
                              dis_v<=instr_dis_high & serv_lvl_arr$tw_group==tw_in)])
      instr_st_count_i = length(which(dis_v>=instr_dis_low &
                            dis_v<=instr_dis_high & serv_lvl_arr$tw_group==tw_in))
      instr_serv_lvl <- rbind(instr_serv_lvl,c(i,tw_in,instr_serv_lvl_i,instr_st_count_i))
    }
  }  
  colnames(instr_serv_lvl) <- c("station_id_index","tw_group","instr_serv_lvl",
                                "instr_st_count")
  instr_serv_lvl <- as.data.frame(instr_serv_lvl)
  instr_serv_lvl$instr_serv_lvl[is.nan(instr_serv_lvl$instr_serv_lvl)]=0
  serv_lvl_arr <- serv_lvl_arr[order(serv_lvl_arr$station_id_index,
                                     serv_lvl_arr$tw_group),]
  instr_serv_lvl <- instr_serv_lvl[order(instr_serv_lvl$station_id_index,
                                         instr_serv_lvl$tw_group),]
  serv_lvl_arr$instr_serv_lvl <- instr_serv_lvl$instr_serv_lvl
  fit <- lm(serv_lvl ~ instr_serv_lvl,serv_lvl_arr)
  predicted_serv_lvl <- coef(fit)[1] + coef(fit)[2]*serv_lvl_arr$serv_lvl
  return(data.frame("predicted_serv_lvl"=predicted_serv_lvl,"station_id_index"
      =serv_lvl_arr$station_id_index,"tw_group"=serv_lvl_arr$tw_group,
      "instr_serv_lvl"=serv_lvl_arr$instr_serv_lvl))  
}

instr_serv_lvl <- function(serv_lvl_arr,wdcMerged) {
  #generate instruments for serv_lvl
  #use average service level of stations between 400 mts and 800 mts
  latlontemp <- unique(wdcMerged[,c("station_id_index","lat","lon")])
  latlontemp <- latlontemp[order(latlontemp$station_id_index),]
  no_st <- nrow(latlontemp)
  tw_group_list <- unique(wdcMerged$tw_group)
  
  instr_serv_lvl <- rep(NA,no_st)
  instr_st_count <- rep(NA,no_st)
  instr_dis_low = 0.4
  instr_dis_high = 0.8
  instr_serv_lvl <- c()
  for(tw_in in tw_group_list) {
    for(i in 1:no_st) {
      lat1 = latlontemp[i,"lat"]
      lon1 = latlontemp[i,"lon"]
      lat2 = latlontemp$lat
      lon2 = latlontemp$lon
      dis_v <- latlondistance(lat1,lon1,lat2,lon2)
      instr_serv_lvl_i = mean(serv_lvl_arr$serv_lvl[which(dis_v>=instr_dis_low &
                                                            dis_v<=instr_dis_high & serv_lvl_arr$tw_group==tw_in)])
      instr_st_count_i = length(which(dis_v>=instr_dis_low &
                                        dis_v<=instr_dis_high & serv_lvl_arr$tw_group==tw_in))
      instr_serv_lvl <- rbind(instr_serv_lvl,c(i,tw_in,instr_serv_lvl_i,instr_st_count_i))
    }
  }  
  colnames(instr_serv_lvl) <- c("station_id_index","tw_group","instr_serv_lvl",
                                "instr_st_count")
  instr_serv_lvl <- as.data.frame(instr_serv_lvl)
  instr_serv_lvl$instr_serv_lvl[is.nan(instr_serv_lvl$instr_serv_lvl)]=0
  serv_lvl_arr <- serv_lvl_arr[order(serv_lvl_arr$station_id_index,
                                     serv_lvl_arr$tw_group),]
  instr_serv_lvl <- instr_serv_lvl[order(instr_serv_lvl$station_id_index,
                                         instr_serv_lvl$tw_group),]
  return(instr_serv_lvl$instr_serv_lvl)
}


instr_serv_lvl_new <- function(serv_lvl_arr,wdcMerged) {
  #generate instruments for serv_lvl
  #use average service level of stations between 400 mts and 800 mts
  latlontemp <- unique(wdcMerged[,c("station_id_index","lat","lon")])
  latlontemp <- latlontemp[order(latlontemp$station_id_index),]
  no_st <- nrow(latlontemp)
  tw_group_list <- unique(wdcMerged$tw_group)
  
  #   instr_serv_lvl <- rep(NA,no_st)
  #   instr_st_count <- rep(NA,no_st)
  instr_dis_low = max_walking_dis
  instr_dis_high = 2*max_walking_dis
  instr_serv_lvl <- c()
  for(tw_in in tw_group_list) {
    for(i in 1:no_st) {
      lat1 = latlontemp[i,"lat"]
      lon1 = latlontemp[i,"lon"]
      lat2 = latlontemp$lat
      lon2 = latlontemp$lon
      dis_v <- latlondistance(lat1,lon1,lat2,lon2)
      instr_serv_lvl_i = mean(serv_lvl_arr$serv_lvl[which(dis_v>=instr_dis_low &
                                                            dis_v<=instr_dis_high & serv_lvl_arr$tw_group==tw_in)])
      instr_st_count_i = length(which(dis_v>=instr_dis_low &
                                        dis_v<=instr_dis_high & serv_lvl_arr$tw_group==tw_in))
      instr_serv_lvl <- rbind(instr_serv_lvl,c(i,tw_in,instr_serv_lvl_i))
    }
  }  
  colnames(instr_serv_lvl) <- c("station_id_index","tw_group","instr_serv_lvl")
  instr_serv_lvl <- as.data.frame(instr_serv_lvl)
  instr_serv_lvl$instr_serv_lvl[is.nan(instr_serv_lvl$instr_serv_lvl)]=0
  if(!identical(order(serv_lvl_arr$station_id_index,serv_lvl_arr$tw_group),
                order(instr_serv_lvl$station_id_index,instr_serv_lvl$tw_group))) stop("not in order")
  return(instr_serv_lvl$instr_serv_lvl)
}


reg_run <- function(wdcMergedin=NULL) {
  wdcMerged <- reg_ready_data(wdcMergedin,filelist)
  #wdcMergedOrg <- wdcMerged
  tw_low=8
  tw_high=19
  wdcMerged <- subset(wdcMerged,tw >=tw_low & tw <=tw_high )
  
  #gen factors
  wdcMerged$day_fac <-  as.ff(as.character(wdcMerged$day))
  wdcMerged$tw_fac <-  as.ff(as.character(wdcMerged$tw))
  #time window group - level at which fixed effects are there  
  wdcMerged$tw_group <- wdcMerged$tw
  wdcMerged$tw_group_fac <-  as.ff(as.character(wdcMerged$tw_group))
  
  wdcMerged$stid_twg_locstate_fac <-  as.ff(as.factor(paste(wdcMerged$station_id_index[],
                                                            wdcMerged$tw_group[],wdcMerged$sto_state_local[]))) 
  
  serv_lvl_arr <- gen_serv_lvls (wdcMerged,tw_low,tw_high)
  
  wdcMerged$sto_state_local_fac <- as.ff(as.factor(wdcMerged$sto_state_local[]))
  #no_st <- length(st_list)
  wdcMerged <- subset(wdcMerged,stocked_out==FALSE )
  wdcMerged <- droplevels(wdcMerged)
  
  #wdcMerged$stid_sto_state_local_fac <-  as.ff(as.factor(paste(wdcMerged$station_id_index[],wdcMerged$sto_state_local[]))) 
    
  agg <- binned_sum(wdcMerged$stid_twg_locstate_fac,wdcMerged$stid_twg_locstate_fac)
  if(length(which(agg[,"count"]==0))>0) stop("levels incorrect agg")
  agg_st_id <- binned_sum(wdcMerged$station_id_index,wdcMerged$stid_twg_locstate_fac)
  if(length(which(agg_st_id[,"count"]==0))>0) stop("levels incorrect agg_st_id")
  agg_st_id <- agg_st_id[,"sum"]/agg_st_id[,"count"]
  agg_twg <- binned_sum(wdcMerged$tw_group,wdcMerged$stid_twg_locstate_fac)
  if(length(which(agg_twg[,"count"]==0))>0) stop("levels incorrect agg_twg")
  agg_twg <- agg_twg[,"sum"]/agg_twg[,"count"]
    
  agg <- ffdf(stid_twg_locstate_fac=as.ff(as.factor(row.names(agg))), 
              stid_twg_locstate_obs_count=as.ff(agg[,"count"]), 
               station_id_index= as.ff(agg_st_id), tw_group= as.ff(agg_twg) ,row.names=NULL) 
  
  agg$groupbyfactor <- as.ff(as.factor(paste(agg$station_id_index[], agg$tw_group[]))) 
  agg <- agg[fforder(agg$groupbyfactor,agg$stid_twg_locstate_obs_count),]
  #select top 10 states
  
  revIndexFUN <- function(x){
    x$groupbyfactor <- droplevels(x$groupbyfactor)  
    x <- split(x, x$groupbyfactor)
    x <- lapply(x, FUN=function(onlyonegroup){      
      ret3 <- c(nrow(onlyonegroup):1)
      data.frame(index=ret3,onlyonegroup)
    })
    x <- do.call(rbind, x)      
    x
  }
    
  agg <- as.data.frame(agg)
  index <- ave(agg$station_id_index, agg$groupbyfactor, FUN=function(x) c(length(x):1))

#   index <- ffdfdply(agg[c("station_id_index","groupbyfactor","stid_twg_locstate_obs_count","stid_twg_locstate_fac")],
#                     agg$groupbyfactor,FUN=function(x) revIndexFUN(x) )
#   if(!identical(index$stid_twg_locstate_fac[],
#                 agg$stid_twg_locstate_fac[])) stop("not in order")
#     
#   agg <- agg[ffwhich(index,index$index<=10),]

  agg <- agg[which(index<=10),]
  agg$station_id_index <- NULL
  agg$tw_group <- NULL
  agg <- as.ffdf(agg)
  wdcMerged <- merge(wdcMerged,agg,by="stid_twg_locstate_fac")
  
  
  
#  list_all <- ffwhich(wdcMerged,tw >=tw_low & tw <=tw_high & sto_state_fac %in% top_states)
#   list_all <- ffwhich(wdcMerged,tw >=tw_low & tw <=tw_high)
#   wdcMerged <- wdcMerged[list_all,]
  
  list <- ffwhich(wdcMerged, stocked_out==FALSE )

  wdcMerged_sub <- wdcMerged[list,]
  #droplevels for entire ffdf tends to exchange some columns, leading to quite weird behaviour 
  wdcMerged_sub$tw_fac <- droplevels(wdcMerged_sub$tw_fac)
  wdcMerged_sub$day_fac <- droplevels(wdcMerged_sub$day_fac)
  wdcMerged_sub$station_id_index_fac <-  as.ff(as.character(wdcMerged_sub$station_id_index))
  
  model="poisson"
  system.time({
    fit <- bigglm(out_dem ~   tw_fac + day_fac +  station_id_index_fac, #   
                  data = wdcMerged_sub, family = poisson(link = "log"),
                  chunksize= 5000, maxit=100)
  })
  #remove the day and time effect
  coef_fit <- coef(fit)
  #tw_predic <- as.numeric(coef_fit[paste0('tw_fac',wdcMerged$tw_fac[list][])])
  day_predic <- as.numeric(coef_fit[paste0('day_fac',wdcMerged$day_fac[list][])])
  #replace NA's by 0 as they correpond to the base factors whose coefficients are 0.
  #tw_predic[which(is.na(tw_predic))] <- 0
  day_predic[which(is.na(day_predic))] <- 0
  day_predic <- day_predic - mean(day_predic)
  fitted <- ffrep.int(0,nrow(wdcMerged))  
  #fitted[list] <-  as.ff(tw_predic + day_predic)
  fitted[list] <-  as.ff(day_predic)
  if(model=="poisson") {
    fitted <- exp(fitted)
  }
  #wdcMerged$out_dem_save <- wdcMerged$out_dem
  wdcMerged[,"out_dem"] <-  wdcMerged[,"out_dem"]/fitted[]

  wdcMerged$stid_twg_locstate_fac <- droplevels(wdcMerged$stid_twg_locstate_fac)
  wdcMerged <- aggdata(wdcMerged)
  return(list("wdcMerged"=wdcMerged,"serv_lvl_arr"=serv_lvl_arr))  
}

# reg_run_serv <- function(wdcMergedin=NULL) {
#   wdcMerged <- reg_ready_data(wdcMergedin)
#   #wdcMergedOrg <- wdcMerged
#   tw_low=8
#   tw_high=19
#   wdcMerged <- subset(wdcMerged,tw >=tw_low & tw <=tw_high )
#   #wdcMerged$sto_state_fac <- as.ff(as.factor(wdcMerged$sto_state[]))
#   no_st <- length(st_list)
#   
#   agg <- binned_sum(wdcMerged$sto_state_fac,wdcMerged$sto_state_fac)
#   agg <- data.frame(sto_state_fac=row.names(agg), sto_state_obs_count=agg[,"count"]/no_st, row.names=NULL) 
#   agg <- agg[order(agg$sto_state_obs_count,decreasing=TRUE),]
#   #select top 100 states
#   top_states <- agg$sto_state_fac[1:min(40,nrow(agg))]
#   
#   wdcMerged <- merge(wdcMerged,as.ffdf(agg),by="sto_state_fac")
#   
#   serv_lvl_arr <- gen_serv_lvls (wdcMerged,tw_low,tw_high)
#   return(serv_lvl_arr)
# }


compute_serv_level <- function(filelistin) {
  wdcMerged <- reg_ready_data(filelist=filelistin)
  #wdcMergedOrg <- wdcMerged
   tw_low=0
   tw_high=24
#   wdcMerged <- subset(wdcMerged,tw >=tw_low & tw <=tw_high )
  
  #gen factors
  wdcMerged$day_fac <-  as.ff(as.character(wdcMerged$day))
  wdcMerged$tw_fac <-  as.ff(as.character(wdcMerged$tw))
  #time window group - level at which fixed effects are there  
  wdcMerged$tw_group <- wdcMerged$tw
  wdcMerged$tw_group_fac <-  as.ff(as.character(wdcMerged$tw_group))
  
  wdcMerged$stid_twg_locstate_fac <-  as.ff(as.factor(paste(wdcMerged$station_id_index[],
                                                            wdcMerged$tw_group[],wdcMerged$sto_state_local[]))) 
  
  serv_lvl_arr <- gen_serv_lvls (wdcMerged,tw_low,tw_high)
  return(serv_lvl_arr)
}




reg_run_bootstrap <- function(wdcMergedin=NULL) {
  save_bootstrap_file <- paste0(csv_dir,"/ffdb/wdcMerged_bootstrap/tract_",paste0(tract,collapse="_"))
  if(file.exists(save_bootstrap_file)) {
    load.ffdf(dir=save_bootstrap_file)
    return(wdcMerged)
  }
  
  wdcMerged <- reg_ready_data(wdcMergedin,filelist)
  #wdcMergedOrg <- wdcMerged
   tw_low=0
   tw_high=24
#   wdcMerged <- subset(wdcMerged,tw >=tw_low & tw <=tw_high )
  
  #gen factors
  wdcMerged$day_fac <-  as.ff(as.character(wdcMerged$day))
  wdcMerged$tw_fac <-  as.ff(as.character(wdcMerged$tw))
  #time window group - level at which fixed effects are there  
  wdcMerged$tw_group <- wdcMerged$tw
  wdcMerged$tw_group_fac <-  as.ff(as.character(wdcMerged$tw_group))
  
  wdcMerged$stid_twg_locstate_fac <-  as.ff(as.factor(paste(wdcMerged$station_id_index[],
                                                            wdcMerged$tw_group[],wdcMerged$sto_state_local[]))) 
  
  #serv_lvl_arr <- gen_serv_lvls (wdcMerged,tw_low,tw_high)
  
  wdcMerged$sto_state_local_fac <- as.ff(as.factor(wdcMerged$sto_state_local[]))
  #no_st <- length(st_list)
  wdcMerged <- subset(wdcMerged,stocked_out==FALSE )
  wdcMerged <- droplevels(wdcMerged)
  
  #wdcMerged$stid_sto_state_local_fac <-  as.ff(as.factor(paste(wdcMerged$station_id_index[],wdcMerged$sto_state_local[]))) 
  #selecting top 30 states
#   agg <- binned_sum(wdcMerged$stid_twg_locstate_fac,wdcMerged$stid_twg_locstate_fac)
#   if(length(which(agg[,"count"]==0))>0) stop("levels incorrect agg")
#   agg_st_id <- binned_sum(wdcMerged$station_id_index,wdcMerged$stid_twg_locstate_fac)
#   if(length(which(agg_st_id[,"count"]==0))>0) stop("levels incorrect agg_st_id")
#   agg_st_id <- agg_st_id[,"sum"]/agg_st_id[,"count"]
#   agg_twg <- binned_sum(wdcMerged$tw_group,wdcMerged$stid_twg_locstate_fac)
#   if(length(which(agg_twg[,"count"]==0))>0) stop("levels incorrect agg_twg")
#   agg_twg <- agg_twg[,"sum"]/agg_twg[,"count"]
#   
#   agg <- ffdf(stid_twg_locstate_fac=as.ff(as.factor(row.names(agg))), 
#               stid_twg_locstate_obs_count=as.ff(agg[,"count"]), 
#               station_id_index= as.ff(agg_st_id), tw_group= as.ff(agg_twg) ,row.names=NULL) 
#   
#   agg$groupbyfactor <- as.ff(as.factor(paste(agg$station_id_index[], agg$tw_group[]))) 
#   agg <- agg[fforder(agg$groupbyfactor,agg$stid_twg_locstate_obs_count),]
#   #select top 10 states
#   
#   revIndexFUN <- function(x){
#     x$groupbyfactor <- droplevels(x$groupbyfactor)  
#     x <- split(x, x$groupbyfactor)
#     x <- lapply(x, FUN=function(onlyonegroup){      
#       ret3 <- c(nrow(onlyonegroup):1)
#       data.frame(index=ret3,onlyonegroup)
#     })
#     x <- do.call(rbind, x)      
#     x
#   }
#   
#   agg <- as.data.frame(agg)
#   index <- ave(agg$station_id_index, agg$groupbyfactor, FUN=function(x) c(length(x):1))
#   
#   #   index <- ffdfdply(agg[c("station_id_index","groupbyfactor","stid_twg_locstate_obs_count","stid_twg_locstate_fac")],
#   #                     agg$groupbyfactor,FUN=function(x) revIndexFUN(x) )
#   #   if(!identical(index$stid_twg_locstate_fac[],
#   #                 agg$stid_twg_locstate_fac[])) stop("not in order")
#   #     
#   #   agg <- agg[ffwhich(index,index$index<=10),]
#   
#   agg <- agg[which(index<=30),]
#   agg$station_id_index <- NULL
#   agg$tw_group <- NULL
#   agg <- as.ffdf(agg)
#   wdcMerged <- merge(wdcMerged,agg,by="stid_twg_locstate_fac")
  
  
  
  #  list_all <- ffwhich(wdcMerged,tw >=tw_low & tw <=tw_high & sto_state_fac %in% top_states)
  #   list_all <- ffwhich(wdcMerged,tw >=tw_low & tw <=tw_high)
  #   wdcMerged <- wdcMerged[list_all,]
  
  list <- ffwhich(wdcMerged, stocked_out==FALSE )
  
  wdcMerged_sub <- wdcMerged[list,]
  #droplevels for entire ffdf tends to exchange some columns, leading to quite weird behaviour 
  wdcMerged_sub$tw_fac <- droplevels(wdcMerged_sub$tw_fac)
  wdcMerged_sub$day_fac <- droplevels(wdcMerged_sub$day_fac)
  wdcMerged_sub$station_id_index_fac <-  as.ff(as.character(wdcMerged_sub$station_id_index))
  
  model="poisson"
  system.time({
    fit <- bigglm(out_dem ~   tw_fac + day_fac +  station_id_index_fac, #   
                  data = wdcMerged_sub, family = poisson(link = "log"),
                  chunksize= 5000, maxit=100)
  })
  #remove the day and time effect
  coef_fit <- coef(fit)
  #tw_predic <- as.numeric(coef_fit[paste0('tw_fac',wdcMerged$tw_fac[list][])])
  day_predic <- as.numeric(coef_fit[paste0('day_fac',wdcMerged$day_fac[list][])])
  #replace NA's by 0 as they correpond to the base factors whose coefficients are 0.
  #tw_predic[which(is.na(tw_predic))] <- 0
  day_predic[which(is.na(day_predic))] <- 0
  day_predic <- day_predic - mean(day_predic)
  fitted <- ffrep.int(0,nrow(wdcMerged))  
  #fitted[list] <-  as.ff(tw_predic + day_predic)
  fitted[list] <-  as.ff(day_predic)
  if(model=="poisson") {
    fitted <- exp(fitted)
  }
  #wdcMerged$out_dem_save <- wdcMerged$out_dem
  wdcMerged[,"out_dem"] <-  wdcMerged[,"out_dem"]/fitted[]
  
  wdcMerged$stid_twg_locstate_fac <- droplevels(wdcMerged$stid_twg_locstate_fac)
  
  save.ffdf(wdcMerged,dir=save_bootstrap_file)
  return(wdcMerged)
}



readtractboundaryfile <- function() {
  
  tractboundaryfile <- readtractboundaryfile_raw()
  tract_list <- unique(tractboundaryfile$tract)
  tractboundaryfile_expand <- data.frame()
  for(tract in tract_list) {
    tractboundaryfile_sub <- tractboundaryfile[which(tractboundaryfile$tract==tract),]
    N = nrow(tractboundaryfile_sub)
    dis = latlondistance(tractboundaryfile_sub$lat,tractboundaryfile_sub$lon,
                         tractboundaryfile_sub$lat[c(2:N,1)],tractboundaryfile_sub$lon[c(2:N,1)])
    points_dis = 0.050
    no_points <- ceiling(dis/points_dis)
    
    tractboundaryfile_sub$next_lat <- tractboundaryfile_sub$lat[c(2:N,1)]
    tractboundaryfile_sub$next_lon <- tractboundaryfile_sub$lon[c(2:N,1)]
    tractboundaryfile_sub$no_points <- no_points
    
    interpolate_points <- function(x) {
      lat1 = x[1]
      lon1 = x[2]
      lat2 = x[3]
      lon2 = x[4]
      np = x[5]
      weights = 0:(np)/(np+1)
      lat = lat1+(weights*(lat2-lat1))
      lon = lon1+(weights*(lon2-lon1))
      return(cbind(lat,lon))
    }
    a=tractboundaryfile_sub[,c("lat","lon","next_lat","next_lon","no_points")]
    ret = apply(a,1,FUN= function(x)interpolate_points(x))
    ret = as.data.frame(do.call(rbind,ret),row.names=NA)
    ret$tract <- tract
    #shrink ret
    median_lat = median(ret$lat)
    median_lon = median(ret$lon)
    ret$lat = median_lat + (ret$lat-median_lat)*0.95
    ret$lon = median_lon + (ret$lon-median_lon)*0.95
    tractboundaryfile_expand <- rbind(tractboundaryfile_expand,ret)  
  }
  return(tractboundaryfile_expand)
}  

reg_run_bootstrap_2 <- function(wdcMergedin=NULL,save_bootstrap_file=NULL) {
  #just changes the saving location
  if(is.null(save_bootstrap_file)) {
    save_bootstrap_file <- paste0(csv_dir,"/ffdb/wdcMerged_bootstrap_2/tract_",paste0(tract,collapse="_"))
  }
  if(file.exists(save_bootstrap_file)) {
    load.ffdf(dir=save_bootstrap_file)
    return(wdcMerged)
  }
  
  wdcMerged <- reg_ready_data(wdcMergedin,filelist)
  #wdcMergedOrg <- wdcMerged
  tw_low=0
  tw_high=24
  #   wdcMerged <- subset(wdcMerged,tw >=tw_low & tw <=tw_high )
  
  #gen factors
  wdcMerged$day_fac <-  as.ff(as.character(wdcMerged$day))
  wdcMerged$tw_fac <-  as.ff(as.character(wdcMerged$tw))
  #time window group - level at which fixed effects are there  
  wdcMerged$tw_group <- wdcMerged$tw
  wdcMerged$tw_group_fac <-  as.ff(as.character(wdcMerged$tw_group))
  
  wdcMerged$stid_twg_locstate_fac <-  as.ff(as.factor(paste(wdcMerged$station_id_index[],
                                                            wdcMerged$tw_group[],wdcMerged$sto_state_local[]))) 
  
  #serv_lvl_arr <- gen_serv_lvls (wdcMerged,tw_low,tw_high)
  
  wdcMerged$sto_state_local_fac <- as.ff(as.factor(wdcMerged$sto_state_local[]))
  #no_st <- length(st_list)
  wdcMerged <- subset(wdcMerged,stocked_out==FALSE )
  wdcMerged <- droplevels(wdcMerged)
  
  #wdcMerged$stid_sto_state_local_fac <-  as.ff(as.factor(paste(wdcMerged$station_id_index[],wdcMerged$sto_state_local[]))) 
  #selecting top 30 states
  #   agg <- binned_sum(wdcMerged$stid_twg_locstate_fac,wdcMerged$stid_twg_locstate_fac)
  #   if(length(which(agg[,"count"]==0))>0) stop("levels incorrect agg")
  #   agg_st_id <- binned_sum(wdcMerged$station_id_index,wdcMerged$stid_twg_locstate_fac)
  #   if(length(which(agg_st_id[,"count"]==0))>0) stop("levels incorrect agg_st_id")
  #   agg_st_id <- agg_st_id[,"sum"]/agg_st_id[,"count"]
  #   agg_twg <- binned_sum(wdcMerged$tw_group,wdcMerged$stid_twg_locstate_fac)
  #   if(length(which(agg_twg[,"count"]==0))>0) stop("levels incorrect agg_twg")
  #   agg_twg <- agg_twg[,"sum"]/agg_twg[,"count"]
  #   
  #   agg <- ffdf(stid_twg_locstate_fac=as.ff(as.factor(row.names(agg))), 
  #               stid_twg_locstate_obs_count=as.ff(agg[,"count"]), 
  #               station_id_index= as.ff(agg_st_id), tw_group= as.ff(agg_twg) ,row.names=NULL) 
  #   
  #   agg$groupbyfactor <- as.ff(as.factor(paste(agg$station_id_index[], agg$tw_group[]))) 
  #   agg <- agg[fforder(agg$groupbyfactor,agg$stid_twg_locstate_obs_count),]
  #   #select top 10 states
  #   
  #   revIndexFUN <- function(x){
  #     x$groupbyfactor <- droplevels(x$groupbyfactor)  
  #     x <- split(x, x$groupbyfactor)
  #     x <- lapply(x, FUN=function(onlyonegroup){      
  #       ret3 <- c(nrow(onlyonegroup):1)
  #       data.frame(index=ret3,onlyonegroup)
  #     })
  #     x <- do.call(rbind, x)      
  #     x
  #   }
  #   
  #   agg <- as.data.frame(agg)
  #   index <- ave(agg$station_id_index, agg$groupbyfactor, FUN=function(x) c(length(x):1))
  #   
  #   #   index <- ffdfdply(agg[c("station_id_index","groupbyfactor","stid_twg_locstate_obs_count","stid_twg_locstate_fac")],
  #   #                     agg$groupbyfactor,FUN=function(x) revIndexFUN(x) )
  #   #   if(!identical(index$stid_twg_locstate_fac[],
  #   #                 agg$stid_twg_locstate_fac[])) stop("not in order")
  #   #     
  #   #   agg <- agg[ffwhich(index,index$index<=10),]
  #   
  #   agg <- agg[which(index<=30),]
  #   agg$station_id_index <- NULL
  #   agg$tw_group <- NULL
  #   agg <- as.ffdf(agg)
  #   wdcMerged <- merge(wdcMerged,agg,by="stid_twg_locstate_fac")
  
  
  
  #  list_all <- ffwhich(wdcMerged,tw >=tw_low & tw <=tw_high & sto_state_fac %in% top_states)
  #   list_all <- ffwhich(wdcMerged,tw >=tw_low & tw <=tw_high)
  #   wdcMerged <- wdcMerged[list_all,]
  
  list <- ffwhich(wdcMerged, stocked_out==FALSE )
  
  wdcMerged_sub <- wdcMerged[list,]
  #droplevels for entire ffdf tends to exchange some columns, leading to quite weird behaviour 
  wdcMerged_sub$tw_fac <- droplevels(wdcMerged_sub$tw_fac)
  wdcMerged_sub$day_fac <- droplevels(wdcMerged_sub$day_fac)
  wdcMerged_sub$station_id_index_fac <-  as.ff(as.character(wdcMerged_sub$station_id_index))
  
  model="poisson"
  system.time({
    fit <- bigglm(out_dem ~   tw_fac + day_fac +  station_id_index_fac, #   
                  data = wdcMerged_sub, family = poisson(link = "log"),
                  chunksize= 5000, maxit=100)
  })
  #remove the day and time effect
  coef_fit <- coef(fit)
  #tw_predic <- as.numeric(coef_fit[paste0('tw_fac',wdcMerged$tw_fac[list][])])
  day_predic <- as.numeric(coef_fit[paste0('day_fac',wdcMerged$day_fac[list][])])
  #replace NA's by 0 as they correpond to the base factors whose coefficients are 0.
  #tw_predic[which(is.na(tw_predic))] <- 0
  day_predic[which(is.na(day_predic))] <- 0
  day_predic <- day_predic - mean(day_predic)
  fitted <- ffrep.int(0,nrow(wdcMerged))  
  #fitted[list] <-  as.ff(tw_predic + day_predic)
  fitted[list] <-  as.ff(day_predic)
  if(model=="poisson") {
    fitted <- exp(fitted)
  }
  #wdcMerged$out_dem_save <- wdcMerged$out_dem
  wdcMerged[,"out_dem"] <-  wdcMerged[,"out_dem"]/fitted[]
  
  wdcMerged$stid_twg_locstate_fac <- droplevels(wdcMerged$stid_twg_locstate_fac)
  
  save.ffdf(wdcMerged,dir=save_bootstrap_file)
  return(wdcMerged)
}

reg_run_bootstrap_simple <- function(wdcMergedin=NULL,save_bootstrap_file) {
  #just changes the saving location and doesnt do poissson regression.
  if(is.null(save_bootstrap_file)) {
    save_bootstrap_file <- paste0(csv_dir,"/ffdb/wdcMerged_bootstrap_2/tract_",paste0(tract,collapse="_"))
  }
  if(file.exists(save_bootstrap_file)) {
    load.ffdf(dir=save_bootstrap_file)
    return(wdcMerged)
  }
  wdcMerged <- reg_ready_data(wdcMergedin,filelist)  
  writefile <- paste0(csv_dir,"/filelogs/observationnumberlog.txt")
  write(paste("after reg_ready_data",nrow(wdcMerged)),file=writefile, append=T)
  write(paste("after reg_ready_data trips",sum(wdcMerged$out_dem)),file=writefile, append=T)
  
  #wdcMergedOrg <- wdcMerged
  tw_low=0
  tw_high=24
  #   wdcMerged <- subset(wdcMerged,tw >=tw_low & tw <=tw_high )
  
  #gen factors
  wdcMerged$day_fac <-  as.ff(as.character(wdcMerged$day))
  wdcMerged$tw_fac <-  as.ff(as.character(wdcMerged$tw))
  #time window group - level at which fixed effects are there  
  wdcMerged$tw_group <- wdcMerged$tw
  wdcMerged$tw_group_fac <-  as.ff(as.character(wdcMerged$tw_group))
  
  wdcMerged$stid_twg_locstate_fac <-  as.ff(as.factor(paste(wdcMerged$station_id_index[],
                                                            wdcMerged$tw_group[],wdcMerged$sto_state_local[]))) 
  
  #serv_lvl_arr <- gen_serv_lvls (wdcMerged,tw_low,tw_high)
  
  wdcMerged$sto_state_local_fac <- as.ff(as.factor(wdcMerged$sto_state_local[]))

  wdcMerged$stid_twg_locstate_fac <- droplevels(wdcMerged$stid_twg_locstate_fac)
  
  #remove data points that are likely reallocations
  writefile <- paste0(csv_dir,"/filelogs/observationnumberlog.txt")
  write(paste("before realloc removal",nrow(wdcMerged)),file=writefile, append=T)
  write(paste("before realloc removal trips",sum(wdcMerged$out_dem)),file=writefile, append=T)
  realloc_trunc = 4
  wdcMerged <- subset(wdcMerged, out_dem <=realloc_trunc)
  wdcMerged <- subset(wdcMerged, in_dem <=realloc_trunc)
  writefile <- paste0(csv_dir,"/filelogs/observationnumberlog.txt")
  write(paste("after realloc removal",nrow(wdcMerged)),file=writefile, append=T)
  write(paste("after realloc removal trips",sum(wdcMerged$out_dem)),file=writefile, append=T)
  save.ffdf(wdcMerged,dir=save_bootstrap_file)
  return(wdcMerged)
}

readtractboundaryfile_presentation <- function() {
  
 return(readtractboundaryfile_finegrain(0.050))
}  

readtractboundaryfile_finegrain <- function(granular_dis = 0.050 ) {
  
  tractboundaryfile <- readtractboundaryfile_raw()
  tract_list <- unique(tractboundaryfile$tract)
  tractboundaryfile_expand <- data.frame()
  for(tract in tract_list) {
    tractboundaryfile_sub <- tractboundaryfile[which(tractboundaryfile$tract==tract),]
    N = nrow(tractboundaryfile_sub)
    dis = latlondistance(tractboundaryfile_sub$lat,tractboundaryfile_sub$lon,
                         tractboundaryfile_sub$lat[c(2:N,1)],tractboundaryfile_sub$lon[c(2:N,1)])
    
    no_points <- ceiling(dis/granular_dis)
    
    tractboundaryfile_sub$next_lat <- tractboundaryfile_sub$lat[c(2:N,1)]
    tractboundaryfile_sub$next_lon <- tractboundaryfile_sub$lon[c(2:N,1)]
    tractboundaryfile_sub$no_points <- no_points
    
    interpolate_points <- function(x) {
      lat1 = x[1]
      lon1 = x[2]
      lat2 = x[3]
      lon2 = x[4]
      np = x[5]
      weights = 0:(np)/(np+1)
      lat = lat1+(weights*(lat2-lat1))
      lon = lon1+(weights*(lon2-lon1))
      return(cbind(lat,lon))
    }
    a=tractboundaryfile_sub[,c("lat","lon","next_lat","next_lon","no_points")]
    ret = apply(a,1,FUN= function(x)interpolate_points(x))
    ret = as.data.frame(do.call(rbind,ret),row.names=NA)
    ret$tract <- tract
    #shrink ret
    median_lat = median(ret$lat)
    median_lon = median(ret$lon)
    #     ret$lat = median_lat + (ret$lat-median_lat)*0.95
    #     ret$lon = median_lon + (ret$lon-median_lon)*0.95
    tractboundaryfile_expand <- rbind(tractboundaryfile_expand,ret)  
  }
  return(tractboundaryfile_expand)
}  

get_sl_alltw <- function(sl1_stations) {
  ret <- c()
  for(tw_in in tw_list) {
    ret <-  c(ret, mean(user_serv_lvl$serv_lvl[which(user_serv_lvl$station_id_index %in% sl1_stations & 
                                              user_serv_lvl$tw==tw_in)]))
  }
  return(ret)
}
