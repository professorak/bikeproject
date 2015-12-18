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
  #   realloc_trunc = 4
  #   wdcMerged <- subset(wdcMerged, out_dem <=realloc_trunc)
  #   wdcMerged <- subset(wdcMerged, in_dem <=realloc_trunc)
  writefile <- paste0(csv_dir,"/filelogs/observationnumberlog.txt")
  write(paste("after realloc removal",nrow(wdcMerged)),file=writefile, append=T)
  write(paste("after realloc removal trips",sum(wdcMerged$out_dem)),file=writefile, append=T)
  save.ffdf(wdcMerged,dir=save_bootstrap_file)
  return(wdcMerged)
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
  wdcRawOrg <- transform(wdcRawOrg,stocked_out=bikes<=stocked_out_threshold)
  wdcRawOrg <- transform(wdcRawOrg,stocked_full=spaces_available==0)
  
  #calculate outgoing demand
  wdcRawOrg <- transform(wdcRawOrg,
     out_dem = (floor(-change_bikes_avg + abs(change_bikes_avg))/2)*(stocked_out==F),
     in_dem = floor(change_bikes_avg + abs(change_bikes_avg))/2 )
  
  return(wdcRawOrg)
}

