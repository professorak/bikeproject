#compare bytw and averaged data


source("GetDataGooglePlaces_10dis_aggmonth_hyperlocalstate_norealloc_stkoutthresh_five_averaged_pr6_test.R")
wdcMerged_10dis_hyperlocal_pr6 <- wdcMerged

source("GetDataGooglePlaces_10dis_aggmonth_bytw_norealloc_stkoutthresh_five_averaged_pr6.R")
wdcMerged_10dis_bytw_pr6 <- wdcMerged


wdcMerged_10dis_hyperlocal_pr6_bytw_totdem <- 
  ave(wdcMerged_10dis_hyperlocal_pr6$out_dem_sum, wdcMerged_10dis_hyperlocal_pr6$station_id_index,
      wdcMerged_10dis_hyperlocal_pr6$tw_group, FUN = sum)

library(data.table)
wdcMerged_10dis_hyperlocal_pr6_dt <- data.table(wdcMerged_10dis_hyperlocal_pr6)
wdcMerged_10dis_hyperlocal_pr6_bytw_dt <- 
  wdcMerged_10dis_hyperlocal_pr6_dt[, list(tot_dem =sum(out_dem_sum), 
                                    tot_obs_weight =sum(obs_weight)), 
                                    by=list(station_id_index, tw_group)]

wdcMerged_10dis_hyperlocal_pr6_bytw_dt$dem_mean <- 
  wdcMerged_10dis_hyperlocal_pr6_bytw_dt$tot_dem/wdcMerged_10dis_hyperlocal_pr6_bytw_dt$tot_obs_weight

#check percentage differnce in out_dem_mean
wdcMerged_10dis_hyperlocal_pr6_bytw_dt <- wdcMerged_10dis_hyperlocal_pr6_bytw_dt[
  order(wdcMerged_10dis_hyperlocal_pr6_bytw_dt$tw_group,wdcMerged_10dis_hyperlocal_pr6_bytw_dt$station_id_index),]

wdcMerged_10dis_bytw_pr6 <- wdcMerged_10dis_bytw_pr6[
  order(wdcMerged_10dis_bytw_pr6$tw_group,wdcMerged_10dis_bytw_pr6$station_id_index),]

perc_dem_diff <- (wdcMerged_10dis_bytw_pr6$out_dem_mean - 
wdcMerged_10dis_hyperlocal_pr6_bytw_dt$dem_mean)/wdcMerged_10dis_bytw_pr6$out_dem_mean * 100

summary(perc_dem_diff)

length(which(abs(perc_dem_diff)>4))
erroneous_dt <- wdcMerged_10dis_hyperlocal_pr6_bytw_dt[which(abs(perc_dem_diff)>4),]
erroneous_dt <- cbind(erroneous_dt,
                      wdcMerged_10dis_bytw_pr6[which(abs(perc_dem_diff)>4),
                  c("out_dem_mean","out_dem_sum","obs_weight")])





