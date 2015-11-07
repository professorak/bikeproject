#compare bytw and averaged data

#test pr3 and pr6 demand number should be same
source("GetDataGooglePlaces_20dis_aggmonth_bytw_norealloc_stkoutthresh_five_averaged_pr3.R")
wdcMerged_20dis_bytw_pr3 <- wdcMerged

source("GetDataGooglePlaces_20dis_aggmonth_bytw_norealloc_stkoutthresh_five_averaged_pr6.R")
wdcMerged_20dis_bytw_pr6 <- wdcMerged

identical(wdcMerged_20dis_bytw_pr6$out_dem_mean, wdcMerged_20dis_bytw_pr3$out_dem_mean)

source("GetDataGooglePlaces_aggmonths_averaged_20districts_norealloc_stkoutthresh_five.R")
wdcMerged_20dis_avg_pr6 <- wdcMerged

source("GetDataGooglePlaces_aggmonths_averaged_20districts_norealloc_stkoutthresh_five_pr3.R")
wdcMerged_20dis_avg_pr3 <- wdcMerged

identical(wdcMerged_20dis_avg_pr6$out_dem_mean, wdcMerged_20dis_avg_pr3$out_dem_mean)

a <- as.numeric(by(wdcMerged_20dis_bytw_pr3$out_dem_mean, wdcMerged_20dis_bytw_pr3$station_id_index, FUN=mean))
plot(a)
summary(wdcMerged_20dis_avg_pr3$out_dem_mean-a)


source("GetDataGooglePlaces_aggmonths_20districts_norealloc_stkoutthresh_five_pr3.R")
wdcMerged_20dis_full_pr3 <- wdcMerged
#count the number of observations by each stationXtw
wdcMerged_20dis_full_pr3_obscount <- wdcMerged_20dis_full_pr3
wdcMerged_20dis_full_pr3_obscount$tot_obs <- ave(wdcMerged_20dis_full_pr3_obscount$obs_weight, )

source("GetDataGooglePlaces_aggmonths_20districts_norealloc_stkoutthresh_five.R")
wdcMerged_20dis_full_pr6 <- wdcMerged


