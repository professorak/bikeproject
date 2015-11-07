source("GetDataGooglePlaces_20dis_aggmonth_bytw_averaged_pr6.R")
wdcMerged_bytw <- wdcMerged

source("GetDataGooglePlaces_20dis_aggmonth_bytw_norealloc_averaged_pr6.R")
wdcMerged_bytw_norealloc <- wdcMerged


#percentage changes
summary((wdcMerged_bytw_norealloc$out_dem_mean-wdcMerged_bytw$out_dem_mean)/wdcMerged_bytw$out_dem_mean*100)
plot(sort((wdcMerged_bytw_norealloc$out_dem_mean-wdcMerged_bytw$out_dem_mean)/wdcMerged_bytw$out_dem_mean*100))





