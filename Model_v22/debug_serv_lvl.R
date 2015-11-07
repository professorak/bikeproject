#checking if there is too high correlation with service level.
#also how good are the instruments

source("GetDataGooglePlaces_20dis_aggmonth_bytw_norealloc_stkoutthresh_five_averaged_pr6.R")

cor(wdcMerged$out_dem_mean, wdcMerged$serv_lvl)
cor(wdcMerged$out_dem_mean, wdcMerged$instr_serv_lvl)
cor(wdcMerged$out_dem_mean, wdcMerged$serv_lvl_neighbours)
cor(wdcMerged$serv_lvl, wdcMerged$instr_serv_lvl)
cor(wdcMerged$serv_lvl, wdcMerged$serv_lvl_neighbours)

serv_lvl_order <- order(wdcMerged$tw, wdcMerged$serv_lvl)
plot(wdcMerged$out_dem_mean[serv_lvl_order])
plot(wdcMerged$serv_lvl[serv_lvl_order])

summary(lm(wdcMerged$out_dem_mean ~ wdcMerged$serv_lvl + wdcMerged$tw_fac))

cor(wdcMerged$metro_den_on_1, wdcMerged$serv_lvl)

wdcMerged$serv_lvl_dis <- floor(wdcMerged$serv_lvl*10)/10
a <- by(wdcMerged$out_dem_mean, wdcMerged[,c("serv_lvl_dis","tw")], FUN=mean)


