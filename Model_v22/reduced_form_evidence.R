source("constants_mw6.R")

#compute substitution percentage
#do weighted linear regression. See the effect of stocked nearby 
#after controlling for st_tw effect.


source("GetDataGooglePlaces_20dis_aggmonth_hyperlocalstate_majoritystates_preweatherreg_moreinstr_withreallocthresh_stkoutthresh_five_averaged_pr6.R")

library(lfe)
fit <- lm(log(out_dem_mean) ~ sto_nearby,
          data=wdcMerged, weights=wdcMerged$obs_weight)

summary(fit)


fit_felm <- felm(log(out_dem_mean) ~ sto_nearby | st_tw,
                 data=wdcMerged, weights=wdcMerged$obs_weight)
summary(fit_felm)

fit_felm <- felm(log(out_dem_mean) ~ sto_nearby | st_tw,
                 data=wdcMerged, weights=sqrt(wdcMerged$obs_weight))
summary(fit_felm)

