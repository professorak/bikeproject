#compare linear fit when using logmetro corrected vs the one before

#run
source("temp_optimalGMM_model5.R")

#distance step at median quartile 0.300 kms
source("eval_func_2_cpp_steps_v0weights.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps_v0weights.R")
source("eval_func_3_cpp_new_2.86_2_steps_v0weights.R")
source("eval_covariates_delta_reg_densitycols_setb36.R") 
source("moredensitycols_functions.R")


get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_aggmonths_averaged_20districts_norealloc_stkoutthresh_five.R")

#log of absolute metro traffic
points_metro_idx <- which(points$type==2)
points$weight[points_metro_idx] <- log(points$weight[points_metro_idx]*22468468)

source("data_estimation_2.6_weather_functions.R")
#regenerate metro moment variables for wdcMerged
wdcMerged_newvars <- get_local_attributes_st_state(wdcMerged[,c("station_id_index","tract","tw","lat","lon","sto_state_local")], points)
list_metro_vars <- c("metro_den_on_1","metro_den_on_2","metro_den_off_1","metro_den_off_2",
                     "metro_den_on_3","metro_den_on_4","metro_den_off_3","metro_den_off_4",
                     "log(metro_den_on_1+1)","log(metro_den_on_2+1)","log(metro_den_on_3+1)","log(metro_den_on_4+1)",
                     "log(metro_den_off_1+1)","log(metro_den_off_2+1)","log(metro_den_off_3+1)","log(metro_den_off_4+1)",
                     "metro_den_on_a","metro_den_on_b","metro_den_on_c",
                     "metro_den_off_a","metro_den_off_b","metro_den_off_c")
wdcMerged[,list_metro_vars] <- wdcMerged_newvars[,list_metro_vars]
rm(wdcMerged_newvars)
wdcMerged_logmetro <- wdcMerged

source("GetDataGooglePlaces_aggmonths_averaged_20districts_norealloc_stkoutthresh_five.R")

source("debug_functions.R")
points <- assign_log_metro_traffic_corrected(points)

source("data_estimation_2.6_weather_functions.R")
#regenerate metro moment variables for wdcMerged
wdcMerged_newvars <- get_local_attributes_st_state(wdcMerged[,c("station_id_index","tract","tw","lat","lon","sto_state_local")], points)
list_metro_vars <- c("metro_den_on_1","metro_den_on_2","metro_den_off_1","metro_den_off_2",
                     "metro_den_on_3","metro_den_on_4","metro_den_off_3","metro_den_off_4",
                     "log(metro_den_on_1+1)","log(metro_den_on_2+1)","log(metro_den_on_3+1)","log(metro_den_on_4+1)",
                     "log(metro_den_off_1+1)","log(metro_den_off_2+1)","log(metro_den_off_3+1)","log(metro_den_off_4+1)",
                     "metro_den_on_a","metro_den_on_b","metro_den_on_c",
                     "metro_den_off_a","metro_den_off_b","metro_den_off_c")
wdcMerged[,list_metro_vars] <- wdcMerged_newvars[,list_metro_vars]
rm(wdcMerged_newvars)

wdcMerged_logmetrocorrected <- wdcMerged

delta_all <- rnorm(nrow(wdcMerged),-10,1)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
list_covariates <- eval_covariates_delta_reg(deltain_stin,theta1,wdcMerged,points)
X <<- list_covariates$X
Z <<- list_covariates$Z


plot(wdcMerged_logmetro$metro_den_on_a)
plot(wdcMerged_logmetrocorrected$metro_den_on_a)

vars <- c("metro_den_on_1","metro_den_on_2","metro_den_on_3","metro_den_on_4",
          "metro_den_on_a","metro_den_on_b","metro_den_on_c")
vars_index <- which(colnames(wdcMerged) %in% vars)
summary(lm(wdcMerged_logmetro$out_dem_mean ~ 0 + Z))


