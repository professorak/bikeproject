source("constants_mw6.R")
#_googlevar2=food

#least squares implementation 
#using incoming demand from previous tw as instrument
#using quadratic formulation for service level

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_20dis_aggmonth_hyperlocalstate_majoritystates_preweatherreg_moreinstr_withreallocthresh_stkoutthresh_five_averaged_pr6.R")

v0_vec <- c(0)
v0_vec_weights <- rep(1,length(v0_vec))
# load(file="v0_vec.RData")
# load(file="v0_vec_weights.RData")

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
#since we are adding dummy, we can bring lowest to 0 so that coefficient on this is positive.
points_metro_idx <- which(points$type==2)
points$weight[points_metro_idx] <- points$weight[points_metro_idx] - min(points$weight[points_metro_idx])

print_iter_values <<- 1

#run
source("temp_optimalGMM_model5.R")

sourceCpp("eval_func_2_hyperlocal_nldis_v0weights_set2.cpp")  
#distance step at median quartile 0.300 kms
source("eval_func_2_cpp_steps_v0weights.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps_v0weights.R")
source("eval_func_3_cpp_new_2.86_2_steps_v0weights.R")
source("eval_covariates_delta_reg_densitycols_setb50.6_googlevar2_wservlvl_metrooff_metrodummies.R")  
source("moredensitycols_functions_metrodummies.R")

colnames_theta1 <<- c("dis coef","rand coef","dis coef",
                      "density ridership","density metro","density intercept",
                      "density metro evening","local_government_office", "food",
                      "lodging","museum","movie_theater",
                      "density bus", "density tourist locations","density metro dummy","density metro evening dummy")
nonden_ceofrange <- c(1:3)
nonden_ceoflength <- 3
density_ridership_col <<- 4
density_metro_col <<- 5
density_intercept_col <<- 6
density_metro_evening_col <<- 7
#density_google_places_count <<- 10
#density_bus_col <<- 10
density_add_den_cols <<- c(13,14,15,16) #bus, tourist location, metro dummies


density_google_places_theta_cols <- c(8,9,10,11,12)
density_google_places_points_cols <- c(which(colnames(points)=="local_government_office"),
                                       which(colnames(points)=="food"),
                                       which(colnames(points)=="lodging"),
                                       which(colnames(points)=="museum"),
                                       which(colnames(points)=="movie_theater")
)




weighing_GMM_mat <<- NULL

length_stklist <- length(which(wdcMerged$stocked_out==F))
theta1 <- c(-2.20961045069201, 0, -19.6523343222199, 0.00523699070208834, 2.12663960984702, 
            0.131087591513272, 1.53182163420034, 0, 0, 0, 0, 0, 
            0, 86.6566512013304, 3.78614232898306, 3.12139107706004)

active_coef <- c(1:length(theta1))[-c(2,density_google_places_theta_cols[-2],
                                      density_add_den_cols[1])]
active_coef_nointercept <- c(1:length(theta1))[-c(2,density_intercept_col,
                                                  density_google_places_theta_cols[-2],density_add_den_cols[1])]



#choose starting values for deltain_stin
#delta_all <- rnorm(nrow(wdcMerged),-10,1)
delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points, 1e-6, 1e-6)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]

list_covariates <- eval_covariates_delta_reg(deltain_stin,theta1,wdcMerged,points)
X <<- list_covariates$X
Z <<- list_covariates$Z
stocked_list <- which(wdcMerged$stocked_out==FALSE)
Z_sqweighted <- (Z * sqrt(wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
print("kappa(t(Z_sqweighted) %*% Z_sqweighted): ")
print(kappa(t(Z_sqweighted) %*% Z_sqweighted))
weighing_GMM_mat <- solve(t(Z_sqweighted) %*% Z_sqweighted)

save(delta_all,file="t_r82.7.1_noRD_getdeltavalues_delta_all.RData")

low_del_idx <- length(which(delta_all <= -11))
sum(wdcMerged$obs_weight[low_del_idx])
sum(wdcMerged$obs_weight[-low_del_idx])

#test for heteroskedasticity
xi <- eval_error_xi_model5(deltain_stin, theta1, wdcMerged, points)$xi
plot(xi)
plot(xi[order(wdcMerged$obs_weight[stocked_list])])

