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
v0_vec <- generate_v0(40)
v0_vec_weights <- rep(1,length(v0_vec))
points$lodging <- 0
points$museum <- 0
points$movie_theater <- 0
points$local_government_office <- 0
points$food <- 0

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


print_iter_values <<- 1

#run
source("temp_optimalGMM_model5.R")

sourceCpp("eval_func_2_new_deltaaveraged_nldis_v0weights_set2.cpp") 
#distance step at median quartile 0.300 kms
source("eval_func_2_cpp_steps_v0weights.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps_v0weights.R")
source("eval_func_3_cpp_new_2.86_2_steps_v0weights.R")
source("eval_covariates_delta_reg_densitycols_setb36.R") 
source("moredensitycols_functions.R")

colnames_theta1 <<- c("dis coef","rand coef","dis coef",
                      "density ridership","density metro","density intercept",
                      "density metro evening","local_government_office", "food",
                      "lodging","museum","movie_theater",
                      "density bus", "density tourist locations")
nonden_ceofrange <- c(1:3)
nonden_ceoflength <- 3
density_ridership_col <<- 4
density_metro_col <<- 5
density_intercept_col <<- 6
density_metro_evening_col <<- 7
#density_google_places_count <<- 10
#density_bus_col <<- 10
density_add_den_cols <<- c(13,14) #bus, tourist location


density_google_places_theta_cols <- c(8,9,10,11,12)
density_google_places_points_cols <- c(which(colnames(points)=="local_government_office"),
                                       which(colnames(points)=="food"),
                                       which(colnames(points)=="lodging"),
                                       which(colnames(points)=="museum"),
                                       which(colnames(points)=="movie_theater")
)




weighing_GMM_mat <<- NULL

length_stklist <- length(which(wdcMerged$stocked_out==F))
theta1 <- c(-3.054947e+00,8.871711e-04, -3.881124e+01,5.191147e-04,1.469369e-01,2.443626e-01 ,0.000000e+00,0.000000e+00,
            0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.579567e+00)

get_total_density(theta1, wdcMerged, points )
sum(wdcMerged$out_dem_mean)/get_total_density(theta1, wdcMerged, points )*100

#choose starting values for deltain_stin
#delta_all <- rnorm(nrow(wdcMerged),-10,1)
delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(theta1,deltain_stin,wdcMerged,
#                                                       points)
length_eta <- ncol(eval_covariates_delta_reg(deltain_stin,theta1,wdcMerged,points)$Z)

#########
list_covariates <- eval_covariates_delta_reg(deltain_stin,theta1,wdcMerged,points)
X <<- list_covariates$X
Z <<- list_covariates$Z

stocked_list <- which(wdcMerged$stocked_out==FALSE)
Z_sqweighted <- (Z * sqrt(wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
print("kappa(t(Z_sqweighted) %*% Z_sqweighted): ")
print(kappa(t(Z_sqweighted) %*% Z_sqweighted))
weighing_GMM_mat <- solve(t(Z_sqweighted) %*% Z_sqweighted)

#########

eta <- rep(0,length_eta)
eta <- c(eval_constraints_moment_conditions(deltain_stin, theta1, eta, wdcMerged, points))
length_theta <- length(theta1)
length_delta <- length(deltain_stin)
params <- c(theta1,deltain_stin,eta)


#eval obj
obj_1 <- eval_obj_GMM_model5_obj(params, wdcMerged, points, length_theta, length(deltain_stin), length_eta)

xi <- eval_error_xi_model5 (deltain_stin, 
                            theta1, wdcMerged,points)$xi




source("GetDataGooglePlaces_aggmonths_averaged_20districts_norealloc_stkoutthresh_five.R")
v0_vec <- generate_v0(40)
v0_vec_weights <- rep(1,length(v0_vec))
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

print_iter_values <<- 1

#run
source("temp_optimalGMM_model5.R")

sourceCpp("eval_func_2_new_deltaaveraged_nldis_v0weights_set2.cpp") 
#distance step at median quartile 0.300 kms
source("eval_func_2_cpp_steps_v0weights.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps_v0weights.R")
source("eval_func_3_cpp_new_2.86_2_steps_v0weights.R")
source("eval_covariates_delta_reg_densitycols_setb36.R") 
source("moredensitycols_functions.R")

colnames_theta1 <<- c("dis coef","rand coef","dis coef",
                      "density ridership","density metro","density intercept",
                      "density metro evening","local_government_office", "food",
                      "lodging","museum","movie_theater",
                      "density bus", "density tourist locations")
nonden_ceofrange <- c(1:3)
nonden_ceoflength <- 3
density_ridership_col <<- 4
density_metro_col <<- 5
density_intercept_col <<- 6
density_metro_evening_col <<- 7
#density_google_places_count <<- 10
#density_bus_col <<- 10
density_add_den_cols <<- c(13,14) #bus, tourist location


density_google_places_theta_cols <- c(8,9,10,11,12)
density_google_places_points_cols <- c(which(colnames(points)=="local_government_office"),
                                       which(colnames(points)=="food"),
                                       which(colnames(points)=="lodging"),
                                       which(colnames(points)=="museum"),
                                       which(colnames(points)=="movie_theater")
)




weighing_GMM_mat <<- NULL

length_stklist <- length(which(wdcMerged$stocked_out==F))
theta1 <- c(-2.087358e+00, 7.160281e-04, -4.287322e+01, 3.533148e-04, 3.357873e-01, 2.471581e-01, 
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 3.675670e+00)


get_total_density(theta1, wdcMerged, points )
sum(wdcMerged$out_dem_mean)/get_total_density(theta1, wdcMerged, points )*100

#choose starting values for deltain_stin
#delta_all <- rnorm(nrow(wdcMerged),-10,1)
delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(theta1,deltain_stin,wdcMerged,
#                                                       points)
length_eta <- ncol(eval_covariates_delta_reg(deltain_stin,theta1,wdcMerged,points)$Z)

#########
list_covariates <- eval_covariates_delta_reg(deltain_stin,theta1,wdcMerged,points)
X <<- list_covariates$X
Z <<- list_covariates$Z

stocked_list <- which(wdcMerged$stocked_out==FALSE)
Z_sqweighted <- (Z * sqrt(wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
print("kappa(t(Z_sqweighted) %*% Z_sqweighted): ")
print(kappa(t(Z_sqweighted) %*% Z_sqweighted))
weighing_GMM_mat <- solve(t(Z_sqweighted) %*% Z_sqweighted)

#########

eta <- rep(0,length_eta)
eta <- c(eval_constraints_moment_conditions(deltain_stin, theta1, eta, wdcMerged, points))
length_theta <- length(theta1)
length_delta <- length(deltain_stin)
params <- c(theta1,deltain_stin,eta)


#eval obj
obj_2 <- eval_obj_GMM_model5_obj(params, wdcMerged, points, length_theta, length(deltain_stin), length_eta)

xi_2 <- eval_error_xi_model5 (deltain_stin, 
                            theta1, wdcMerged,points)$xi



