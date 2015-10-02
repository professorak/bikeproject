source("constants_mw6.R")
source("debug_functions.R")
#model with averaged station demand, no service level. Trying to get distance
#coefficient right.

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_aggmonth_averaged_pr6.R")

print_iter_values <<- 1

#run
source("temp_optimalGMM_model5.R")

sourceCpp("eval_func_2_new_deltaaveraged_stepfunction_steps_set7.cpp") 
#distance step at median quartile 0.300 kms
source("eval_func_2_cpp_steps50_100_200.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps50_100_200.R")
source("eval_func_3_cpp_new_2.86_2_steps50_100_200.R")
source("eval_covariates_delta_reg_densitycols_setb17.R") #all metro, googleplaces etc. density variables minus census density
#and catchement area instruments + service level single step
source("moredensitycols_functions.R")

colnames_theta1 <<- c("dis coef","rand coef","dis coef",
                      "density ridership","density metro","density intercept",
                      "density metro evening","local_government_office", "food",                      
                      "density bus", "density tourist locations")
nonden_ceofrange <- c(1:3)
nonden_ceoflength <- 3
density_ridership_col <<- 4
density_metro_col <<- 5
density_intercept_col <<- 6
density_metro_evening_col <<- 7
#density_google_places_count <<- 10
#density_bus_col <<- 10
density_add_den_cols <<- c(10,11) #bus, tourist location


density_google_places_theta_cols <- c(8,9)
density_google_places_points_cols <- c(which(colnames(points)=="local_government_office"),
                                       which(colnames(points)=="food"))
)


weighing_GMM_mat <<- NULL

length_stklist <- length(which(wdcMerged$stocked_out==F))
theta <- c(-15.3761554117194, -7.75959563777998,
           0, 0, 0.29086, 0, 0, 0, 0, 0)

active_coef <- c(1:length(theta))[-c(density_add_den_cols[1]-1,density_metro_evening_col-1,density_ridership_col-1)]
active_coef_nointercept <- c(1:length(theta))[-c(density_intercept_col-1,density_add_den_cols[1]-1,density_metro_evening_col-1,
                                                 density_ridership_col-1)]


#choose starting values for deltain_stin
delta_all <- rep(-6.424,nrow(wdcMerged))
#delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(c(theta[1],0,theta[-1]),deltain_stin,wdcMerged,
#                                                       points)

get_total_density(theta, wdcMerged, points)

dem <- eval_lambda_full(deltain_stin, get_theta1fromtheta(theta), wdcMerged, points)

plot(dem, type="l")


