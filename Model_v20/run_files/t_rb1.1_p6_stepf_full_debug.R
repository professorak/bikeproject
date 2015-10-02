source("constants_mw6.R")

#model with averaged station demand, no service level. Trying to get distance
#coefficient right.

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_aggmonth_averaged_pr6.R")

print_iter_values <<- 1

#run
source("temp_optimalGMM_model5.R")

sourceCpp("eval_func_2_new_deltaaveraged_stepfunction_steps_set4.cpp") 
#distance step at median quartile 0.116 kms
source("eval_func_2_cpp_steps50_100_200.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps50_100_200.R")
source("eval_func_3_cpp_new_2.86_2_steps50_100_200.R")
source("eval_covariates_delta_reg_densitycols_setb1.R") #all metro, googleplaces etc. density variables 
#and catchement area instruments + service level single step
source("moredensitycols_functions.R")

colnames_theta1 <<- c("dis coef","rand coef","dis coef",
                      "density ridership","density metro","density intercept",
                      "density metro evening","local_government_office", "food",
                      "density bus")
nonden_ceofrange <- c(1:3)
nonden_ceoflength <- 3
density_ridership_col <<- 4
density_metro_col <<- 5
density_intercept_col <<- 6
density_metro_evening_col <<- 7
#density_google_places_count <<- 10
density_bus_col <<- 10

density_google_places_theta_cols <- c(8,9)
density_google_places_points_cols <- c(10,11)



weighing_GMM_mat <<- NULL

length_stklist <- length(which(wdcMerged$stocked_out==F))
#theta <- c(-4 ,-4, 1, 146.83, 20, 0, 0.42, 0.42,0)
theta <- c(-29.6575325067742, -1.19717472440932, 0.499068754637771, 
           526.222361566426, 1.06544076377622, 0, 1.33679461513106, 
           0.191897784044254, 0)

active_coef <- c(1:length(theta))[-c(density_bus_col-1,density_metro_evening_col-1)]
active_coef_nointercept <- c(1:length(theta))[-c(density_intercept_col-1,density_bus_col-1,density_metro_evening_col-1)]


#choose starting values for deltain_stin
#delta_all <- rnorm(nrow(wdcMerged),-10,1)
delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(c(theta[1],0,theta[-1]),deltain_stin,wdcMerged,
#                                                       points)
length_eta <- ncol(eval_covariates_delta_reg(deltain_stin,c(theta[1],0,theta[-1]),wdcMerged,points)$Z)

#########
list_covariates <- eval_covariates_delta_reg(deltain_stin,c(theta[1],0,theta[-1]),wdcMerged,points)
X <<- list_covariates$X
Z <<- list_covariates$Z

stocked_list <- which(wdcMerged$stocked_out==FALSE)
Z_sqweighted <- (Z * sqrt(wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
print("kappa(t(Z_sqweighted) %*% Z_sqweighted): ")
print(kappa(t(Z_sqweighted) %*% Z_sqweighted))
weighing_GMM_mat <- solve(t(Z_sqweighted) %*% Z_sqweighted)

#########

eta <- rep(0,length_eta)
eta <- c(eval_constraints_moment_conditions(deltain_stin, c(theta[1],0,theta[-1]), eta, wdcMerged, points))
length_theta <- length(theta)
length_delta <- length(deltain_stin)
params <- c(theta,deltain_stin,eta)


#eval obj
eval_obj_GMM_model5_obj(params, wdcMerged, points, length(theta), length(deltain_stin), length_eta)

xi <- eval_error_xi_model5 (deltain_stin, 
                                  c(theta[1],0,theta[-1]), wdcMerged,points)$xi


#########
#comparing with the case when theta discoef is -2,-2

theta_case2 <- c(-2, -2, 0.499068754637771, 
           526.222361566426, 1.06544076377622, 0, 1.33679461513106, 
           0.191897784044254, 0)
delta_all_case2 <- compute_delta_list_cntrt_map_new(c(theta_case2[1],0,theta_case2[-1]),wdcMerged, points)
deltain_stin_case2 <- delta_all_case2[which(wdcMerged$stocked_out==F)]

eta_case2 <- rep(0,length_eta)
eta_case2 <- c(eval_constraints_moment_conditions(deltain_stin_case2, 
            c(theta_case2[1],0,theta_case2[-1]), eta_case2, wdcMerged, points))
params_case2 <- c(theta_case2,deltain_stin_case2,eta_case2)


eval_obj_GMM_model5_obj(params_case2, wdcMerged, points, length(theta), length(deltain_stin), length_eta)

xi_case2 <- eval_error_xi_model5 (deltain_stin_case2, 
              c(theta_case2[1],0,theta_case2[-1]), wdcMerged,points)$xi
##############
#comparing
eta_df <- cbind(eta, eta_case2, diag(weighing_GMM_mat))
rownames(eta_df) <- colnames(Z)
#are there outliers in xi values:


outliers_1 <- wdcMerged[which(xi-xi_case2<=-0.2),
          c("lat","lon","station_id","station_id_index","tract","out_dem_mean","catchment_area_step2","metro_den_on_1","metro_den_on_2",
            "googleplaces_den_1","googleplaces_food_den_1","googleplaces_grocery_den_1")]

summary(wdcMerged$out_dem_sum)
summary(wdcMerged$metro_den_on_1)
summary(wdcMerged$googleplaces_den_1)
summary(wdcMerged$googleplaces_food_den_1)
summary(wdcMerged$googleplaces_grocery_den_1)

#find the quantile of density values of above outliers:
which(sort(wdcMerged$metro_den_on_1)>max(outliers_1$metro_den_on_1))[1]/nrow(wdcMerged)
which(sort(wdcMerged$googleplaces_food_den_1)>max(outliers_1$googleplaces_food_den_1))[1]/nrow(wdcMerged)
which(sort(wdcMerged$googleplaces_grocery_den_1)>max(outliers_1$googleplaces_grocery_den_1))[1]/nrow(wdcMerged)
