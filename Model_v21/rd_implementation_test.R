source("constants_mw6.R")

#model with averaged station demand, no service level. Trying to get distance
#coefficient right.

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_aggmonth_averaged_pr6.R")
#v0_vec <- generate_v0(5)
v0_vec <- c(0,0)
  
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
theta1 <- c(-4 ,0,-4, 0, 146.83, 20, 0, rep(0.42, 5),0,1)

active_coef <- c(1:length(theta1))[-c(density_add_den_cols[1],density_metro_evening_col,density_ridership_col)]
active_coef_nointercept <- c(1:length(theta1))[-c(density_intercept_col,density_add_den_cols[1],density_metro_evening_col,
                                                 density_ridership_col)]



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
length_theta1 <- length(theta1)
length_delta <- length(deltain_stin)
params <- c(theta1,deltain_stin,eta)
get_total_density(params, wdcMerged, points)
tot_density <<- 100000
target_density_constraint <<- 100000

eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
eval_hess_struct <- eval_hessian_structure(params,wdcMerged, points, length_theta1, length(deltain_stin), length_eta)



v_eval_g <- eval_g(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
summary(v_eval_g)
v_eval_jac_g <- eval_jac_g(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
v_eval_jac_g[1:10]
summary(v_eval_jac_g)
length(v_eval_jac_g)

v_eval_obj_GMM_model5_obj <- eval_obj_GMM_model5_obj(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
v_eval_obj_GMM_model5_obj

v_eval_obj_GMM_model5_grad <- eval_obj_GMM_model5_grad(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
summary(v_eval_obj_GMM_model5_grad)

set.seed(1)
obj_factor <- rnorm(1,0,1)  
hessian_lambda <- rnorm(length(v_eval_g),0,1)  
v_eval_hessian <- eval_hessian(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta,
                               obj_factor,hessian_lambda)
v_eval_hessian[1:10]
summary(v_eval_hessian)
length(v_eval_hessian)

