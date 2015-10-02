#get data
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_tiny.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse.R")

#run
print_iter_values <- 1
source("temp_optimalGMM_model5.R")
sourceCpp("eval_func_2_new_deltaaveraged_stepfunction_steps50_100_200.cpp")
source("eval_func_2_cpp_steps50_100_200.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps50_100_200.R")
source("eval_func_3_cpp_new_2.86_2_steps50_100_200.R")
source("eval_covariates_delta_reg_step_servlvl.R")

colnames_theta1 <<- c("dis coef","rand coef","dis coef","dis coef","dis coef",
                      "density ridership","density metro","density intercept",
                      "density metro evening", "density google places count",
                      "density bus")
nonden_ceofrange <- c(1:5)
nonden_ceoflength <- 5
density_ridership_col <<- 6
density_metro_col <<- 7
density_intercept_col <<- 8
density_metro_evening_col <<- 9
density_google_places_count <<- 10
density_bus_col <<- 11

weighing_GMM_mat <<- NULL

#########
wdcMerged$out_dem_sum <- wdcMerged$out_dem_sum/wdcMerged$obs_weight
wdcMerged$obs_weight <- 1
#########

length_stklist <- length(which(wdcMerged$stocked_out==F))
theta <- c(-4 ,-4,-4,-4, 1, 146.83, 20, 203.69, 0.42, 9.41)
#choose starting values for deltain_stin
# delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
# deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
delta_all <- rnorm(nrow(wdcMerged),-3,1)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(c(theta[1],0,theta[-1]),deltain_stin,wdcMerged,
#                                                       points)
length_eta <- ncol(eval_covariates_delta_reg(deltain_stin,c(theta[1],0,theta[-1]),wdcMerged,points)$Z)

#########
X <<- NULL
Z <<- NULL
Z <- eval_covariates_delta_reg(deltain_stin,c(theta[1],0,theta[-1]),wdcMerged,points)$Z
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
get_total_density(params, wdcMerged, points)
tot_density <<- 100000
target_density_constraint <<- 100000

lb = c(-rep(100,4),0,rep(0,5),rep(-20, length_delta), rep(-1.0e19,length_eta))
ub = c(rep(100,4),1000000,1000000,1000000,1000000,1000000,1000000,rep(10, length_delta), rep(1.0e19,length_eta))
constraint_lb <- c(rep(0, length_delta+length_eta), target_density_constraint) 
constraint_ub <- c(rep(0, length_delta+length_eta), target_density_constraint)

eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length(theta), length(deltain_stin), length_eta)
eval_hess_struct <- eval_hessian_structure(params,wdcMerged, points, length(theta), length(deltain_stin), length_eta)

#W_optimal <<- eval_weighting_mat_new(c(-10),wdcMerged,points) 

opts <- list("tol"=1.0e-6,
             "print_info_string"='yes',
             "derivative_test"="second-order",
             "max_iter"=1
)

print_iter_values <- 0

time <- proc.time()
res <- ipoptr( x0=params,
               lb=lb,
               ub=ub,
               eval_f=eval_obj_GMM_model5_obj, 
               eval_grad_f=eval_obj_GMM_model5_grad,
               eval_g=eval_g,
               eval_jac_g=eval_jac_g,
               eval_jac_g_structure=eval_jac_g_structure_val,
               eval_h=eval_hessian,
               eval_h_structure=eval_hess_struct,
               constraint_lb=constraint_lb,
               constraint_ub=constraint_ub,
               wdcMerged=wdcMerged, 
               points=points, 
               length_theta=length(theta),
               length_delta=length_delta,
               length_eta=length_eta,
               opts=opts
) 
time_1 <- proc.time() - time
print("time: ")
print(time_1)


