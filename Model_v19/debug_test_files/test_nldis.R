source("constants_mw6.R")



#get data
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_tiny.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse.R")

#run
print_iter_values <- 1
source("temp_optimalGMM_model5.R")

weighing_GMM_mat <<- NULL

#########
wdcMerged$out_dem_sum <- wdcMerged$out_dem_sum/wdcMerged$obs_weight
wdcMerged$obs_weight <- 1
#########

length_stklist <- length(which(wdcMerged$stocked_out==F))
theta <- c(-4, 0, 146.83, 20, 203.69, 0.42, 9.41)
#choose starting values for deltain_stin
# delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
# deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
delta_all <- rnorm(nrow(wdcMerged),-10,1)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(c(theta[1],0,theta[-1]),deltain_stin,wdcMerged,
#                                                       points)
length_eta <- ncol(eval_covariates_delta_reg(deltain_stin,c(theta[1],0,theta[-1]),wdcMerged,points)$Z)

#########
Z <- eval_covariates_delta_reg(deltain_stin,c(theta[1],0,theta[-1]),wdcMerged,points)$Z
stocked_list <- which(wdcMerged$stocked_out==FALSE)
Z_weighted <- (Z * (wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
print("kappa(t(Z_weighted) %*% Z_weighted): ")
print(kappa(t(Z_weighted) %*% Z_weighted))
weighing_GMM_mat <- solve(t(Z_weighted) %*% Z_weighted)

#########

eta <- rep(0,length_eta)

eta <- c(eval_constraints_moment_conditions(deltain_stin, c(theta[1],0,theta[-1]), eta, wdcMerged, points))

length_theta <- length(theta)
length_delta <- length(deltain_stin)
params <- c(theta,deltain_stin,eta)



###########
#test
obj_factor <- 1
hessian_lambda <- rep(1, length_delta+length_eta+1  )
sourceCpp("eval_func_2_new_deltaaveraged_nldis.cpp")
colnames_theta1 <<- c("dis coef","rand coef","nldis_coef","density ridership","density metro","density intercept",
                      "density metro evening", "density google places count",
                      "density bus")
density_ridership_col <<- 4
density_metro_col <<- 5
density_intercept_col <<- 6
density_metro_evening_col <<- 7
density_google_places_count <<- 8
density_bus_col <<- 9

theta <- c(-4 ,0, 0, 146.83, 20, 203.69, 0.42, 9.41)
length_theta <- length(theta)
params <- c(theta,deltain_stin,eta)
# a <- eval_hessian(params, wdcMerged, points, length_theta, length_delta, length_eta, obj_factor, hessian_lambda)
# eval_hessian_struct <- eval_hessian_structure(params, wdcMerged, points, length_theta, length_delta, length_eta)
# a_full <- getfullfromsparsematrix (eval_hessian_struct, a) 
a <- eval_jac_g(params, wdcMerged, points, length_theta, length_delta, length_eta)
eval_jac_struct <- eval_jac_g_structure(params, wdcMerged, points, length_theta, length_delta, length_eta)
(length(a)==length(unlist(eval_jac_struct)))
a_full <- getfullfromsparsematrix (eval_jac_struct, a) 



sourceCpp("eval_func_2_new_deltaaveraged.cpp")
colnames_theta1 <<- c("dis coef","rand coef","density ridership","density metro","density intercept",
                      "density metro evening", "density google places count",
                      "density bus")
density_ridership_col <<- 3
density_metro_col <<- 4
density_intercept_col <<- 5
density_metro_evening_col <<- 6
density_google_places_count <<- 7
density_bus_col <<- 8

theta <- c(-4 , 0, 146.83, 20, 203.69, 0.42, 9.41)
length_theta <- length(theta)
params <- c(theta,deltain_stin,eta)
# b <- eval_hessian(params, wdcMerged, points, length_theta, length_delta, length_eta, obj_factor, hessian_lambda)
# eval_hessian_struct <- eval_hessian_structure(params, wdcMerged, points, length_theta, length_delta, length_eta)
# b_full <- getfullfromsparsematrix (eval_hessian_struct, b) 

b <- eval_jac_g(params, wdcMerged, points, length_theta, length_delta, length_eta)
eval_jac_struct <- eval_jac_g_structure(params, wdcMerged, points, length_theta, length_delta, length_eta)
b_full <- getfullfromsparsematrix (eval_jac_struct, b) 

#identical(a_full[-2,-2], b_full)
identical(a_full[,-2], b_full)



####################################################
####################################################
#test gradient:
sourceCpp("eval_func_2_new_deltaaveraged_nldis.cpp")
colnames_theta1 <<- c("dis coef","rand coef","nldis_coef","density ridership","density metro","density intercept",
                      "density metro evening", "density google places count",
                      "density bus")
density_ridership_col <<- 4
density_metro_col <<- 5
density_intercept_col <<- 6
density_metro_evening_col <<- 7
density_google_places_count <<- 8
density_bus_col <<- 9

theta <- c(-4 ,1, 0, 146.83, 20, 203.69, 0.42, 9.41)
#theta <- c(-10.7374110745,0.8860478015,403.4177015742,2.8258972200, 200)
set.seed(34675)
deltain_stin <- rnorm(length(which(!wdcMerged$stocked_out)),-3, 0.5)
params <- c(theta,deltain_stin,eta)

a1_eval_constraints <- eval_g(params, wdcMerged, points, length_theta, length_delta, length_eta)
a1_eval_jac_constraints <- eval_jac_g(params, wdcMerged, points, length_theta, length_delta, length_eta)
eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length_theta, length_delta, length_eta)
a1_eval_jac_constraints <- getfullfromsparsematrix (eval_jac_g_structure_val, a1_eval_jac_constraints) 

idx <- 1
diff <- 0.0001

params_2 <- params
params_2[idx] <- params_2[idx] + diff  

a2_eval_constraints <- eval_g(params_2,wdcMerged, points, length_theta, length_delta, length_eta)
a2_eval_jac_constraints <- eval_jac_g(params_2,wdcMerged, points, length_theta, length_delta, length_eta)
a2_eval_jac_constraints <- getfullfromsparsematrix (eval_jac_g_structure_val, a2_eval_jac_constraints) 

(sum(a2_eval_constraints)-sum(a1_eval_constraints))/diff
sum(a2_eval_jac_constraints[,idx])
sum(a1_eval_jac_constraints[,idx])
identical(round((a2_eval_constraints-a1_eval_constraints)/diff,2),
          round(a2_eval_jac_constraints[,idx],2))

(er_id <- which(round((a2_eval_constraints-a1_eval_constraints)/diff,2)!=
          round(a2_eval_jac_constraints[,idx],2)))
round((a2_eval_constraints-a1_eval_constraints)/diff,2)[er_id][1:10]
round(a2_eval_jac_constraints[,idx],2)[er_id][1:10]



