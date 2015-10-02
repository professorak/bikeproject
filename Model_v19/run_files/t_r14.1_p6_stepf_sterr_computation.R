source("constants_mw6.R")

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_aggmonth_pr6.R")

print_iter_values <<- 1

#run
source("temp_optimalGMM_model5.R")

sourceCpp("eval_func_2_new_deltaaveraged_stepfunction_steps_set5.cpp") #distance step at 1st quartile 0.076810 kms
source("eval_func_2_cpp_steps50_100_200.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps50_100_200.R")
source("eval_func_3_cpp_new_2.86_2_steps50_100_200.R")
source("eval_covariates_delta_reg_stepset4_servlvl_densitycols_set3.R")
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

theta <- c(-32.15637679, -66.31767490,   0.22106725, 291.70979719,   1.96944162,
           291.60925357,   1.00669628,   0.05990265,   0.00000000)
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

params_save <- params
res_save <- res
theta <- params[c(1:length(theta))]
theta_save <- theta
deltain_stin <- params[length(theta)+c(1:length_delta)]
eta <- params[length(theta)+length_delta+c(1:length_eta)]
theta1 <- c(theta[1],0,theta[-1])  
eval_xi <- eval_error_xi_model5 (deltain_stin, theta1, wdcMerged,points)
theta2_save <- eval_xi$theta2  
theta2_save

# var_covar_theta <- eval_theta_variance(theta1,deltain_stin,wdcMerged,points)
source("temp_optimalGMM_model5_allstderr.R")
active_coef <- c(1:length_theta)[-c(density_bus_col-1)]
active_coef_nointercept <- c(1:length_theta)[-c(density_intercept_col-1,density_bus_col-1)]
var_covar_theta_save <- eval_theta_variance(theta1,deltain_stin,wdcMerged,points,active_coef)
round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta_save)),2)
var_covar_theta_2_save <- eval_theta2_variance(theta1,deltain_stin,wdcMerged,points,active_coef)


weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(theta1,deltain_stin,wdcMerged,
                                                      points)


##step 2
theta <- c(-54.4458359,4.5244449,0.7346146,456.7708460,0.5844075,
           1201.8140911,4.2887575,0.1506534,0.0000000)
#choose starting values for deltain_stin
#delta_all <- rnorm(nrow(wdcMerged),-10,1)
delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
eta <- rep(0,length_eta)
eta <- c(eval_constraints_moment_conditions(deltain_stin, c(theta[1],0,theta[-1]), eta, wdcMerged, points))
params <- c(theta,deltain_stin,eta)




theta <- params[c(1:length(theta))]
deltain_stin <- params[length(theta)+c(1:length_delta)]
theta1 <- c(theta[1],0,theta[-1])  
eval_xi <- eval_error_xi_model5 (deltain_stin, theta1, wdcMerged,points)
theta2 <- eval_xi$theta2  
theta2
theta3 <- eval_xi$theta3  
theta3



var_covar_theta <- eval_theta_variance(theta1,deltain_stin,wdcMerged,points,active_coef)
round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta)),2)

st_err <- sqrt(diag(var_covar_theta))
round(correlation_theta_estimates <- diag(1/st_err) %*% var_covar_theta %*% diag(1/st_err),3)

var_covar_theta_2 <- eval_theta2_variance(theta1,deltain_stin,wdcMerged,points,active_coef)
serv_cols <- which(rownames(theta2) %like% "serv_lvl")
round(theta2_save[serv_cols]/sqrt(diag(var_covar_theta_2)[serv_cols]),2)
st_err_theta2 <- sqrt(diag(var_covar_theta_2))
round(correlation_theta_2_estimates <- diag(1/st_err_theta2[serv_cols]) %*% var_covar_theta_2[serv_cols,serv_cols] %*% 
        diag(1/st_err_theta2[serv_cols]),3)

serv_lvl_scale <- scalecols_X[rownames(theta2)[serv_cols]]
##########
print("Step 1 output: ")
print(theta_save)
print(res_save$objective)
print(theta2_save[which(rownames(theta2_save) %like% "serv_lvl")]/serv_lvl_scale)
print("theta1 t-stat:")
print(round(theta_save[active_coef_nointercept]/sqrt(diag(var_covar_theta_save)),2))
print("theta2 t-stat:")
print(round(theta2_save[serv_cols]/sqrt(diag(var_covar_theta_2_save)[serv_cols]),2))
print(res_save$iterations)
print(time_1)
print("Step 2 output: ")
print(theta)
print(res$objective)
print(theta2[which(rownames(theta2) %like% "serv_lvl")]/serv_lvl_scale)
print("theta1 t-stat:")
print(round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta)),2))
print("theta2 t-stat:")
print(round(theta2[serv_cols]/sqrt(diag(var_covar_theta_2)[serv_cols]),2))
print(res$iterations + res_save$iterations)
print(time + time_1)
print(round(correlation_theta_estimates <- diag(1/st_err) %*% var_covar_theta %*% diag(1/st_err),3))
print(round(correlation_theta_2_estimates,3))


