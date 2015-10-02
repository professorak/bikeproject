source("temp_optimalGMM_model5_allstderr.R")
active_coef <- c(1:length_theta)[-c(density_bus_col-1)]
active_coef_nointercept <- c(1:length_theta)[-c(density_intercept_col-1,density_bus_col-1)]
var_covar_theta_save <- eval_theta_variance(theta1,deltain_stin,wdcMerged,points,active_coef)
round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta_save)),2)
var_covar_theta_2_save <- eval_theta2_variance(theta1,deltain_stin,wdcMerged,points,active_coef)
print(round(theta2_save[serv_cols]/sqrt(diag(var_covar_theta_2_save)[serv_cols]),2))

diag(var_covar_theta_save)
diag(var_covar_theta_2_save)

source("debug_optimalGMM_model5_allstderr_1.R")

var_covar_thetafull_save <- eval_thetafull_variance(theta1,deltain_stin,wdcMerged,points,active_coef)
var_covar_theta_save_new <- var_covar_thetafull_save[c(1:length(active_coef_nointercept)),
                                                     c(1:length(active_coef_nointercept))]
var_covar_theta_2_save_new <- var_covar_thetafull_save[-c(1:length(active_coef_nointercept)),
                                                       -c(1:length(active_coef_nointercept))]

round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta_save_new)),2)
print(round(theta2_save[serv_cols]/sqrt(diag(var_covar_theta_2_save_new)[serv_cols]),2))


