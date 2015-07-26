#test_gradients eval_obj_func_GMM_model4
source("eval_obj_func_GMM_model4.R")
source("GetDataGooglePlaces.R")



######
#test
#test gradient manually for fewer dimensions
#theta <- c(-10.7374110745,0.8860478015,403.4177015742,2.8258972200, 200,10,10)
theta <- c(1.91,0.76,148.54,0.0768,206.8,0.25,7.799,0.036)

set.seed(34675)
deltain_stin <- rnorm(length(which(!wdcMerged$stocked_out)),-3, 2)
# delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
# deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
params <- c(theta,deltain_stin)  

#test objective gradient
a1 <- eval_obj_GMM_model4_obj(params,wdcMerged, points, length(theta))
a1_g <- eval_obj_GMM_model4_grad(params,wdcMerged, points, length(theta))

idx <- 10
diff <- 0.01

params_2 <- params
params_2[idx] <- params_2[idx] + diff
a2 <- eval_obj_GMM_model4_obj(params_2,wdcMerged, points, length(theta))
a2_g <- eval_obj_GMM_model4_grad(params_2,wdcMerged, points, length(theta))

(a2-a1)/diff
a2_g[idx]
a1_g[idx]

#test constraints gradient
a1_eval_constraints <- eval_g(params,wdcMerged, points, length(theta))
get_total_density(params[1:length(theta)], wdcMerged, points)
a1_eval_jac_constraints <- eval_jac_g(params,wdcMerged, points, length(theta))
eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length(theta))
if(length(a1_eval_jac_constraints)!=length(unlist(eval_jac_g_structure_val))) stop("gradient lengths dont match a1_eval_jac_constraints")
a1_eval_jac_constraints <- getfullfromsparsematrix (eval_jac_g_structure_val, a1_eval_jac_constraints) 


idx <- 6
diff <- 0.0001

params_2 <- params
params_2[idx] <- params_2[idx] + diff  

a2_eval_constraints <- eval_g(params_2,wdcMerged, points, length(theta))
get_total_density(params_2[1:length(theta)], wdcMerged, points)
a2_eval_jac_constraints <- eval_jac_g(params_2,wdcMerged, points, length(theta))
a2_eval_jac_constraints <- getfullfromsparsematrix (eval_jac_g_structure_val, a2_eval_jac_constraints) 

(sum(a2_eval_constraints)-sum(a1_eval_constraints))/diff
sum(a2_eval_jac_constraints[,idx])
sum(a1_eval_jac_constraints[,idx])
identical(round((a2_eval_constraints-a1_eval_constraints)/diff,2),
          round(a2_eval_jac_constraints[,idx],2))

((a2_eval_constraints-a1_eval_constraints)/diff)[1:10]
(a2_eval_jac_constraints[1:10,idx])
a1_eval_jac_constraints[1:10,idx]
num_grad <- (a2_eval_constraints-a1_eval_constraints)/diff
identical(round((a2_eval_constraints-a1_eval_constraints)/diff,2),
          round(a2_eval_jac_constraints[,idx],2))
num_grad[which(round(num_grad,2)!=
                 round(a2_eval_jac_constraints[,idx],2))][1:10]
a2_eval_jac_constraints[,idx][which(round(num_grad,2)!=
                                      round(a2_eval_jac_constraints[,idx],2))][1:10]



#########################
#test gradient wrt dis coef, delta optimized
theta <- c(1.91,0.76,148.54,0.0768,206.8,0.25,7.799)
delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
params <- c(theta,deltain_stin)  
a1 <- eval_obj_GMM_model4_obj(params,wdcMerged, points, length(theta))
a1_eval_constraints <- eval_g(params,wdcMerged, points, length(theta))

theta_2 <- c(-2,0.76,148.54,0.0768,206.8,0.25,7.799)
delta_all_2 <- compute_delta_list_cntrt_map_new(c(theta_2[1],0,theta_2[-1]),wdcMerged, points)
deltain_stin_2 <- delta_all_2[which(wdcMerged$stocked_out==F)]
params_2 <- c(theta_2,deltain_stin_2)  
a1_2 <- eval_obj_GMM_model4_obj(params_2,wdcMerged, points, length(theta_2))


#########################
#test why xi values are not the same for a constnt change in delta.
theta <- c(1.91,0.76,148.54,0.0768,206.8,0.25,7.799)
theta1 <- c(theta[1],0,theta[-1])
delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
params <- c(theta,deltain_stin)  
ret <- eval_error_xi_model4(theta1,deltain_stin,wdcMerged,points)
xi <- ret$xi

deltain_stin_2 <- deltain_stin+ 1
#params <- c(theta,deltain_stin)  
ret_2 <- eval_error_xi_model4(theta1,deltain_stin_2,wdcMerged,points)
xi_2 <- ret_2$xi







