#test_gradients eval_obj_func_GMM_model4
#source("eval_obj_func_GMM_model4.R")
source("GetDataGooglePlaces_tiny.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse.R")


######
#test
#test gradient manually for fewer dimensions
#theta <- c(-10.7374110745,0.8860478015,403.4177015742,2.8258972200, 200,10,10)
theta <- c(-4 ,0, 146.83, 0.052, 203.69, 0.42, 9.41)

set.seed(34675)
deltain_stin <- rnorm(length(which(!wdcMerged$stocked_out)),-3, 2)
# delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
# deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
eta <- rep(0,12)
length_eta <- length(eta)
length_theta <- length(theta)
length_delta <- length(deltain_stin)
params <- c(theta,deltain_stin, eta)
# a <- eval_jac_g(params, wdcMerged, points, length(theta), length(deltain_stin), length_eta)
# length(c(a))
# length(which(a!=0))
# length(which(a!=0))/length(c(a))*100


#test constraints gradient
a1_eval_constraints <- eval_g(params,wdcMerged, points, length(theta), length(deltain_stin), length_eta)
get_total_density(params[1:length(theta)], wdcMerged, points)
a1_eval_jac_constraints <- eval_jac_g(params,wdcMerged, points, length(theta), length(deltain_stin), length_eta)

#idx <- 6
idx <- length(theta)+9
#idx <- length(theta)+length(deltain_stin)+9
diff <- 0.000001

params_2 <- params
params_2[idx] <- params_2[idx] + diff  

a2_eval_constraints <- eval_g(params_2,wdcMerged, points, length(theta), length(deltain_stin), length_eta)
get_total_density(params_2[1:length(theta)], wdcMerged, points)
a2_eval_jac_constraints <- eval_jac_g(params_2,wdcMerged, points, length(theta), length(deltain_stin), length_eta)

(sum(a2_eval_constraints)-sum(a1_eval_constraints))/diff
sum(a2_eval_jac_constraints[,idx])
sum(a1_eval_jac_constraints[,idx])
which(round((a2_eval_constraints-a1_eval_constraints)/diff,2)!=
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


