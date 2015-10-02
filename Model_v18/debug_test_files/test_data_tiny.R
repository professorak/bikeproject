#test time
#data
source("eval_obj_func_GMM_model4.R")
source("GetDataGooglePlaces_tiny.R")

############################################################
#run
theta <- c(-2.66047062428, 0.29046474623, 1681.87287960098, 3.43903446529, 
           836.49733133688, 0.31363526842, 0.01000001459)
length_theta <- length(theta)
set.seed(34675)
deltain_stin <- rnorm(length(which(!wdcMerged$stocked_out)),-3, 2)
params <- c(theta,deltain_stin)  

a1_eval_constraints <- eval_g_tiny(params,wdcMerged, points, length(theta))
a1_eval_jac_constraints <- eval_jac_g_tiny(params,wdcMerged, points, length(theta))
