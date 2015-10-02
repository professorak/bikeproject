source("constants_mw6.R")
sourceCpp("eval_func_2_new_deltaaveraged.cpp", verbose = T)

source("GetDataGooglePlaces_tiny.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse.R")
source("add_eval_func_2_cpp_MPEC_model5_nonsparse.R")

source("test_hessian_implementation_functions.R")
source("eval_obj_GMM_model5.R")

weighing_GMM_mat <<- NULL

theta <- c(-4 ,1, 146.83, 0.052, 203.69, 0.42, 9.41)

set.seed(34675)
deltain_stin <- rnorm(length(which(!wdcMerged$stocked_out)),-10, 1)
eta <- rep(0,12)
length_eta <- length(eta)
length_theta <- length(theta)
length_delta <- length(deltain_stin)
params <- c(theta,deltain_stin, eta)

deltain <- rep(-30, nrow(wdcMerged))
deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
theta1 <- c(theta[1],0,theta[-1])
stocked_list <- which(wdcMerged$stocked_out==F)

lambda_multiplers_in <- rep(0,length(stocked_list))
lambda_multiplers_in[c(100,101)] <- 1

# tw_group_list <- unique(wdcMerged$tw_group)
# i<- 1
# tw_groupin = tw_group_list[i]
# wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
# deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]
# lambda_multiplers_full <- rep(0, nrow(wdcMerged))
# lambda_multiplers_full[which(wdcMerged$stocked_out==F)] <- lambda_multiplers_in
# lambda_multiplers_tw = lambda_multiplers_full[which(wdcMerged$tw_group==tw_groupin)]


################################################################
################################################################

#check if the linearization of hessian is happening properly.
#Do this by first creating linear hessian in a full-proof way, ie. by creating
#entire matrix and then linearizing.


# full_hess_lambda <- test_eval_hessian_lambda_constraints_MPEC(deltain_stin, theta1, 
#                                                        eta, wdcMerged, points, lambda_multiplers_in)
# 
# a1_hess_lambda <- eval_hessian_lambda_constraints_MPEC(deltain_stin, theta1, 
#                                      eta, wdcMerged, points, lambda_multiplers_in)

eval_hess_struct <- eval_hessian_structure(params,wdcMerged, points, length(theta), length(deltain_stin), length_eta)

obj_factor <- 0
constraint_multipliers <- rep(0,length_delta+length_eta+1)
constraint_multipliers[c(1)] <- 1
diff <- 0.01

a1_hess <- eval_hessian(params,wdcMerged, points, length(theta), length(deltain_stin), length_eta, obj_factor, constraint_multipliers)
if(length(a1_hess)!=length(unlist(eval_hess_struct))) stop("hessian lengths dont match a1_eval_jac_constraints")
a1_hess_full <- getfullfromsparsematrix (eval_hess_struct, a1_hess) 

index1 <- c(1,10,30,50,100,1000,2000,2790)
index2 <- c(1,10,30,50,100,1000,2000,2790)
a1_num_hess_at_index <- matrix(0, length(index1), length(index2))
for(i in 1:length(index1)) {
  for(j in 1:length(index2)) {
    a1_num_hess_at_index[i,j] <- compute_hess_atindex(params, index1[i],index2[j], obj_factor, constraint_multipliers)
  }
}

a1_hess_full[index1,index2]
a1_num_hess_at_index
a <- rbind(a1_hess_full[index1,index2],a1_num_hess_at_index)

#compute hessian at [1,1] individually for each constraint

#for(i in 1:c(length_delta+length_eta+1)) {
obj_factor <- 0
diff <- 0.00001
hessian_compare_values <- c()
for(i in 1:10) {
  constraint_multipliers <- rep(0,length_delta+length_eta+1)
  constraint_multipliers[c(i)] <- 1
  
  a1_hess <- eval_hessian(params,wdcMerged, points, length(theta), length(deltain_stin), length_eta, obj_factor, constraint_multipliers)
  if(length(a1_hess)!=length(unlist(eval_hess_struct))) stop("hessian lengths dont match a1_eval_jac_constraints")
  a1_hess_full <- getfullfromsparsematrix (eval_hess_struct, a1_hess) 
  
  a1_num_hess_at_index <-  compute_hess_atindex(params, 1, 1, obj_factor, constraint_multipliers)
  hessian_compare_values <- rbind(hessian_compare_values, c(a1_hess_full[1,1],a1_num_hess_at_index))
}

