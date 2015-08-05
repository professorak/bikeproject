sourceCpp("eval_func_2_new_deltaaveraged.cpp")

source("GetDataGooglePlaces_tiny.R")


weighing_GMM_mat <<- NULL

length_stklist <- length(which(wdcMerged$stocked_out==F))
theta <- c(-4 ,0, 146.83, 20, 203.69, 0.42, 9.41)
#choose starting values for deltain_stin
delta_all <- rnorm(nrow(wdcMerged), -10,1)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(c(theta[1],0,theta[-1]),deltain_stin,wdcMerged,
#                                                       points)
params <- c(theta,deltain_stin)


deltain <- rep(-30, nrow(wdcMerged))
deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
theta1 <- c(theta[1],0,theta[-1])
stocked_list <- which(wdcMerged$stocked_out==F)

tw_group_list <- unique(wdcMerged$tw_group)
i<- 1
tw_groupin = tw_group_list[i]
wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]


################################################################
################################################################


compute_hess_beta1sq <- function(constraint_index, lambda_multiplers_in) {
  a1 <- eval_lambda_new(deltain_tw,theta1, wdcMergedday, points, tw_groupin)[constraint_index]
  diff_vec <- rep(0, length(theta1))
  diff_vec[1] <- diff
  #gradient wrt delta1_index at base value
  a1_2 <- eval_lambda_new(deltain_tw,theta1+diff_vec, wdcMergedday, points, tw_groupin)[constraint_index]
  a1_grad <- (a1_2-a1)/diff
  
  a2 <- eval_lambda_new(deltain_tw, theta1+diff_vec, wdcMergedday, points, tw_groupin)[constraint_index]
  #a2_hess <- eval_hessian_delta_list_new(deltain_tw_3, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in)
  
  a2_2 <- eval_lambda_new(deltain_tw, theta1+2*diff_vec, wdcMergedday, points, tw_groupin)[constraint_index]
  a2_grad <- (a2_2-a2)/diff
  
  
  ###
  a1_hess_num <- (a2_grad-a1_grad)/diff
  return(a1_hess_num)
  
}

compute_hess_theta1sq_atindex <- function(constraint_index, theta1_index, theta2_index, lambda_multiplers_in) {
  a1 <- eval_lambda_new(deltain_tw,theta1, wdcMergedday, points, tw_groupin)[constraint_index]
  diff1_vec <- rep(0, length(theta1))
  diff1_vec[theta1_index] <- diff
  diff2_vec <- rep(0, length(theta1))
  diff2_vec[theta2_index] <- diff
  #gradient wrt delta1_index at base value
  a1_2 <- eval_lambda_new(deltain_tw,theta1+diff1_vec, wdcMergedday, points, tw_groupin)[constraint_index]
  a1_grad <- (a1_2-a1)/diff
  
  a2 <- eval_lambda_new(deltain_tw, theta1+diff2_vec, wdcMergedday, points, tw_groupin)[constraint_index]
  #a2_hess <- eval_hessian_delta_list_new(deltain_tw_3, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in)
  
  a2_2 <- eval_lambda_new(deltain_tw, theta1+diff2_vec+diff1_vec, wdcMergedday, points, tw_groupin)[constraint_index]
  a2_grad <- (a2_2-a2)/diff
  
  
  ###
  a1_hess_num <- (a2_grad-a1_grad)/diff
  return(a1_hess_num)
  
}

#test hessian 

constraint_index <- 4
diff <- 0.001
lambda_multiplers_in <- rep(0,nrow(wdcMergedday))
lambda_multiplers_in[constraint_index] <- 1
theta1_index <- 1
theta2_index <- 1

a1_hess <- eval_hessian_beta1_sq_tw(deltain_tw, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in)
a1_hess
a1_hess_num <- compute_hess_beta1sq(constraint_index, lambda_multiplers_in)
a1_hess_num

a2_hess <- eval_hessian_theta1_sq_tw(deltain_tw, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in)
a2_hess
a2_hess_num <- compute_hess_theta1sq_atindex(constraint_index, theta1_index, theta2_index, lambda_multiplers_in)
a2_hess_num

a2_hess_num <- matrix(0,length(theta1),length(theta1))
for(i in 1:length(theta1)) {
  for(j in 1:i) {
    theta1_index <- i
    theta2_index <- j
    a2_hess_num[i,j] <- compute_hess_theta1sq_atindex(constraint_index, theta1_index, theta2_index, lambda_multiplers_in)
    a2_hess_num[j,i] <- a2_hess_num[i,j]
  }
}





