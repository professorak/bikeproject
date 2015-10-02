eval_constraints_MPEC_tw_groupin <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  stocked_list <- which(wdcMergedday$stocked_out==FALSE)
  dem_T <- eval_lambda_new(deltain, theta1, wdcMergedday, points, tw_groupin)[stocked_list]
  ind <- which(dem_T<= 1e-32)
  dem_T[ind] <- 1e-32
  
  dem_hat_T <- wdcMergedday$out_dem_sum[stocked_list]/wdcMergedday$obs_weight[stocked_list]

  obj <- log(dem_T)- log(dem_hat_T)
  return(obj)
}

eval_grad_constraints_MPEC_tw_groupin <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  grad_all <- eval_grad_lambda_new(theta1, deltain, wdcMergedday, points, tw_groupin)
  grad_all <- grad_all[,-2] #remove rand coef col
  
  #could alternatively use version which returns lambda and grad_lambda_theta together.
  stocked_list <- which(wdcMergedday$stocked_out==FALSE)
  dem_T <- eval_lambda_new(deltain, theta1, wdcMergedday, points, tw_groupin)[stocked_list]
  
  #reduce the rows and columns corresponding to only stocked in obervations of wdcMergedday
  length_theta <- length(theta1)-1
  keep_cols <- c(c(1:length_theta),length_theta+stocked_list)
  grad_all <- grad_all[stocked_list,keep_cols]
  #accounting for log constraints
  grad_all <- grad_all/dem_T
  return(linearize_sparsify(grad_all))
}
