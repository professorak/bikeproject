eval_obj_GMM_model5_obj <- function(params, wdcMerged, points, length_theta, length_delta, length_eta) {   
#   print("objective comput")  
  theta1 <- params[c(1:length_theta)]
  deltain_stin <- params[c((length_theta+1):(length_theta+length_delta))]
  eta <- params[c((length_theta+length_delta+1):(length_theta+length_delta+length_eta))]
  
  if(print_iter_values) {
    print("params: theta, summary of delta, summary of eta")
    print(theta1)
    print(as.numeric(summary(deltain_stin,digits = 10)))
    print(as.numeric(summary(eta,digits = 10)))
    print("")    
  }
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(length_eta)    
  } else {
    A_N <- weighing_GMM_mat
  }
  
#   print("obj:")
#   return(print(eta %*% A_N %*% eta))
    return((eta %*% A_N %*% eta))
}  

eval_obj_GMM_model5_grad <- function(params, wdcMerged, points, length_theta, length_delta, length_eta) {    
#   print("gradint comput")
  theta1 <- params[c(1:length_theta)]
  deltain_stin <- params[c((length_theta+1):(length_theta+length_delta))]
  eta <- params[c((length_theta+length_delta+1):(length_theta+length_delta+length_eta))]
#   print("theta:")
#   print(theta)
  
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(length_eta)    
  } else {
    A_N <- weighing_GMM_mat
  }

  
  grad_eta <- c(2*A_N %*% eta)
  grad <- c(rep(0,length_theta+length_delta), grad_eta)
  return(grad)
} 
