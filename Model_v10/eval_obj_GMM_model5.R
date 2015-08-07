eval_obj_GMM_model5_obj <- function(params, wdcMerged, points, length_theta, length_delta, length_eta) {   
  print("objective comput")  
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[c((length_theta+1):(length_theta+length_delta))]
  eta <- params[c((length_theta+length_delta+1):(length_theta+length_delta+length_eta))]
  theta1 <- c(theta[1],0,theta[-1])  
  print("theta:")
  print(theta)
  
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(length_eta)    
  } else {
    A_N <- weighing_GMM_mat
  }
  
  print("obj:")
  return(print(eta %*% A_N %*% eta))
}  

eval_obj_GMM_model5_grad <- function(params, wdcMerged, points, length_theta, length_delta, length_eta) {    
  print("gradint comput")
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[c((length_theta+1):(length_theta+length_delta))]
  eta <- params[c((length_theta+length_delta+1):(length_theta+length_delta+length_eta))]
  theta1 <- c(theta[1],0,theta[-1])  
  print("theta:")
  print(theta)
  
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(length_eta)    
  } else {
    A_N <- weighing_GMM_mat
  }

  
  grad_eta <- c(2*A_N %*% eta)
  grad <- c(rep(0,length_theta+length_delta), grad_eta)
  return(grad)
} 
