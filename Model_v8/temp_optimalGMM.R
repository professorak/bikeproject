#optimal weighing matrix
eval_GMM_optimal_weighing_matrix <- function(theta1,deltain,wdcMerged,points) {
  ret <- eval_error_xi_model4(theta1,deltain,wdcMerged,points)
  Z = ret$Z
  weights = ret$weights
  xi = ret$xi
  Z_scaled = Z * weights
  #we need to calculate S_0^(-1)=1/N \sum(G_i*G_i') 
  #where G_i is column vector of moments for each observation
  #i.e. in my setting it is Z_scaled[i,]*xi[i]
  #\sum(G_i*G_i') can be also calculated as product of a 
  #matrix[G1,G2,..] i.e. where G1,G2 etc are stacked next to each on other in ith columns
  #..and its transpose.
  #we construct t(moment_vec) which is the matrix above.
  #Note the corresponding row weights are multipled to Z already in creating Z_scaled
  moment_vec <- c()
  for(i in 1:ncol(Z_scaled)) {
    moment_vec <- cbind(moment_vec,Z_scaled[,i]*xi)
  }  
  S0 = (t(moment_vec) %*% moment_vec) / sum(weights)
  S0_inv <- solve(S0)
  return(S0_inv)  
}


