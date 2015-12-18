ret <- eval_error_xi_model5(deltain,theta1,wdcMerged,points)
Z = ret$Z
weights = ret$weights
xi = ret$xi
Z_scaled = Z * weights

#moment conditoins:
G <- t(Z_scaled) %*% xi/ sum(weights)


#moment variance
moment_vec <- c()
for(i in 1:ncol(Z_scaled)) {
  moment_vec <- cbind(moment_vec,Z_scaled[,i]*xi)
}  
S0 = (t(moment_vec) %*% moment_vec) / sum(weights)

G[1]
S0[1,1]
G[1]*G[1]/S0[1,1]

S0_inv <- solve(S0)
t(G) %*% S0_inv %*% G
##################################################
#If instead the wieghts were such that they had mean of 1

ret <- eval_error_xi_model5(deltain,theta1,wdcMerged,points)
Z = ret$Z
weights = ret$weights
xi = ret$xi

weights <- weights/mean(weights)
Z_scaled = Z * weights

#moment conditoins:
G <- t(Z_scaled) %*% xi/ sum(weights)


#moment variance
moment_vec <- c()
for(i in 1:ncol(Z_scaled)) {
  moment_vec <- cbind(moment_vec,Z_scaled[,i]*xi)
}  
S0 = (t(moment_vec) %*% moment_vec) / sum(weights)

G[1]
S0[1,1]
G[1]*G[1]/S0[1,1]

S0_inv <- solve(S0)
t(G) %*% S0_inv %*% G


