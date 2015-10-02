
# scale_G_hat_theta1 <<- sqrt(diag((t(G_hat_theta1) %*% G_hat_theta1)))
# G_hat_theta1_save <- G_hat_theta1
# 
# G_hat_theta1 <- scalecols(G_hat_theta1, scale_G_hat_theta1)
# kappa(t(G_hat_theta1) %*% G_hat_theta1)

G_hat_theta1 <- grad_delta %*% grad_delta_theta_full
G_hat_theta2 <- - t(Z_weighted) %*% X
G_hat <- cbind(G_hat_theta1, G_hat_theta2)
print(paste0("kappa t(G_hat) %*% A_N %*% G_hat: ",kappa(t(G_hat) %*% A_N %*% G_hat)))

var_p1 <- solve(t(G_hat) %*% A_N %*% G_hat)
thetafull_variance <- (var_p1 %*% t(G_hat) %*% A_N %*% S0 %*% A_N %*% 
                         G_hat %*% var_p1)/sum(weights)

thetafull_variance[2,2]

var_p1 <- solve(t(G_hat_theta1) %*% A_N %*% G_hat_theta1)
theta1_variance <- (var_p1 %*% t(G_hat_theta1) %*% A_N %*% S0 %*% A_N %*% 
                      G_hat_theta1 %*% var_p1)/sum(weights)

theta1_variance[1,1]


G_hat_save <- G_hat
scale_G_hat <<- sqrt(diag((t(G_hat) %*% G_hat)))

G_hat <- scalecols(G_hat, scale_G_hat)
print(paste0("kappa t(G_hat) %*% A_N %*% G_hat: ",kappa(t(G_hat) %*% A_N %*% G_hat)))
kappa(t(G_hat) %*% G_hat)

var_p1 <- solve(t(G_hat) %*% A_N %*% G_hat)
thetafull_variance_2 <- (var_p1 %*% t(G_hat) %*% A_N %*% S0 %*% A_N %*% 
                         G_hat %*% var_p1)/sum(weights)

thetafull_variance_2[2,2]/(scale_G_hat[2]*scale_G_hat[2])

