

#######
#test gradients of delta
theta <- c(-4 ,0, 146.83, 20, 203.69, 0.42, 9.41)
theta1 <- c(theta[1],0,theta[-1])
deltain <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
grad_delta_theta_full <- eval_grad_delta_theta_full(theta1, deltain, wdcMerged, points)

idx <- 5
diff <- 0.0001

theta_2 <- theta
theta_2[idx] <- theta_2[idx]+diff
theta1_2 <- c(theta_2[1],0,theta_2[-1])
deltain_2 <- compute_delta_list_cntrt_map_new(theta1_2,wdcMerged, points)
grad_delta_theta_full_2 <- eval_grad_delta_theta_full(theta1_2, deltain_2, wdcMerged, points)

deltain_list <- cbind(deltain_2, deltain)

summary(grad_delta_theta_full_numerical <- (deltain_2-deltain)/diff)
summary(grad_delta_theta_full_2[,idx])
summary(grad_delta_theta_full[,idx])
length(which(abs(grad_delta_theta_full_numerical-grad_delta_theta_full_2[,idx])>=1e-4))


