theta2 <- c(-3.68,0,1,0)
a2 <- eval_obj_GMM_list_extended_new (theta2[-2], wdcMerged, points)  
deltain2 <- prev_deltain 
ret2 <- eval_error_xi(theta2, wdcMerged, points)

xi2=ret2$xi
xi_org2 <- xi2/sqrt(wdcMerged$obs_weight)
grad_delta_theta2 = ret2$grad_delta_theta

obj_2 = c(t(xi2 ) %*% xi2 )/length(xi2)
idx2 <- which(wdcMerged$stocked_out==F)
png(file="AnalysisPlots/test_xi_values_xi2.png")
  plot(xi2[idx2], cex=0.1, xlab="index",ylab="xi")
dev.off()
png(file="AnalysisPlots/test_xi_values_xi_org2.png")
  plot(xi_org2[idx2], cex=0.1, xlab="index",ylab="xi_org")
dev.off()
png(file="AnalysisPlots/test_xi_values_deltain2.png")
  plot(deltain2[idx2], cex=0.3, xlab="index",ylab="deltain")
dev.off()
png(file="AnalysisPlots/test_xi_values_obs_weight2.png")
  plot(wdcMerged$obs_weight[idx2], cex=0.3, xlab="index",ylab="obs_weight")
dev.off()




theta1 <- c(-3.68,0,18.33,110)
a1 <- eval_obj_GMM_list_extended_new (theta1[-2], wdcMerged, points)  
deltain1 <- prev_deltain 
ret1 <- eval_error_xi(theta1, wdcMerged, points)

xi1=ret1$xi
xi_org1 <- xi1/sqrt(wdcMerged$obs_weight)
grad_delta_theta1 = ret1$grad_delta_theta

obj_1 = c(t(xi1 ) %*% xi1 )/length(xi1)
idx1 <- which(wdcMerged$stocked_out==F)
png(file="AnalysisPlots/test_xi_values_xi1.png")
plot(xi1[idx1], cex=0.1, xlab="index",ylab="xi")
dev.off()
png(file="AnalysisPlots/test_xi_values_xi_org1.png")
plot(xi_org1[idx1], cex=0.1, xlab="index",ylab="xi_org")
dev.off()
png(file="AnalysisPlots/test_xi_values_deltain1.png")
plot(deltain1[idx1], cex=0.3, xlab="index",ylab="deltain")
dev.off()
png(file="AnalysisPlots/test_xi_values_obs_weight1.png")
plot(wdcMerged$obs_weight[idx1], cex=0.3, xlab="index",ylab="obs_weight")
dev.off()





