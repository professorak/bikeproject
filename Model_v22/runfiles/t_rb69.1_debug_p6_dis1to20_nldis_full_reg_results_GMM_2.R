moment_list <- c(1:ncol(Z))

moment_list <- c(1:ncol(Z))[-which(colnames(Z) %in% c("metro_den_on_3","metro_den_on_4","touristlocs_den_3","touristlocs_den_4"))]
colnames(Z)[moment_list]

moment_list <- c(1:ncol(Z))[which(colnames(Z) %in% c("touristlocs_den_1", "touristlocs_den_2", "touristlocs_den_3", "touristlocs_den_4"))]
moment_list <- c(1:ncol(Z))[which(colnames(Z) %in% c("metro_den_on_1","metro_den_on_2","metro_den_on_3","metro_den_on_4"))]

t(eta[moment_list]) %*% weighing_GMM_mat[moment_list,moment_list] %*% eta[moment_list]
t(eta_case2[moment_list]) %*% weighing_GMM_mat[moment_list,moment_list] %*% eta_case2[moment_list]

t(xi) %*% wdcMerged$metro_den_on_a
t(xi_case2) %*% wdcMerged$metro_den_on_a

t(xi) %*% wdcMerged$metro_den_on_b
t(xi_case2) %*% wdcMerged$metro_den_on_b

t(xi) %*% wdcMerged$metro_den_on_c
t(xi_case2) %*% wdcMerged$metro_den_on_c

t(xi) %*% wdcMerged$metro_den_on_
t(xi_case2) %*% wdcMerged$metro_den_on_d


t(xi) %*% wdcMerged$metro_den_on_a
t(xi_case2) %*% wdcMerged$metro_den_on_a

t(xi) %*% wdcMerged$metro_den_on_1
t(xi_case2) %*% wdcMerged$metro_den_on_1

#linear reg with 
summary(lm(wdcMerged$out_dem_mean ~ Z + (wdcMerged$metro_den_on_1>0)))
