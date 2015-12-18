idx <- which(wdcMerged$out_dem_mean==0.00001)
wdcMerged$out_dem_mean[idx] <- min(wdcMerged$out_dem_mean[-idx])

#OLS with weight=obs_weight
fit1 <- lm( log(wdcMerged$out_dem_mean)[stocked_list] ~ 
             X[,"Intercept"] +  X[,other_vars] + 
             X[,tract_tw_vars] + X[,distance_vars] + X[,density_vars]
           + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0, 
           weights=wdcMerged$obs_weight[stocked_list])
summary(fit1)


#OLS with weight=1
fit2 <- lm( log(wdcMerged$out_dem_mean)[stocked_list] ~ 
              X[,"Intercept"] +  X[,other_vars] + 
              X[,tract_tw_vars] + X[,distance_vars] + X[,density_vars]
            + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0)
summary(fit2)


#OLS with weight=obs_weight^2
fit3 <- lm( log(wdcMerged$out_dem_mean)[stocked_list] ~ 
              X[,"Intercept"] +  X[,other_vars] + 
              X[,tract_tw_vars] + X[,distance_vars] + X[,density_vars]
            + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0, 
            weights=(wdcMerged$obs_weight^2)[stocked_list])
summary(fit3)




#OLS with weight=obs_weight^3
fit4 <- lm( log(wdcMerged$out_dem_mean)[stocked_list] ~ 
              X[,"Intercept"] +  X[,other_vars] + 
              X[,tract_tw_vars] + X[,distance_vars] + X[,density_vars]
            + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0, 
            weights=log(wdcMerged$obs_weight[stocked_list]))
summary(fit4)


#OLS with weight=log(obs_weight)
fit5 <- lm( log(wdcMerged$out_dem_mean)[stocked_list] ~ 
              X[,"Intercept"] +  X[,other_vars] + 
              X[,tract_tw_vars] + X[,distance_vars] + X[,density_vars]
            + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0, 
            weights=log(wdcMerged$obs_weight[stocked_list]))
summary(fit5)

#only obs with weight more than 100
idx <- which(wdcMerged$obs_weight>100)
fit6 <- lm( log(wdcMerged$out_dem_mean)[idx] ~ 
              X[idx,"Intercept"] +  X[idx,other_vars] + 
              X[idx,tract_tw_vars] + X[idx,distance_vars] + X[idx,density_vars]
            + Instr_df[idx,c("serv_lvl","serv_lvl_sq")] + 0, 
            weights=wdcMerged$obs_weight[idx])
summary(fit6)

#only obs with weight more than 100
idx <- which(wdcMerged$obs_weight>100)
fit7 <- lm( log(wdcMerged$out_dem_mean)[idx] ~ 
              X[idx,"Intercept"] +  X[idx,other_vars] + 
              X[idx,tract_tw_vars] + X[idx,distance_vars] + X[idx,density_vars]
            + Instr_df[idx,c("serv_lvl","serv_lvl_sq")] + 0)
summary(fit7)

#weights above 1000 truncated to 1000
weights <- wdcMerged$obs_weight
weights[which(weights>1000)] <- 1000
fit8 <- lm( log(wdcMerged$out_dem_mean)[stocked_list] ~ 
              X[,"Intercept"] +  X[,other_vars] + 
              X[,tract_tw_vars] + X[,distance_vars] + X[,density_vars]
            + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0, 
            weights=log(wdcMerged$obs_weight[stocked_list]))
summary(fit8)



