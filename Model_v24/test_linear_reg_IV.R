#linear regressions

library("AER")

#get data
source("GetDataGooglePlaces_20dis_aggmonth_hyperlocalstate_moreinstr_withreallocthresh_stkoutthresh_five_averaged_pr6.R")
source("eval_covariates_linear_reg.R")

# #collapse data into aggregated observations (could have also used bytw files but 
# #better to stick to a single source)
# wdcMerged_linear_reg_df$tot_obs_weight <- 
#   ave(wdcMerged$obs_weight, wdcMerged$)

stocked_list <- which(wdcMerged$stocked_out==FALSE)

lin_covariates <- eval_covariates_linear_reg(wdcMerged)
X <- lin_covariates$X
#Z <- lin_covariates$Z
Instr_df <- wdcMerged[stocked_list,c("serv_lvl", "instr_serv_lvl", "serv_lvl_neighbours", "in_dem_rate",
   "out_dem_rate", "in_dem_rate_lagtw", "out_dem_rate_lagtw", "diff_dem_rate_lagtw",
   "diffdummy_dem_rate_lagtw")]
Instr_df$serv_lvl_sq = Instr_df$serv_lvl*Instr_df$serv_lvl
Instr_df$serv_lvl_neighbours_sq = Instr_df$serv_lvl_neighbours*Instr_df$serv_lvl_neighbours
Instr_df$in_dem_rate_sq = Instr_df$in_dem_rate * Instr_df$in_dem_rate
Instr_df$in_dem_rate_lagtw_sq = Instr_df$in_dem_rate_lagtw * Instr_df$in_dem_rate_lagtw

Instr_df$diff_dem_rate_lagtw_sq <- Instr_df$diff_dem_rate_lagtw * Instr_df$diff_dem_rate_lagtw
Instr_df <- as.matrix(Instr_df)
#   X <- cbind(Xbase, serv_lvl=wdcMerged$serv_lvl[stocked_list], 
#              serv_lvl_sq=wdcMerged$serv_lvl[stocked_list]*wdcMerged$serv_lvl[stocked_list])


# reg_vars <- c("metro_den_on_1","metro_den_on_2","metro_den_off_1","metro_den_off_2",
#               "serv_lvl")
reg_vars <- colnames(X)
#reg_vars_Z <- colnames(Z)

#80.7 no instruments
summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
               + Instr_df[,c("serv_lvl")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))

summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
               + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))


#76.7 service level neighobours
summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
  + Instr_df[,c("serv_lvl")] + 0|
  X[,reg_vars]  + Instr_df[,c("serv_lvl_neighbours")] + 0, 
  weights=wdcMerged$obs_weight[stocked_list]))

summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
     + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0 |
     X[,reg_vars]  + Instr_df[,c("serv_lvl_neighbours","serv_lvl_neighbours_sq")] + 0, 
     weights=wdcMerged$obs_weight[stocked_list]))

#78.7 in dem rate
summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
               + Instr_df[,c("serv_lvl")] + 0|
                 X[,reg_vars]  + Instr_df[,c("in_dem_rate")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))

summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
               + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0 |
                 X[,reg_vars]  + Instr_df[,c("in_dem_rate","in_dem_rate_sq")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))

#74.7 service level neighobours
summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
               + Instr_df[,c("serv_lvl")] + 0|
                 X[,reg_vars]  + Instr_df[,c("serv_lvl_neighbours","in_dem_rate")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))

summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
               + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0 |
                 X[,reg_vars]  + Instr_df[,c("serv_lvl_neighbours","in_dem_rate",
                "serv_lvl_neighbours_sq","in_dem_rate_sq")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))

#82.7 previous tw incoming demand rate
summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
               + Instr_df[,c("serv_lvl")] + 0|
                 X[,reg_vars]  + Instr_df[,c("in_dem_rate_lagtw")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))

summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
               + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0 |
                 X[,reg_vars]  + Instr_df[,c("in_dem_rate_lagtw","in_dem_rate_lagtw_sq")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))

#instrument strength
summary(lm( Instr_df[,c("serv_lvl")] ~ X[,reg_vars] 
               + Instr_df[,c("in_dem_rate_lagtw")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))
summary(lm( Instr_df[,c("serv_lvl")] ~ X[,reg_vars] 
            + Instr_df[,c("in_dem_rate_lagtw","in_dem_rate_lagtw_sq")] + 0, 
            weights=wdcMerged$obs_weight[stocked_list]))
summary(lm( Instr_df[,c("serv_lvl_sq")] ~ X[,reg_vars] 
            + Instr_df[,c("serv_lvl","in_dem_rate_lagtw","in_dem_rate_lagtw_sq")] + 0, 
            weights=wdcMerged$obs_weight[stocked_list]))



#84.7 previous tw diff incoming outgoing demand rate
summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
               + Instr_df[,c("serv_lvl")] + 0|
                 X[,reg_vars]  + Instr_df[,c("diff_dem_rate_lagtw")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))

summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
               + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0 |
                 X[,reg_vars]  + Instr_df[,c("diff_dem_rate_lagtw","diff_dem_rate_lagtw_sq")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))


#86.7 previous tw diff dummy incoming outgoing demand rate
summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ X[,reg_vars] 
               + Instr_df[,c("serv_lvl")] + 0|
                 X[,reg_vars]  + Instr_df[,c("diffdummy_dem_rate_lagtw")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))



#record the results


#do tests for strenght of instruments and weak instruments





