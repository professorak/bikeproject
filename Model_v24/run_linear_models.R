#linear regressions

source("constants_mw6.R")
library("AER")

#least squares implementation 
#using incoming demand from previous tw as instrument
#using quadratic formulation for service level

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_20dis_aggmonth_hyperlocalstate_majoritystates_preweatherreg_moreinstr_withreallocthresh_stkoutthresh_five_averaged_pr6.R")

v0_vec <- c(0)
v0_vec_weights <- rep(1,length(v0_vec))
# load(file="v0_vec.RData")
# load(file="v0_vec_weights.RData")
points$lodging <- 0
points$museum <- 0
points$movie_theater <- 0
points$local_government_office <- 0
points$food <- 0
#log of absolute metro traffic
points_metro_idx <- which(points$type==2)
points$weight[points_metro_idx] <- log(points$weight[points_metro_idx]*22468468)

source("data_estimation_2.6_weather_functions.R")
#regenerate metro moment variables for wdcMerged
wdcMerged_newvars <- get_local_attributes_st_state(wdcMerged[,c("station_id_index","tract","tw","lat","lon","sto_state_local")], points)
list_metro_vars <- c("metro_den_on_1","metro_den_on_2","metro_den_off_1","metro_den_off_2",
                     "metro_den_on_3","metro_den_on_4","metro_den_off_3","metro_den_off_4",
                     "log(metro_den_on_1+1)","log(metro_den_on_2+1)","log(metro_den_on_3+1)","log(metro_den_on_4+1)",
                     "log(metro_den_off_1+1)","log(metro_den_off_2+1)","log(metro_den_off_3+1)","log(metro_den_off_4+1)",
                     "metro_den_on_a","metro_den_on_b","metro_den_on_c",
                     "metro_den_off_a","metro_den_off_b","metro_den_off_c")
wdcMerged[,list_metro_vars] <- wdcMerged_newvars[,list_metro_vars]
rm(wdcMerged_newvars)
#since we are adding dummy, we can bring lowest to 0 so that coefficient on this is positive.
points_metro_idx <- which(points$type==2)
points$weight[points_metro_idx] <- points$weight[points_metro_idx] - min(points$weight[points_metro_idx])


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
all_X_vars <- colnames(X)
density_vars <- colnames(X)[c(which(colnames(X) %like% "census_density"),
                  which(colnames(X) %like% "metro_den"),                  
                  which(colnames(X) %like% "touristlocs_den"))]
density_1_vars <- colnames(X)[c(which(colnames(X) %like% "census_density"),
                                which(colnames(X) %like% "metro_den_on_1"),
                                which(colnames(X) %like% "metro_den_off_1"),
                                which(colnames(X) %like% "touristlocs_den_1"))]
tract_tw_vars <- colnames(X)[which(colnames(X) %like% "tract_tw")]
distance_vars <- colnames(X)[which(colnames(X) %in% c("dis_nearest","dis_nearest_sq"))]
distance_5_vars <- colnames(X)[which(colnames(X) %in% c("dis_nearest_5","dis_nearest_5_sq"))]
other_vars <- colnames(X)[which(colnames(X) %like% "dis_center")]

if(length(all_X_vars)!=
     (1+length(density_vars)+length(distance_5_vars)+length(tract_tw_vars)+
        length(distance_vars)+length(other_vars))) stop("var length doesnt match")

#functions
get_var_tract_tw_effects <- function(fit) {
  coef_fit <- coef(fit)
  coef_tract_tw <- coef_fit[which(names(coef_fit) %like% "tract_tw_fac")]
  coef_tract_tw_df <- data.frame(coef_tract_tw=coef_tract_tw)
  pattern <- "(\\d+)"
  m <- regexpr(pattern, names(coef_tract_tw))
  coef_tract_tw_df$tract <- as.numeric(regmatches(names(coef_tract_tw), m))
  coef_tract_tw_0 <- -as.numeric(by(coef_tract_tw_df$coef_tract_tw, coef_tract_tw_df$tract, FUN=sum))
  coef_tract_tw_all <- c(coef_tract_tw_df$coef_tract_tw, coef_tract_tw_0)
  sum(coef_tract_tw_all)
  return(var(coef_tract_tw_all))
}


#OLS no density variables
summary(lm( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ 
              X[,"Intercept"] + X[,other_vars] + 
              X[,tract_tw_vars] + X[,distance_vars] 
            + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list]))
#OLS with density variables
fit <- lm( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ 
              X[,"Intercept"] +  X[,other_vars] + 
              X[,tract_tw_vars] + X[,distance_vars] + X[,density_vars]
            + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0, 
            weights=wdcMerged$obs_weight[stocked_list])
summary(fit)
get_var_tract_tw_effects(fit)

#IV with density variables
fit <- ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ 
                 X[,"Intercept"] +  X[,other_vars] + 
                 X[,tract_tw_vars] + X[,distance_vars] + X[,density_vars]
               + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0 |                 
                 X[,"Intercept"] +  X[,other_vars] + 
                 X[,tract_tw_vars] + X[,distance_vars] + X[,density_vars]
               + Instr_df[,c("in_dem_rate_lagtw","in_dem_rate_lagtw_sq")] + 0, 
               weights=wdcMerged$obs_weight[stocked_list])
summary(fit)
get_var_tract_tw_effects(fit)

#changing distance definition to nearest 5 avg.
summary(lm( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ 
              X[,"Intercept"] +  X[,other_vars] + 
              X[,tract_tw_vars] + X[,distance_5_vars] 
            + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0, 
            weights=wdcMerged$obs_weight[stocked_list]))

summary(lm( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ 
              X[,"Intercept"] +  X[,other_vars] + 
              X[,tract_tw_vars] + X[,distance_5_vars] + X[,density_vars]
            + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0, 
            weights=wdcMerged$obs_weight[stocked_list]))

summary(ivreg( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ 
                 X[,"Intercept"] +  X[,other_vars] + 
                 X[,tract_tw_vars] + X[,distance_5_vars] + X[,density_vars]
               + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0 |                 
                 X[,"Intercept"] +  X[,other_vars] + 
                 X[,tract_tw_vars] + X[,distance_5_vars] + X[,density_vars]
               + Instr_df[,c("in_dem_rate_lagtw","in_dem_rate_lagtw_sq")] + 0, 
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

