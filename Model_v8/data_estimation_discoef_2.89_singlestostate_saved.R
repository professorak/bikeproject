#source("constants_mw6.R")
#v3.0 two stockstate observations per stations, if results are consistent, can experiment with other fixed effects etc.
point_range <- 3
no_bootstrap_runs <- 1
reg_runs_save_file <- "reg_runs_discoef_2.89_1states_mw6_bootstrap.txt"


source("eval_func_4_GMM_new_temp_2.89.R")
source("data_prep_functions_ffdf_limstostate.R")
source("data_prep_inte_points_limstostate.R")
source("eval_func_2_cpp_cntrt_map_2.86_2.R")
source("eval_func_3_cpp_new_2.86_2.R")


load(file="wdcMerged_singlestostate_lin.RData")
load(file="current_serv_lvl_singlestostate_lin.RData")
load(file="user_serv_lvl_singlestostate_lin.RData")
load(file="points_singlestostate_lin.RData")
v0_vec <- c(0)

x0_start <<- NULL
prev_theta <<- NULL


# #save commands
# save(wdcMerged, file="wdcMerged_singlestostate_lin.RData")
# save(current_serv_lvl, file="current_serv_lvl_singlestostate_lin.RData")
# save(user_serv_lvl, file="user_serv_lvl_singlestostate_lin.RData")
# save(points, file="points_singlestostate_lin.RData")


