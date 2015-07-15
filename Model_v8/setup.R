
#PLEASE RUN constants.R before this to setup path variables.

#file to do intial setup

setwd(work_dir)


#source("data_prep_functions_ffdf.R")
#source("data_prep_inte_points.R")
source("eval_func_2_cpp.R")
#source("eval_func_2.R")
#source("eval_func_3_cpp.R")
#source("eval_func_4.R")
#source("eval_func_5_cpp.R")
source("util_functions.R")
#source("eval_func_7.R")

#source("data_prep_main.R")

#sourceCpp(file="eval_func_2.cpp")
#sourceCpp(file="eval_func_2_new.cpp")
sourceCpp(file="eval_func_2_new_deltaaveraged.cpp")
#source("eval_func_2_cpp_log.R")
#source("eval_func_2_cpp_cntrt_map.R")

#source("eval_func_4_GMM.R")
source("counterfactual_functions_new.R")

#source("eval_func_3_cpp_new.R")
#source("eval_func_4_GMM_new.R")
