################################


#points$density = points$density/100
wdcMerged$out_dem_sum[which(wdcMerged$out_dem_sum<=0.0001 & 
                              wdcMerged$stocked_out==FALSE)] <- 0.0001

wdcMerged <- wdcMerged[order(wdcMerged$tw_group,wdcMerged$station_id_index,wdcMerged$sto_state_local),]
current_serv_lvl <- current_serv_lvl[order(current_serv_lvl$tw_group,current_serv_lvl$station_id_index),]
user_serv_lvl <- current_serv_lvl
###
#x0 <- c(-10,-1)
x0 <- c(-4,  -11)
if(!is.null(prev_theta)) {
  x0 <- prev_theta[1,3]
}
lb = c(-40,-40)
ub = c(40,100)

x0_start <<- NULL
prev_theta <<- NULL
W_optimal <<- NULL
#W_optimal <<- eval_weighting_mat_new(c(-10),wdcMerged,points) 

opts <- list("tol"=1.0e-4,
             #"print_options_documentation"='yes'
             "print_info_string"='yes'
             #,"print_level"=12
)
#W_optimal <- diag(c(1,1,1))
#covariates_st <- gen_covariates(wdcMerged)
time <- proc.time()
res <- ipoptr( x0=x0,
               lb=lb,
               ub=ub,
               eval_f=eval_obj_GMM_list_extended_new_obj, 
               eval_grad_f=eval_obj_GMM_list_extended_new_grad,
               wdcMerged=wdcMerged, 
               points=points, 
               opts=opts
) 
theta1 <- res$solution 
print(theta1)
time <- proc.time() - time
theta1 <- c(theta1[1],0,theta1[2])
theta2 <- eval_error_xi(theta1,wdcMerged,points)$theta2
if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
  delta_list <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
  prev_theta <<- theta1
  prev_deltain <<- deltain
} else {
  delta_list <- prev_deltain
}
ddc  <- demand_distance_change_new_nls(1.1,theta1,delta_list, wdcMerged,points)

sdc <- demand_serv_lvl_change_new2(theta1,wdcMerged,points,0.9,
                                   theta2[which(rownames(theta2) == "serv_lvl_covar")])
subs_perc <- average_demand_susbtitution_demModel(theta1,wdcMerged,points,delta_list)
# dis_scaling <- c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.20, 2.00, 4.00)

# serv_lvl_values <- c(0.499578416947, 0.667279607877, 0.750356654283, 0.798240923192, 0.832845203737, 0.857831671781, 0.889291282838, 0.899833962270, 0.917188255333, 0.929012726490, 0.936217702815, 0.944900500743, 0.951098883953, 0.957015177987, 0.960315335619, 0.964528097363, 0.967168731425, 0.976991182653, 0.991486347336, 0.997805086199)
# serv_lvl_base <- serv_lvl_values[which(dis_scaling==1)]
# ddc_vec <- c()
# sdc_vec <- c()

# for(i in 1:length(dis_scaling)) {
	# ddc_vec <- c(ddc_vec, demand_distance_change_new_nls(1/dis_scaling[i],theta1,delta_list, wdcMerged,points)[1])
	# sdc_vec <- c(sdc_vec,demand_serv_lvl_change_new2(theta1,wdcMerged,points,serv_lvl_values[i]/serv_lvl_base,theta2[which(rownames(theta2) == "serv_lvl_covar")])[1])
# }


#var_covar <- eval_covariace_new (theta1, wdcMerged, points)
#st_err <- c(sqrt(var_covar[1,1]),sqrt(var_covar["serv_lvl_covar","serv_lvl_covar"]))  
  
if(!exists("reg_runs_save_file")) {  
    reg_runs_save_file <- "reg_runs_discoef_1.3.txt"
}

write("tract: ",reg_runs_save_file, append=TRUE)
write(paste0(tract,collapse="_") ,reg_runs_save_file, append=TRUE)
write(theta1,reg_runs_save_file, append=TRUE)
write("ddc",reg_runs_save_file, append=TRUE)
write(ddc,reg_runs_save_file, append=TRUE)
write("sdc",reg_runs_save_file, append=TRUE)
write(sdc,reg_runs_save_file, append=TRUE)
write("subs_perc",reg_runs_save_file, append=TRUE)
write(subs_perc,reg_runs_save_file, append=TRUE)
write("points density",reg_runs_save_file, append=TRUE)
write(points$density[1],reg_runs_save_file, append=TRUE)
tw_list_print <- paste0(unique(wdcMerged$tw), collapse="_")
write("tw_list",reg_runs_save_file, append=TRUE)
write(tw_list_print,reg_runs_save_file, append=TRUE)
write("theta2",reg_runs_save_file, append=TRUE)
write(theta2,reg_runs_save_file, append=TRUE)
write("time",reg_runs_save_file, append=TRUE)
write(time["elapsed"],reg_runs_save_file, append=TRUE)
# write("ddc_vec",reg_runs_save_file, append=TRUE)
# write(ddc_vec,reg_runs_save_file, append=TRUE)
# write("sdc_vec",reg_runs_save_file, append=TRUE)
# write(sdc_vec,reg_runs_save_file, append=TRUE)
#write("st_err",reg_runs_save_file, append=TRUE)
#write(st_err,reg_runs_save_file, append=TRUE)
write("\n\n\n",reg_runs_save_file, append=TRUE)


