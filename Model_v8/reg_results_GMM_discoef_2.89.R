################################


#points$density = points$density/100
library("ipoptr")
wdcMerged$out_dem_sum[which(wdcMerged$out_dem_sum<=0.0001 & 
                              wdcMerged$stocked_out==FALSE)] <- 0.0001

wdcMerged <- wdcMerged[order(wdcMerged$tw_group,wdcMerged$station_id_index),]
current_serv_lvl <- current_serv_lvl[order(current_serv_lvl$tw_group,current_serv_lvl$station_id_index),]
user_serv_lvl <- current_serv_lvl
###
colnames_theta1 <<- c("dis coef","rand coef","density ridership","density metro","density intercept")
density_ridership_col <<- 3
density_metro_col <<- 4
density_intercept_col <<- 5

x0 <- c(-4.044482381,1.1,1.1,1.1)
if(!is.null(prev_theta)) {
  x0 <- prev_theta[-2]
}
lb = c(-20,0.5,1,0.1)
ub = c(0,1000000,10000000,10000000)

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

##market share rough
wdcMerged_nostckout <- wdcMerged[which(wdcMerged$stocked_out==F),]
st_demand <- by(wdcMerged_nostckout$out_dem_sum, wdcMerged_nostckout$station_id_index, FUN=sum)
st_weight <- by(wdcMerged_nostckout$obs_weight, wdcMerged_nostckout$station_id_index, FUN=sum)
st_demand <- st_demand/st_weight
sum(st_demand)
tot_density <- sum(points$weight[which(points$type==1)])*theta1[2] + 
  sum(points$weight[which(points$type==2)])*theta1[3] + length(which(points$type==1))*theta1[4]
print(paste0("market share: ", round(sum(st_demand)/tot_density*100,8),"%"))


# theta2 <- eval_error_xi(c(theta1[1],0),wdcMerged,points)$theta2
# delta_list <- compute_delta_list_cntrt_map_new(c(theta1[1],0),wdcMerged,points)
# 
# ddc  <- demand_distance_change_new(1.1,c(theta1[1],0),delta_list, wdcMerged,points)
# sdc <- demand_serv_lvl_change_new2(c(theta1[1],0),wdcMerged,points,0.9,
#                                    theta2[which(rownames(theta2)=="serv_lvl_covar")])
# ret_demand_lost_st <- demand_lost_st(c(theta1[1],0),wdcMerged,points,delta_list)
# perc_frac_dem_lost <- -(ret_demand_lost_st[[1]]/ret_demand_lost_st[[2]])*100
# subs_perc <- mean(perc_frac_dem_lost)
#   
# if(!exists("reg_runs_save_file")) {  
#     reg_runs_save_file <- "reg_runs_discoef_1.3.txt"
# }
# 
# write("tract: ",reg_runs_save_file, append=TRUE)
# write(paste0(tract,collapse="_") ,reg_runs_save_file, append=TRUE)
# write(theta1,reg_runs_save_file, append=TRUE)
# write("ddc",reg_runs_save_file, append=TRUE)
# write(ddc,reg_runs_save_file, append=TRUE)
# write("sdc",reg_runs_save_file, append=TRUE)
# write(sdc,reg_runs_save_file, append=TRUE)
# write("subs_perc",reg_runs_save_file, append=TRUE)
# write(subs_perc,reg_runs_save_file, append=TRUE)
# write("points density",reg_runs_save_file, append=TRUE)
# write(points$density[1],reg_runs_save_file, append=TRUE)
# tw_list_print <- paste0(unique(wdcMerged$tw), collapse="_")
# write("tw_list",reg_runs_save_file, append=TRUE)
# write(tw_list_print,reg_runs_save_file, append=TRUE)
# write("theta2",reg_runs_save_file, append=TRUE)
# write(theta2,reg_runs_save_file, append=TRUE)
# write("time",reg_runs_save_file, append=TRUE)
# write(time["elapsed"],reg_runs_save_file, append=TRUE)
# #write("st_err",reg_runs_save_file, append=TRUE)
# #write(st_err,reg_runs_save_file, append=TRUE)
# write("\n\n\n",reg_runs_save_file, append=TRUE)
# 
# 
