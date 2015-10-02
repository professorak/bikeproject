source("debug_functions.R")
theta1 <- get_theta1fromtheta(theta)
theta_save <- theta
theta_dis <- theta
theta_dis[-c(1:2)] <- 0
params_save <- params

tot_lambda <- c()
for(i in 3:length(theta)) {
  theta_i <- theta_dis
  theta_i[i] <- theta[i]
  theta1_i <- get_theta1fromtheta(theta_i)
  tot_lambda <- c(tot_lambda, 
    sum(eval_lambda_full(deltain_stin, theta1_i, wdcMerged, points)*wdcMerged$obs_weight))
  
}


paste(tot_lambda/sum(tot_lambda)*100, collapse = ",")


#what does the coefficients imply?
mean_delta <- mean(deltain_stin)
step_dis <- 116.7 #mts

#function to return demand at distance x in mts
demand_distancex <- function(x, step_dis) {
  
  util <- theta[1]*min(x,step_dis)/1000 + theta[2]*max(0,x-step_dis)/1000
  return(exp(util)/(1+exp(util)))
}

x_vec <- c(0:60)*10 #0 to 600mts at steps of 10

demand_distancex_vec <- c()
for(i in 1:length(x_vec)) {
  demand_distancex_vec <- c(demand_distancex_vec,
                            demand_distancex(x_vec[i],step_dis))
}

plot(x_vec, demand_distancex_vec)

