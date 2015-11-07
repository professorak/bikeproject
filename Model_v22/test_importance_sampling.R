load(file="v0_vec.RData")
v0_vec <- rep(1, 160)
v0_vec_weights <- rep(1, length(v0_vec))

system.time({
  a <- eval_lambdasum_v0_vec(deltain_stin, theta1, wdcMerged, points)
})

system.time({
  eval_lambda_full(deltain_stin, theta1, wdcMerged, points)
})


system.time({
  i <- 3
  random_coef_acceptance_prob (v0_vec[i],
                                 deltain_stin, theta1, wdcMerged, points)     
  
})

a[1]/a[3]


source("Importance_sampling_functions.R")
v0_vec <- rnorm(20)
random_coef_acceptance_prob (v0_vec,
                             deltain_stin, theta1, wdcMerged, points)     
a <- generate_v0_importance_sampling(40, deltain_stin, theta1, wdcMerged, points)



