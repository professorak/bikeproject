dis <- 0.010
  
theta_dis <- c(-2.210, -19.652)  

#theta_dis <- c(-0.251, -50.376)
100*(1/2-(exp(theta_dis[1]*dis + theta_dis[2]*dis^2)/(1+exp(theta_dis[1]*dis + theta_dis[2]*dis^2))))/(1/2)

100*(1-exp(theta_dis[1]*dis + theta_dis[2]*dis^2))

