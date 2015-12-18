#Counterfactual Level 2
#Trade-off Curve
#read results_counterfactual_l2.csv
results_counterfactual_l2_tab <- read.csv("PaperPlots/results_counterfactual_l2.csv",header = T)
results_counterfactual_l2_tab <- results_counterfactual_l2_tab[order(1/results_counterfactual_l2_tab$scale_distance_factor_vec),]
results_counterfactual_l2_tab <- results_counterfactual_l2_tab[which(results_counterfactual_l2_tab$beta_list>0.25),]
scale_distance_inv <- 1/results_counterfactual_l2_tab$scale_distance_factor_vec
serv_lvl <- results_counterfactual_l2_tab$serv_lvl_factor_vec *
  results_counterfactual_l2_tab$base_serv_lvl[1]

lo <- lm(serv_lvl ~ 
           poly(scale_distance_inv,8))


pdf("PaperPlots/TradeoffCurve.pdf",width=3.25,height=3.25,
    pointsize=9*0.7)
plot(scale_distance_inv,
     serv_lvl, col="white", bty="n", xlab="Station Density", xaxt="n", ylab="Average Station Bike Availability")
lines(scale_distance_inv,predict(lo),col="black",lwd=1)
#station_density_name <- round((1/scale_distance_inv)^2, digits=2)
station_density_name <- c(0.25,1,2,4,8,16,32)
station_density_name_pos <- (station_density_name)^(-0.5)
axis(1,at=station_density_name_pos,labels=station_density_name)
dev.off()

####################################
#Counterfactual System-use curve
results_counterfactual_l2_tab <- read.csv("PaperPlots/results_counterfactual_l2.csv",header = T)
results_counterfactual_l2_tab <- results_counterfactual_l2_tab[order(1/results_counterfactual_l2_tab$scale_distance_factor_vec),]
results_counterfactual_l2_tab <- results_counterfactual_l2_tab[which(results_counterfactual_l2_tab$beta_list>0.25),]
scale_distance_inv <- 1/results_counterfactual_l2_tab$scale_distance_factor_vec

relative_system_use_change <- (results_counterfactual_l2_tab$demand_vec/ 
  results_counterfactual_l2_tab$demand_vec[which(results_counterfactual_l2_tab$beta_list==1)] - 1)*100

lo <- lm(relative_system_use_change ~ 
           poly(scale_distance_inv,5))


pdf("PaperPlots/CounterfactualL2Curve.pdf",width=3.25,height=3.25,
    pointsize=9*0.7)
plot(scale_distance_inv,
     relative_system_use_change, col="white", bty="n", xlab="Station Density", xaxt="n", ylab="Relative System-Use Change (%)")
lines(scale_distance_inv,predict(lo),col="black",lwd=1)
#station_density_name <- round((1/scale_distance_inv)^2, digits=2)
station_density_name <- c(0.25,1,2,4,8,16,32)
station_density_name_pos <- (station_density_name)^(-0.5)
axis(1,at=station_density_name_pos,labels=station_density_name)
dev.off()




