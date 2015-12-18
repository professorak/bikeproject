fit <- lm( log(wdcMerged$out_dem_mean+0.00001)[stocked_list] ~ 
             X[,"Intercept"] +  X[,other_vars] + 
             X[,tract_tw_vars] + X[,distance_vars] + X[,density_vars]
           + Instr_df[,c("serv_lvl","serv_lvl_sq")] + 0, 
           weights=wdcMerged$obs_weight[stocked_list])





x <- c("A and B", "A, B and C", "A, B, C and D", "foobar")
pattern <- "[[:space:]]*(,|and)[[:space:]]"
## Match data from regexpr()
m <- regexpr(pattern, x)
regmatches(x, m)

x <- "tract_tw_fac17_4"
#pattern <- "tract_tw_fac(\\d+)*"



