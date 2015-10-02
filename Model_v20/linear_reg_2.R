#linear regression:
# regress wdcMerged$out_dem_mean with weights as wdcMerged$obs_weight on 
# 1. metro coefficients, on and off.
# metro_den_on_a,b,c,d,1,2 metro_den_off_a,b,c,d,1,2 
# metro_den_off_1 + metro_den_off_2 + 
#   metro_den_off_a + metro_den_off_b + metro_den_off_c + metro_den_off_d
# 2. census density
# 3. google places count
# 4. bus density
# 5. Later include tw fixed effect
# 6. Later include tract_tw fixed effect, while not including tract fixed effect

#First only regression 8_1.
#Then on all data

source("debug_functions.R")

###### Load saved file if exists
file <- "total_station_demand_linearreg"
filepath <- paste0(csv_dir,"/Paris/",file,".RData")
if(!file.exists(filepath)) {
  wdcMerged <- get_total_station_demand()  
  save(wdcMerged,file=filepath)
}
###
load(filepath)  

wdcMerged$metro_nearby <- (wdcMerged$metro_den_on_1>0)


# summary(lm(out_dem_mean ~ metro_den_on_1 + metro_den_on_2 + metro_den_on_a + metro_den_on_b + 
#      metro_den_on_c  + census_density + googleplaces_den_1 + googleplaces_den_2 +
#      bus_den_1 + bus_den_2, data=wdcMerged, weights=wdcMerged$obs_weight))


#All tws
library("AER")

# list_density_vars <- paste0("census_density + metro_nearby + metro_den_on_1 + metro_den_on_2",
#                             "+ catchment_area_step1 + catchment_area_step2",
#                             "+ catchment_area_step3 + catchment_area_step4", 
#                             "+ metro_den_off_1 + metro_den_off_2 + bus_den_1 + bus_den_2 + googleplaces_den_1 + googleplaces_den_2",
#                             "+ metro_den_on_3 + metro_den_on_4 + metro_den_off_3 + metro_den_off_4 + bus_den_3",
#                             "+ bus_den_4 + googleplaces_den_3 + googleplaces_den_4 + metro_den_on_a + metro_den_on_b + metro_den_on_c + ",
#                             "  metro_den_off_a + metro_den_off_b + metro_den_off_c + googleplaces_den_a + googleplaces_den_b",
#                             "+ googleplaces_den_c + googleplaces_food_den_1 +  googleplaces_food_den_2 +  googleplaces_food_den_3 +  googleplaces_food_den_4",
#                             "+ googleplaces_food_den_a +  googleplaces_food_den_b +  googleplaces_food_den_c + googleplaces_grocery_den_1 +  googleplaces_grocery_den_2",
#                             "+  googleplaces_grocery_den_3 +  googleplaces_grocery_den_4 + googleplaces_grocery_den_a +  googleplaces_grocery_den_b +  googleplaces_grocery_den_c",
#                             "+ googleplaces_government_den_1 +  googleplaces_government_den_2 +  googleplaces_government_den_3 +  googleplaces_government_den_4",
#                             "+ googleplaces_government_den_a +  googleplaces_government_den_b +  googleplaces_government_den_c")
# 

list_density_vars <- paste0("census_density + metro_nearby + metro_den_on_1 + metro_den_on_2",
                            "+ catchment_area_step2 + census_density:catchment_area_step2", 
                            "+ metro_den_off_1 + metro_den_off_2 + bus_den_1 + bus_den_2 + googleplaces_den_1 + googleplaces_den_2",
                            "+ metro_den_on_3 + metro_den_on_4 + metro_den_off_3 + metro_den_off_4 + bus_den_3",
                            "+ bus_den_4 + googleplaces_den_3 + googleplaces_den_4 + metro_den_on_a + metro_den_on_b + metro_den_on_c + ",
                            "  metro_den_off_a + metro_den_off_b + metro_den_off_c + googleplaces_den_a + googleplaces_den_b",
                            "+ googleplaces_den_c + googleplaces_food_den_1 +  googleplaces_food_den_2 +  googleplaces_food_den_3 +  googleplaces_food_den_4",
                            "+ googleplaces_food_den_a +  googleplaces_food_den_b +  googleplaces_food_den_c + googleplaces_grocery_den_1 +  googleplaces_grocery_den_2",
                            "+  googleplaces_grocery_den_3 +  googleplaces_grocery_den_4 + googleplaces_grocery_den_a +  googleplaces_grocery_den_b +  googleplaces_grocery_den_c",
                            "+ googleplaces_government_den_1 +  googleplaces_government_den_2 +  googleplaces_government_den_3 +  googleplaces_government_den_4",
                            "+ googleplaces_government_den_a +  googleplaces_government_den_b +  googleplaces_government_den_c")


reg_formula <- as.formula(paste0("log(out_dem_mean+1) ~ ",list_density_vars))

summary(lm(reg_formula, data=wdcMerged))

summary(lm(out_dem_mean ~ census_density:catchment_area_step4, data=wdcMerged))

reg_formula_tractfe <- as.formula(paste0("log(out_dem_mean+1) ~ ",list_density_vars, " + tract_tw_fac "))


summary(ivreg(reg_formula_tractfe, data=wdcMerged, weights=wdcMerged$obs_weight))

#minimal vars
list_density_vars_minimal <- paste0("metro_den_on_1 + metro_den_off_1 + bus_den_1  + googleplaces_den_1 + ",
                                    "googleplaces_food_den_1 + googleplaces_grocery_den_1 + googleplaces_government_den_1")

reg_formula <- as.formula(paste0("log(out_dem_mean+1) ~ ",list_density_vars_minimal, "+ serv_lvl | ", list_density_vars_minimal,
                                 "+ instr_serv_lvl + serv_lvl_neighbours"))

reg_formula_tractfe <- as.formula(paste0("log(out_dem_mean+1) ~ ",list_density_vars_minimal, "+ serv_lvl + tract_tw_fac | ", list_density_vars_minimal,
                                         "+ instr_serv_lvl + serv_lvl_neighbours + tract_tw_fac"))

summary(ivreg(reg_formula, data=wdcMerged, weights=wdcMerged$obs_weight))

summary(ivreg(reg_formula_tractfe, data=wdcMerged, weights=wdcMerged$obs_weight))


