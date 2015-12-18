output_theta_csvformattted <- function(thetaden, thetanonden, thetatstat_den, thetatstat_nonden, theta2coef,
                                       theta2coeftstat, obj, itr, runtime) {
  outcoef <- cbind(thetaden, thetatstat_den)
  if(length(theta2coef)==0) {
    outcoef <- rbind(outcoef, matrix("",8,2))
  } else {
    outcoef <- rbind(outcoef, matrix("",2,2))
    outcoef <- rbind(outcoef, cbind(theta2coef,theta2coeftstat))
    outcoef <- rbind(outcoef, matrix("",8-2-length(theta2coef),2))
  }
  outcoef <- rbind(outcoef, cbind(thetanonden, thetatstat_nonden))
  outcoef <- rbind(outcoef, c("",""))
  outcoef <- rbind(outcoef, c(obj,""))
  outcoef <- rbind(outcoef, c(runtime,""))
  outcoef <- rbind(outcoef, c(itr,""))
  return(outcoef)
}


printoutput_theta_csvformattted <- function(thetaden_s1, thetanonden_s1, thetatstat_den_s1, thetatstat_nonden_s1, theta2coef_s1, 
  theta2coeftstat_s1, obj_s1, itr_s1, runtime_s1, thetaden, thetanonden, thetatstat_den, thetatstat_nonden, theta2coef,
  theta2coeftstat, obj, itr, runtime) {
  
  out_1 <- output_theta_csvformattted(thetaden_s1, thetanonden_s1, thetatstat_den_s1, thetatstat_nonden_s1, theta2coef_s1, 
                                      theta2coeftstat_s1, obj_s1, itr_s1, runtime_s1)
  out_2 <- output_theta_csvformattted(thetaden, thetanonden, thetatstat_den, thetatstat_nonden, theta2coef,
                                      theta2coeftstat, obj, itr, runtime)
  
  out_12 <- apply(cbind(out_1,out_2),1,FUN=function(x){paste0(x, collapse = ", ")})
  out_12 <- matrix(out_12,,1)
  print(out_12)  
  
  printdemand_division_df <- t(t(colSums(demand_division_df)))/sum(demand_division_df)*100
  printdemand_division_df <- cbind(rownames(printdemand_division_df),printdemand_division_df)
  rownames(printdemand_division_df) <- NULL
  printdemand_division_df <- apply(printdemand_division_df, 
                                   1,FUN=function(x){paste0(x, collapse = ", ")})
  printdemand_division_df <- matrix(printdemand_division_df,,1)
  print(printdemand_division_df)  
}







