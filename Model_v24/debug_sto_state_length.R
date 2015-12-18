source("GetDataGooglePlaces_aggmonths_averaged_20districts_norealloc_stkoutthresh_five.R")

length_sto_state <- function(str,splitchar="_"){  
  return(length(strsplit(as.character(str)
                                     , split=splitchar)[[1]]))
}

length_stostate_vec <- tapply(wdcMerged$sto_state_local,c(1:nrow(wdcMerged)),FUN=length_sto_state)
summary(length_stostate_vec)
