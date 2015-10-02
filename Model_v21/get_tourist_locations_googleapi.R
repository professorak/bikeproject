library(rjson)
library(RCurl)



tourist_list_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/ParisTourismeData/tourist_locations_list.csv") 
tourist_list <- read.csv(tourist_list_file_name)

#crawl
tourist_list_results <- c()
for(i in 1:nrow(tourist_list)) {
  print(i)
  urlname <- paste0("https://maps.googleapis.com/maps/api/geocode/json?address=",
                    gsub(" ", "+", tourist_list$name[i]),"+Paris+France&key=AIzaSyA7iRBhzklZhLWqpvUTVkgxptSfrH1wbVQ")
  
  results <- fromJSON(paste(getURL(urlname), collapse="")[1])$results[[1]]
  tourist_list_results <- rbind(tourist_list_results,
                              c(results$geometry$location$lat,results$geometry$location$lng,
                                results$geometry$location_type,paste0(results$types,collapse=";"),
                                results$formatted_address))
}

tourist_list_results_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/ParisTourismeData/tourist_results_list.csv") 
colnames(tourist_list_results) <- c("lat","lon","result_type","locations_types","formatted_address")

tourist_list_results <- cbind(tourist_list,tourist_list_results)

write.csv(tourist_list_results,file=tourist_list_results_file_name)



