library(rjson)
library(RCurl)



metro_list_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/ParisMetroData/metro_stations_list.csv") 
metro_list <- read.csv(metro_list_file_name)

#only keeping paris stations
metro_list <- metro_list[!is.na(metro_list$paris_district),]

#crawl
metro_list_results <- c()
for(i in 1:nrow(metro_list)) {
  print(i)
  urlname <- paste0("https://maps.googleapis.com/maps/api/geocode/json?address=",
                gsub(" ", "+", metro_list$station_name[i]),"+Ile-de-France+France&key=AIzaSyA7iRBhzklZhLWqpvUTVkgxptSfrH1wbVQ")
  
  results <- fromJSON(paste(getURL(urlname), collapse="")[1])$results[[1]]
  if(!("train_station" %in% results$types | "subway_station" %in% results$types | 
       "transit_station" %in% results$types)) {
    urlname <<- paste0("https://maps.googleapis.com/maps/api/geocode/json?address=",
                  gsub(" ", "+", metro_list$station_name[i]),
                  "+",metro_list$station_type[i],
                  "+Ile-de-France+France&key=AIzaSyA7iRBhzklZhLWqpvUTVkgxptSfrH1wbVQ")
    
    results_1 <<- fromJSON(paste(getURL(urlname), collapse="")[1])$results[[1]]
    if("train_station" %in% results_1$types | "subway_station" %in% results_1$types | 
           "transit_station" %in% results_1$types) {
      results <<- results_1 
    }
  }
  metro_list_results <- rbind(metro_list_results,
    c(results$geometry$location$lat,results$geometry$location$lng,
    results$geometry$location_type,paste0(results$types,collapse=";"),
    results$formatted_address,as.numeric("train_station" %in% results$types),
    as.numeric("subway_station" %in% results$types),as.numeric("transit_station" %in% results$types)
    ))
}

metro_list_results_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/ParisMetroData/metro_stations_results_list.csv") 
colnames(metro_list_results) <- c("lat","lon","result_type","locations_types","formatted_address",
                                  "train_station","subway_station","transit_station")

metro_list_results <- cbind(metro_list,metro_list_results)

write.csv(metro_list_results,file=metro_list_results_file_name)



