#read paris metro locations data

locations_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/Metropolitain.io/out_stations.csv") 
metro_locations <- read.csv(file=locations_file_name,fileEncoding="latin1")
metro_locations$Station <- as.character(metro_locations$Station)
metro_locations <- 
  metro_locations[which(!duplicated(metro_locations[,])),]
metro_locations <- subset(metro_locations, !(Station %in% c("COURONNES","COURCELLES","ARTS ET METIERS")))
additional_locations <- rbind(
  c(NA,48.879268,2.303304, 'COURCELLES', NA, NA),
  c(NA,48.869274, 2.380155, 'COURONNES' , NA, NA),          
  c(NA,48.768713, 2.464298, 'CRETEIL-POINTE DU LAC'  , NA, NA),
  c(NA,48.906675, 2.365796, 'FRONT POPULAIRE', NA, NA),   
  c(NA,48.885227, 2.342544, 'FUNICULAIRE', NA, NA),  
  c(NA, 48.93016, 2.284074, 'LES COURTILLES'  , NA, NA),
  c(NA, 48.86529, 2.356378, 'ARTS ET METIERS', NA, NA), 
  c(NA,48.817889, 2.31922, 'MAIRIE DE MONTROUGE', NA, NA)
)
colnames(additional_locations) <- colnames(metro_locations)
metro_locations <- rbind(metro_locations, additional_locations)
metro_locations$Station <- trim(metro_locations$Station) 
metro_locations <- metro_locations[order(metro_locations$Station),]


#read paris metro incoming traffic data
metrotraffic_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/RATPOpenData/AnnualIncomingTrafficPerStation2014/trafic_entrant2014.csv") 
metro_traffic <- read.csv(file=metrotraffic_file_name,fileEncoding="latin1")
metro_traffic$Station <- as.character(metro_traffic$Station)
metro_traffic$Traffic <- as.numeric(gsub(",","", metro_traffic$Traffic))
metro_traffic$Station <- trim(metro_traffic$Station) 
metro_traffic <- metro_traffic[order(metro_traffic$Station),]

#merge
metro_traffic_locations <- merge(metro_traffic, metro_locations, all.x=T, by="Station")
metro_traffic_locations$lat <- as.numeric(metro_traffic_locations$lat)
metro_traffic_locations$lon <- as.numeric(metro_traffic_locations$lon)
metro_traffic_locations <- metro_traffic_locations[,c('Station','Traffic','lat','lon')]
#missing stations location data
metrotrafficloc_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/ParisTransportData/metro_traffic_locations.csv") 
write.csv(metro_traffic_locations, file=metrotrafficloc_file_name)




