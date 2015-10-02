#read bus stops data
bus_stops_data_filename <- paste0(csv_dir,"/../ParisData/RATPOpenData/Accessibility RATP bus stops/accessibilite_arrets_bus.csv")
bus_stops_data <- read.csv(bus_stops_data_filename, header = T) 

#remove NA entries
bus_stops_data <- bus_stops_data[which(!is.na(bus_stops_data$IDPTAR)),]
bus_stops_data <- bus_stops_data[order(bus_stops_data$CODE.INSEE),]

#keep INSEE codes 75101 to 75110 for districts 1 to 20 in paris
bus_stops_data <- subset(bus_stops_data, CODE.INSEE>=75101 & CODE.INSEE<=75110)

bus_stops_data_write_filename <- paste0(csv_dir,"/../ParisData/RATPOpenData/Accessibility RATP bus stops/bus_stops_data_dis1to10.csv")

write.csv(bus_stops_data, file=bus_stops_data_write_filename)

######reading converted data
bus_stops_data_converted_filename <- paste0(csv_dir,"/../ParisData/RATPOpenData/Accessibility RATP bus stops/bus_stops_data_dis1to10_converted.csv")
bus_stops_data_converted <- read.csv(bus_stops_data_converted_filename, header = T) 



