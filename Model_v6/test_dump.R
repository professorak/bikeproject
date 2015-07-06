path <- "I:/ET_Data/All\ Rides\ All\ Countries"
filelist <- c("temp_AllRidesAugust2014.csv", "temp_AllRides_Apr2014.csv", "temp_AllRides_Apr2015.csv", 
"temp_AllRides_Before2014.csv", "temp_AllRides_Dec2014.csv", "temp_AllRides_Feb2014.csv", 
"temp_AllRides_Feb2015.csv", "temp_AllRides_Jan2014.csv", "temp_AllRides_Jan2015.csv", 
"temp_AllRides_July2014.csv", "temp_AllRides_Jun2014.csv", "temp_AllRides_Jun_1_12_2015.csv",
"temp_AllRides_Mar2014.csv", "temp_AllRides_Mar2015.csv", "temp_AllRides_May2015.csv", 
"temp_AllRides_Nov2014.csv", "temp_AllRides_Oct2014.csv", "temp_AllRides_Sep2014.csv",
"temp_May2014.csv")

for(i in c(1:length(filelist))) {
  filename <- paste0(path,"/",filelist[i])
  a <- read.csv(filename)
  date <- as.POSIXct(a$requested_at, format="%Y-%m-%d %H:%M:%S", tz="Asia/Singapore")
  day <- as.POSIXlt(date)$mday
  mon <- as.POSIXlt(date)$mon
  mon_day <- paste0(mon,"-",day)
  print(filelist[i])
  print(unique(mon_day[order(day)]))  
}

