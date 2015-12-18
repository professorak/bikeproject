a <- merge(wdcMerged,as.ffdf(weather_data), by=c("day","time_halfhour_int"), all.x=T)
ffwhich(a,is.na(a$Conditions_num))[1]
a[ffwhich(a,is.na(a$Conditions_num))[1],]


dim(weather_data)
b <- weather_data[which(weather_data$day==161),]
c <- weather_data[which(weather_data$day>=150 & weather_data$day<=180),]

161, 2013-06-11
166, 2013-06-16
167, 2013-06-17
170 2013-06-20




a <- merge(wdcMerged,as.ffdf(weather_scale_df), by=c("day","time_halfhour_int"), all.x=T)    
ffwhich(a, is.na(outdem_weather_scale))[1]
a[ffwhich(a, is.na(outdem_weather_scale))[1],]
table(weather_scale_df$day)


a <- wdcMerged_weather_df[which(wdcMerged_weather_df$day_time_halfhour_int_fac %like% "161 "),]
