#get paris weather data for 2014
year <- 2013
months <- c(1:12)
month_days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
for(i in c(1:12)) {
  no_days <- month_days[i]
  for(day in 1:no_days) {
    date <- paste0(year,"-",months[i],"-",day)
    print(date)
    try({
      url <- paste0("http://www.weatherbase.com/weather/weatherhourly.php3?s=94170&cityname=Paris-%CEle-de-France-France&date=",
                  date,"&units=metric")
      content <- readLines(url)
      filename <- paste0("WeatherData/data_",date,".html")
      write(content,filename)
    })      
  }
}
  



