a <- read.csv(file=filename,head=FALSE,sep=",",
              nrows=1000000)
a <- subset(a,V3!=-1)
a[which(a$V1==1001)[1],]


setClass("myDate")
setAs("character","myDate", function(from) as.POSIXct(from, format='%Y-%m-%d %H:%M:%S',tz="Europe/London") )

a <- read.csv.ffdf(file=filename,head=FALSE,sep=",",
                        nrows=10,colClasses=c(V1="integer",V2="myDate",
                                              V3="integer",V4="integer",V5="integer"))  

colnames(a) <- c("station_id" ,"time" ,"bikes" ,"spaces_available" ,"unbalanced")
a$time[1]


as.POSIXlt(wdcRawOrg$time[1],"%Y-%m-%d %H:%M:%S",
           tz="Europe/London")
wdcRawOrg$time[1]

as.POSIXlt(wdcRawOrg$time[1],"%Y-%m-%d %H:%M:%S",
           tz="Europe/Paris")$hour
a_BST <- as.POSIXlt(wdcRawOrg$time[1],"%Y-%m-%d %H:%M:%S",
                    tz="Europe/London")
a_CEST  <- as.POSIXlt(wdcRawOrg$time[1],"%Y-%m-%d %H:%M:%S",
                      tz="Europe/Paris")
as.POSIXlt(a_BST)$hour
as.POSIXlt(a_CEST)$hour
as.POSIXlt(a_BST)$yday
as.POSIXlt(a_CEST)$yday





setClass("myDate")
setAs("character","myDate", function(from) as.POSIXlt(as.POSIXct(from, format='%Y-%m-%d %H:%M:%S',tz="Europe/London"), 
                                                      format='%Y-%m-%d %H:%M:%S', tz="Europe/Paris") )

wdcRaw <- read.csv.ffdf(file=filename,head=FALSE,sep=",",
                        nrows=100,colClasses=c(V1="integer",V2="myDate",
                                              V3="integer",V4="integer",V5="integer"))  
colnames(wdcRaw) <- c("station_id" ,"time" ,"bikes" ,"spaces_available" ,"unbalanced")


a_CEST_format <- format(as.character(a_BST), tz="Europe/Paris", usetz=T)

as.ff(as.POSIXlt(wdcRawOrg$time[])$hour)

attr(a_BST, "tzone") <- "Europe/Paris"


as.POSIXlt(as.POSIXct(from, format='%Y-%m-%d %H:%M:%S',tz="Europe/London"), format='%Y-%m-%d %H:%M:%S', tz="Europe/Paris")


a <- wdcRawOrg$time[1]
a
attr(a, "tzone") <- "Europe/London"

a


attr(wdcRawOrg_all$time[1], "tzone") <-  "Europe/Paris"

