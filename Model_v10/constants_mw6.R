#constants file

rm(list=ls())

if(.Platform$OS.type=="windows"  & Sys.info()["user"]=="Ashish") {
  csv_dir <- "C:\\Users\\Ashish\\Desktop\\Velib Data\\CSVData" 
  work_dir <- "C:/Users/Ashish/Dropbox/Ashish_Karan/Sustainable Transportation System/Data/Code/R/Model_v10"  
  dropbox_dir <- "C:/Users/Ashish/Dropbox/Ashish_Karan"
} else {
  if(.Platform$OS.type=="windows" & Sys.info()["user"]=="A") {
    #incidates my XPS DELL Machine
    work_dir <- "C:/Users/A/Dropbox/Ashish_Karan/Sustainable Transportation System/Data/Code/R/Model_v10"
    csv_dir <- "C:/Users/A/Dropbox/VelibData/CSVData"
    dropbox_dir <- "C:/Users/A/Dropbox/Ashish_Karan"   
  } 
  else {
  if((Sys.info()["user"])=="ashish") {
    #incidates my linux machine
    work_dir <- "/home/ashish/Dropbox/Ashish_Karan/Sustainable Transportation System/Data/Code/R/Model_v10"
    csv_dir <- "/home/ashish/Dropbox/VelibData/CSVData"
    dropbox_dir <- "/home/ashish/Dropbox/Ashish_Karan"    
    #options(fftempdir="/home/ashish/fftempdir")
  }
  if(.Platform$OS.type=="unix"  & (Sys.info()["user"])=="ashishkabra") {
    #incidates my macbook
    work_dir <- "/Users/ashishkabra/Dropbox/Ashish_Karan/Sustainable Transportation System/Data/Code/R/Model_v10"
    csv_dir <- "/Users/ashishkabra/Dropbox/VelibData/CSVData"
    dropbox_dir <- "/Users/ashishkabra/Dropbox/Ashish_Karan"    
    #options(fftempdir="/home/ashish/fftempdir")
  }
  if((Sys.info()["user"])=="ak") {
    #incidates my linux machine
    work_dir <- "/home/ak/Dropbox/Ashish_Karan/Sustainable Transportation System/Data/Code/R/Model_v10"
    csv_dir <- "/home/ak/Desktop/Velib_Data/CSVData"
    dropbox_dir <- "/home/ak/Dropbox/Ashish_Karan"    
  }
  if((Sys.info()["user"])=="rstudio") {
    csv_dir <- "/home/rstudio/workspace/CSVData"
    work_dir <- "/home/rstudio/Data/Model_v10"
    dropbox_dir <- "/home/ubuntu/Dropbox/Ashish_Karan"    
  }
}
}
setwd(work_dir)

st_id_file <- "Paris/st_id_data_01_07_13.txt"

max_walking_dis <- 0.6

library('nloptr')
library('ffbase')
library('biglm')
library('sfsmisc')

library("Rcpp")
library("RcppArmadillo")
library("ipoptr")
library("plm")
library("sandwich")


source("setup.R")

options(digits=10)
#for ffdf chunk size
options(ffbatchbytes=167772160)

require("colorspace")


#obj.size <- napply(names, function(x) round(object.size(x)/1024^2,2))
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  #obj.size <- napply(names, object.size)
  obj.size <- napply(names, function(x) round(object.size(x)/1024^2,2))
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}
# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
