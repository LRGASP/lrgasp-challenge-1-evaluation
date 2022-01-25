#####################################
##### SQANTI3 report generation ######
#####################################



### Francisco J Pardo-Palacios
### Last Modified: 05/03/2020 by francisco.pardo.palacios@gmail.com

#********************** Taking arguments from python script

args <- commandArgs(trailingOnly = TRUE)
class.file <- args[1]
junc.file <- args[2]
name <- args[3]
utilities.path <- args[4]
platform <- args[5]
rdata <- args[6]
if (length(args)>6){
  sim_prefix <- args[7]
}

report.prefix <- strsplit(class.file, "_classification.txt")[[1]][1];
report.file <- paste(report.prefix, "Evaluation_report.html", sep="_");


#********************** Packages (install if not found)

list_of_packages <- c("ggplot2", "scales", "knitr","rmarkdown", "data.table")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

library(ggplot2)
library(scales)
library(knitr)
library(rmarkdown)
library(data.table)
library(tidyverse)
#********************* Run Calculation scripts

setwd(utilities.path)

if (exists("sim_prefix")){
  source("LRGASP_calculations.simulation.R")
  
  LRGASP_calculations(NAME = name , out.dir = rdata,
                      class.file=class.file, junc.file=junc.file,
                      functions.dir = utilities.path, sim_prefix = sim_prefix)
  
  RMD = paste(utilities.path, "Evaluation_metrics.simulation.Rmd", sep = "/")
}else{
  source("LRGASP_calculations.R")
  
  LRGASP_calculations(NAME = name , out.dir = rdata,
                      class.file=class.file, junc.file=junc.file,
                      functions.dir = utilities.path)
  
  RMD = paste(utilities.path, "Evaluation_metrics.Rmd", sep = "/")
}


rmarkdown::render(RMD, params = list(
  output.directory = rdata,
  Name = name,
  Platform = platform  ), output_file = report.file
)

