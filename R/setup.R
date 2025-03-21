
# load required packages ----
library(devtools)
if(identical(find.package("splott", quiet=TRUE), character(0))) {
  # if splott not installed, install it
  install_github("paigejo/splott")
}
library(splott)
library(fields)
library(INLA)
library(inlabru)
library(akima)
library(ks)

# set R and package options ----
inla.setOption(num.threads=1) # consider raising
options(error=recover)

# set working directory ----
setwd("~/git/PrefDrill/")

# source R code ----
source('R/modSPDE.R')
source('R/modSuccessiveSPDE.R')
source('R/modWellDensity.R')
source('R/utilityFuns.R')
source('R/simStudy.R')
source('R/getTestDat.R')

# set up global variables ----

# directories
rDir <<- "~/git/PrefDrill/R/" # R functions, setup.R script
figDir <<- "~/git/PrefDrill/figures/" # figures
dataDir <<- "~/git/PrefDrill/data/" # for storing real world datasets
outputDir <<- "~/git/PrefDrill/savedOutput/" # for storing output from the R code
globalDir <<- "~/git/PrefDrill/savedOutput/global/" # stores datasets used globally
scriptDir <<- "~/git/PrefDrill/scripts"

# other variables
simStudyXlims <<- c(0, 10000) # TODO: CHECK THIS
simStudyYlims <<- c(0, 10000) # TODO: CHECK THIS


