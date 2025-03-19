
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

# set R and package options ----
inla.setOption(num.threads=1) # consider raising
options(error=recover)

# set working directory ----
setwd("~/git/PrefDrill/")

# source R code ----
source('R/modSPDE.R')
source('R/utilityFuns.R')

# set up global variables ----

# directories
rDir <<- "~/git/PrefDrill/R/" # R functions, setup.R script
figDir <<- "~/git/PrefDrill/figures/" # figures
dataDir <<- "~/git/PrefDrill/data/" # for storing real world datasets
outputDir <<- "~/git/PrefDrill/savedOutput/" # for storing output from the R code
globalDir <<- "~/git/PrefDrill/savedOutput/global/" # stores datasets used globally
scriptDir <<- "~/git/PrefDrill/scripts"

# other variables
simStudyXlims <<- c(0, 1000) # TODO: CHECK THIS
simStudyYlims <<- c(0, 1000) # TODO: CHECK THIS


