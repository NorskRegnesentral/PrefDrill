
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
library(gstat)
library(sf)
library(Matrix)
library(terra)

# set R and package options ----
inla.setOption(num.threads=1) # consider raising
options(error=recover)

# set working directory ----
setwd("~/git/PrefDrill/")

# source R code ----
source('R/modSPDE.R')
source('R/modSuccessiveSPDE.R')
source('R/modWatson.R')
source('R/modWellDensity.R')
source('R/utilityFuns.R')
source('R/simStudy.R')
source('R/getTestDat.R')
source("R/readRMS.R")

# set up global variables ----

# directories
rDir <<- "~/git/PrefDrill/R/" # R functions, setup.R script
figDir <<- "~/git/PrefDrill/figures/" # figures
dataDir <<- "~/git/PrefDrill/data/" # for storing real world datasets
outputDir <<- "~/git/PrefDrill/savedOutput/" # for storing output from the R code
globalDir <<- "~/git/PrefDrill/savedOutput/global/" # stores datasets used globally
scriptDir <<- "~/git/PrefDrill/scripts/"

# globals folder
out = load(paste0(globalDir, "testDat_truthPref.RData")) # "truthTestDat", "seismicTestDat", "wellTestDat" 
truthTestDat_truthPref = truthTestDat
seismicTestDat_truthPref = seismicTestDat
wellTestDat_truthPref = wellTestDat
out = load(paste0(globalDir, "testSuccessiveDat_ipp_uniform.RData")) # "truthTestDat", "seismicTestDat", "wellTestDat" 
truthTestDat_successiveIPP = truthTestDat
seismicTestDat_successiveIPP = seismicTestDat
wellTestDat_successiveIPP = wellTestDat
out = load(paste0(globalDir, "testDat.RData")) # "truthTestDat", "seismicTestDat", "wellTestDat" 

# other variables
simStudyXlims <<- c(-12.5, 15012.5)
simStudyYlims <<- c(-8.3335, 5008.4336)

# simStudyXlims <<- c(0, 10000) # TODO: CHECK THIS
# simStudyYlims <<- c(0, 10000) # TODO: CHECK THIS



