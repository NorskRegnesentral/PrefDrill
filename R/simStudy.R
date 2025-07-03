# functions for the simulation study

# construct a rectangular prediction grid over a domain
# xLims: a 2-vector with the min and max x value
# yLims: a 2-vector with the min and max y value
getSimStudyPredGrid = function(xLims = simStudyXlims, yLims = simStudyXlims) {
  make.surface.grid(list(x=xLims, y=yLims))
}

getFitModFuns = function() {
  funs = list()
  
  # basic SPDE (no seismic)
  # funs = c(funs, list(function(...) {fitSPDEsimDat(addKDE=TRUE, ...)}))
  
  # SPDE with seismic
  funs = c(funs, list(fitSPDEsimDat))
  
  # SPDE with kde covariate
  funs = c(funs, list(function(...) {fitSPDEsimDat(addKDE=TRUE, ...)}))
  
  # Diggle model
  funs = c(funs, list(fitDigglesimDat))
  
  # Watson et al. model
  funs = c(funs, list(fitWatsonSimDat))
  
  funs
}


# proportion of domain area -> repel radius bandwidth
repAreaToDist = function(repArea=.01) {
  # A = pi R^2
  # R = sqrt(A/pi)
  domainArea = diff(simStudyXlims) * diff(simStudyYlims)
  sqrt(repArea * domainArea / pi)
}

# samples the next drilling location.
# 
# If betaP == gammaP == 0, sampling is uniform, 
# otherwise with prob proportional to the exponential of seismic and spatial data 
# times betaP and gammaP respectively. Regardless, samples uniformly within 
# selected grid cell. Calculates seismic estimate at sampled location based on 
# bilinear interpolation.
# 
# Inputs:
# seismicDat: data.frame with columns:
#   east
#   north
#   seismicEst: estimate of sand volume fraction
#             NOTE: Doesn't need to be actual seismic data, just what is sampled 
#             preferentially with respect to.
# spatDat: data.frame with columns:
#   east
#   north
#   spatEst: estimate of the spatial component of the model. If 0, only seismic 
#            data is used to select the next drill location
# betaP: preferentiality parameter with respect to seismic estimates
# gammaP: preferentiality parameter with respect to spatial estimates
# seed: if not NULL, sets random number seed
# 
# Output: 
# a 3-vector with values: east, north, seismicEst
getPrefDrillLoc = function(seismicDat, 
                           spatDat=cbind(seismicDat$east, seismicDat$north, spatEst=0), 
                           betaP=0, gammaP=0, seed=NULL, 
                           repulsionType=c("none", "uniform"), bandwidth=10) {
  
  repulsionType = match.arg(repulsionType)
  
  # set random number seed
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # calculate grid width, delta (assume north and east grids are equal along with delta)
  eastGrid = sort(unique(seismicDat$east))
  delta = eastGrid[2] - eastGrid[1]
  
  # calculate linear predictor, sampling probabilities
  etaP = logit(seismicDat$seismicEst) * betaP + spatDat$spatEst * gammaP
  sampleProbs = exp(etaP)
  sampleProbs = sampleProbs/sum(sampleProbs)
  
  # sample centroid location
  sampleI = sample(1:nrow(seismicDat), 1, replace=FALSE, prob = sampleProbs)
  centEast = seismicDat$east[sampleI]
  centNorth = seismicDat$north[sampleI]
  
  # sample point location uniformly within grid cell
  east = centEast = runif(1) * delta - delta/2
  north = centNorth = runif(1) * delta - delta/2
  
  # use bilinear interpolation to get seismic estimate at that point
  seismicEst = bilinearInterp(rbind(c(east, north)), seismicDat)
  
  # return sampled drill location and associated seismic estimate
  c(east=east, north=north, seismicEst=seismicEst)
}

# Main simulation study functions -----

# NOTE:
# repI: ID of truth replicate. 100 total.
# sampleParI: ID of parameters related to well sampling. 4 total. Depends on:
#   repelAreaProp (1-2)
#   propVarCase (1-2)
# wellDatI: ID of full well dataset. 400 total. Depends on:
#   repelAreaProp (1-2)
#   propVarCase (1-2)
#   repI (1-100)
# fitModFunI: ID of model used to make predictions. 4+1 total including seismic only case.
# ModelFitI: ID of parameters for model fit. 4800+100 total. Depends on:
#   repelAreaProp (1-2) (doesn't affect seismic data only case)
#   propVarCase (1-2) (doesn't affect seismic data only case)
#   repI (1-100)
#   n (1-3) (doesn't affect seismic data only case)
#   fitModFunI (1-4 + seismic data only case)

# Inputs:
# Seed: random seed
# inputListFile: filename to save output as
# Outputs:
# simParList: List of simulation parameters with elements including:
#   parI: ID of simulation parameters (excluding n)
#   repI: ID of simulation replicate for this parI
#   n: number of well observations
#   fitModFunI: index of the model fitter
#   repelAreaProp: proportion of domain area repulsion around a single pt covers. Either 0 or 0.01
#   propVarCase: either sequential or clustered sampling
#   
setupSimStudy = function(seed=123, inputListFile="savedOutput/simStudy/simParList.RData") {
  set.seed(seed)
  
  fitModFunI = 1:length(getFitModFuns())
  repelAreaProp = c(0, 0.01)
  n = c(10, 20, 30)
  propVarCase = c("sequential", "clustered")
  repEffect = Inf
  nuggetVar = 0.1^2
  nsim = 100
  sigmaSq = 1
  prefPar = 3
  
  # Generate all combinations of varying parameters
  sampleParCombs = expand.grid(
    repelAreaProp = repelAreaProp,
    propVarCase = propVarCase,
    repEffect = repEffect, 
    nuggetVar = nuggetVar, 
    sigmaSq = sigmaSq, 
    prefPar = prefPar, 
    stringsAsFactors = FALSE
  )
  wellDatCombs = expand.grid(
    repelAreaProp = repelAreaProp,
    propVarCase = propVarCase,
    repI = 1:nsim, 
    repEffect = repEffect, 
    nuggetVar = nuggetVar, 
    sigmaSq = sigmaSq, 
    prefPar = prefPar, 
    stringsAsFactors = FALSE
  )
  modelFitCombs = expand.grid(
    repelAreaProp = repelAreaProp,
    propVarCase = propVarCase,
    fitModFunI = fitModFunI, 
    repI = 1:nsim, 
    n = n, 
    repEffect = repEffect, 
    nuggetVar = nuggetVar, 
    sigmaSq = sigmaSq, 
    prefPar = prefPar, 
    stringsAsFactors = FALSE
  )
  
  # add in IDs into the combination lists
  sampleParCombs$sampleParI = 1:nrow(sampleParCombs)
  wellDatCombs$wellDatI = 1:nrow(wellDatCombs)
  modelFitCombs$modelFitI = 1:nrow(modelFitCombs)
  
  # add in IDs from corresponding combination lists
  wellDatCombs = merge(wellDatCombs, sampleParCombs)
  modelFitCombs = merge(modelFitCombs, wellDatCombs)
  wellDatCombs$n = max(n)
  
  # Generate unique seeds for each replicated simulation and simulated dataset
  nDatasets = nrow(wellDatCombs)
  nFits = nrow(modelFitCombs)
  seedVec = sample.int(.Machine$integer.max, size = nDatasets + nFits, replace = FALSE)
  wellDatCombs$seed = seedVec[1:nDatasets]
  modelFitCombs$seed = seedVec[(nDatasets+1):(nDatasets + nFits)]
  
  # construct list of lists from each data.frame
  sampleParCombsList = dfToListOfLists(sampleParCombs)
  wellDatCombsList = dfToListOfLists(wellDatCombs)
  modelFitCombsList = dfToListOfLists(modelFitCombs)
  
  save(sampleParCombsList, sampleParCombs, wellDatCombsList, wellDatCombs, 
       modelFitCombsList, modelFitCombs, file=inputListFile)
}

# generates the well data for the simulation study
simStudySequentialSampler = function(i=1, regenData=FALSE) {
  
  # parameters
  out = load("savedOutput/simStudy/simParList.RData")
  thisPar = wellDatCombsList[[i]]
  prefPar = thisPar$prefPar
  wellDatI = thisPar$wellDatI
  repI = thisPar$repI
  sigmaSqErr = thisPar$nuggetVar
  repelAmount = thisPar$repEffect
  nWells = thisPar$n
  modelFitter = getFitModFuns()[[1]] # use SPDE model for predictions when sampling
  repelDist = repAreaToDist(thisPar$repelAreaProp)
  seed = thisPar$seed
  
  # if the well data already exists and we don't want to regenerate it, don't
  wellDatFile = paste0("savedOutput/simStudy/wellDat_par", wellDatI, "_rep", repI, ".RData")
  if(file.exists(wellDatFile) && !regenData) {
    return(invisible(NULL))
  }
  
  # get truth and data
  
  # seismic data
  out = readSurfaceRMS(paste0("../../synthetic_model/RegularizedPred_", repI, ".txt"))
  seismicDat = out$surfFrame
  
  # truth
  out = readSurfaceRMS(paste0("../../synthetic_model/RegularizedSand_", repI, ".txt"))
  truthDat = out$surfFrame
  
  # set repulsion parameters
  repelType = ifelse(repelDist == 0, "none", "rbf")
  if(repelDist != 0) {
    bwRepel = repelDist
  } else {
    bwRepel = NULL
  }
  
  
  # sample the wells
  browser()
  wellDat = wellSampler(truthDat, seismicDat, modelFitter, nWells=nWells, minN=4, 
                        predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
                        transform=logit, invTransform=expit, prefPar=prefPar, 
                        samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr, 
                        repelType=repelType, bwRepel=bwRepel, 
                        repelAmount=repelAmount, seed=seed, 
                        int.strategy="eb", strategy="gaussian")
  
  browser()
  save(wellDat, simPar=simParList[[maxJ]], 
       file=paste0("savedOutput/simStudy/wellDat_par", wellDatI, "_rep", repI, ".RData"))
}

runSimStudyI = function(i, significance=c(.8, .95), rerunModel=FALSE) {
  
  out = load("savedOutput/simStudy/simParList.RData")
  simPar = modelFitCombsList[[i]]
  
  sampleParI = simPar$sampleParI
  wellDatI = simPar$wellDatI
  modelFitI = simPar$modelFitI
  repI = simPar$repI
  fitModFunI = simPar$fitModFunI
  n = simPar$n
  repelAreaProp = simPar$repelAreaProp
  propVarCase = simPar$propVarCase
  prefPar = simPar$prefPar
  repEffect = simPar$repEffect
  nuggetVar = simPar$nuggetVar
  seed = simPar$seed
  
  parString = paste(
    "sampI", sampleParI,
    "datI", wellDatI,
    "mFitI", modelFitI, 
    "mod", fitModFunI,
    "n", n,
    "scen", substr(propVarCase, 1, 3),
    "pref", prefPar,
    "repP", repelAreaProp,
    "repE", repEffect,
    "nugV", nuggetVar,
    sep = "_"
  )
  
  browser() # check length of string
  fitModFun = getFitModFuns()[fitModFunI]
  
  # get seismic + well data and truth
  # list(xgrid=xgrid, ygrid=ygrid, surfMat=surfMat, 
  #      nx=nx, ny=ny, xstart=xstart, ystart=ystart, xend=xend, yend=yend)
  
  # seismic data
  out = readSurfaceRMS(paste0("../../synthetic_model/RegularizedPred_", repI, ".txt"))
  seismicDat = out$surfFrame
  
  # truth
  out = readSurfaceRMS(paste0("../../synthetic_model/RegularizedSand_", repI, ".txt"))
  truth = out$surfFrame
  
  # well data
  out = load(paste0("savedOutput/simStudy/wellDat_par", wellDatI, "_rep", repI, ".RData"))
  
  # interpolate truth to well points
  truthWells = bilinearInterp(seismicDat[,1:2], truth[,3], transform=logit, invTransform=expit)
  
  # Fit model if need be
  if(!file.exists() || rerunModel) {
    out = fitModFun(wellDat, seismicDat)
    predMat = out$predMat # doesn't include nugget
    predAggMat = out$predAggMat # doesn't include nugget?
    obsMat = out$obsMat # doesn't include nuggget
    
    save(predMat, predAggMat, truthWells, obsMat, 
         file=paste0("savedOutput/simStudy/modeRes_", resI, ".RData"))
  } else {
    out = load(paste0("savedOutput/simStudy/modeRes_", resI, ".RData"))
  }
  
  # calculate scoring rules and metrics based on predictions
  pwScoresMean = getScores(truth, estMat=predMat, significance=significance, doFuzzyReject=FALSE)
  pwScoresWorst = getScores(truth, estMat=predMat, significance=significance, doFuzzyReject=FALSE, aggFun=getWorst)
  aggScores = getScores(truth, estMat=predAggMat, significance=significance, doFuzzyReject=FALSE)
  
  # calculate informative summary statistics, first wrt wells, then over grid
  browser() # check wellDat variables
  corSeisTruthWells = cor(truthWells, wellDat$seismicEst)
  corSeisTruthTrue = cor(seismicDat[,3], truth[,3])
  varTruth = var(truth[,3])
  varSeis = var(seismicDat[,3])
  
  ests = rowMeans(predMat)
  estsWells = rowMeans(obsMat)
  varEst = var(ests)
  corEstTruthWells = cor(estsWells, truthWells)
  corEstTruthTrue = cor(ests, truth[,3])
  
  # save results
  
  
  save(pwScoresMean, pwScoresWorst, aggScores, 
       corSeisTruthWells, corSeisTruthTrue, varTruth, varSeis, 
       varEst, corEstTruthWells, corEstTruthTrue, file=paste0("scores", parString, ".RData"))
  
  invisible(NULL)
}











