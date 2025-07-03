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
  # combinations = expand.grid(
  #   fitModFuns = fitModFuns, 
  #   prefPar = prefPar,
  #   repelAreaProp = repelAreaProp,
  #   propVarCase = propVarCase,
  #   stringsAsFactors = FALSE
  # )
  combinations = expand.grid(
    fitModFunI = fitModFunI, 
    repelAreaProp = repelAreaProp,
    propVarCase = propVarCase,
    stringsAsFactors = FALSE
  )
  
  nComb = nrow(combinations)
  totalSim = nComb * nsim * length(n)
  
  # Generate unique seeds for each replicated simulation
  seedVec = sample.int(.Machine$integer.max, size = totalSim, replace = FALSE)
  
  # Repeat each combination nsim times
  simParList = vector("list", length = totalSim)
  idx = 1
  
  for(i in 1:nsim) {
    for(k in 1:length(n)) {
      for(j in 1:nComb) {
        simParList[[idx]] = list(
          parI = j, 
          repI = i, 
          fitModFunI = combinations$fitModFunI[j],
          n = n[k],
          repelAreaProp = combinations$repelAreaProp[j],
          propVarCase = combinations$propVarCase[j],
          prefPar = prefPar,
          repEffect = repEffect,
          nuggetVar = nuggetVar,
          seed = seedVec[idx]
        )
        idx = idx + 1
      }
    }
  }
  
  save(simParList, file=inputListFile)
}

# generates the well data for the simulation study
simStudySequentialSampler = function(i=1, seed=1, regenData=FALSE) {
  
  # parameters
  out = load("savedOutput/simStudy/simParList.RData")
  thisPar = simParList[[i]]
  prefPar = thisPar$prefPar
  parI = thisPar$parI
  repI = thisPar$repI
  sigmaSqErr = thisPar$nuggetVar
  repelAmount = thisPar$repEffect
  nWells = thisPar$n
  modelFitter = getFitModFuns()[[thisPar$fitModFunI]]
  repelDist = repAreaToDist(thisPar$repelAreaProp)
  
  # Generate unique seeds for each replicated simulation
  set.seed(seed)
  seedVec = sample.int(.Machine$integer.max, size = length(simParList), replace = FALSE)
  thisSeed = seedVec[i]
  set.seed(seed)
  
  # if the well data already exists and we don't want to regenerate it, don't
  wellDatFile = paste0("savedOutput/simStudy/wellDat_par", parI, "_rep", repI, ".RData")
  if(file.exists(wellDatFile) && !regenData) {
    return(invisible(NULL))
  }
  
  # get truth and data
  
  # seismic data
  out = readSurfaceRMS(paste0("/nr/sand/user/jpaige/synthetic_model/RegularizedPred_", repI, ".txt"))
  seismicDat = out$surfFrame
  
  # truth
  out = readSurfaceRMS(paste0("/nr/sand/user/jpaige/synthetic_model/RegularizedSand_", repI, ".txt"))
  truthDat = out$surfFrame
  
  # set repulsion parameters
  repelType = ifelse(repelDist == 0, "none", "rbf")
  if(repelDist != 0) {
    bwRepel = repelDist
  } else {
    bwRepel = NULL
  }
  
  # get max n fo this set of parameters
  maxJ = NULL
  for(j in 1:length(simParList)) {
    tempPar = simParList[[j]]
    tempPar$n = NULL
    thisTempPar = thisPar
    thisTempPar$n = NULL
    thisTempPar$seed = tempPar$seed
    
    if(identical(tempPar, thisTempPar)) {
      if(is.null(maxJ)) {
        maxJ = j
      } else if(simParList[[j]]$n > thisPar$n) {
        maxJ = j
      }
    }
  }
  maxN = simParList[[maxJ]]$n
  thisSeed = simParList[[maxJ]]$seed
  
  # sample the wells
  wellDat = wellSampler(truthDat, seismicDat, modelFitter, nWells=maxN, minN=4, 
                        predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
                        transform=logit, invTransform=expit, prefPar=prefPar, 
                        samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr, 
                        repelType=repelType, bwRepel=bwRepel, 
                        repelAmount=repelAmount, seed=thisSeed)
  
  save(wellDat, simPar=simParList[[maxJ]], file=paste0("savedOutput/simStudy/wellDat_par", parI, "_rep", repI, ".RData"))
}

runSimStudyI = function(i, significance=c(.8, .95), rerunModel=FALSE) {
  
  out = load("savedOutput/simStudy/simParList.RData")
  simPar = simParList[[i]]
  
  parI = simPar$parI
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
    "i", parI,
    "j", repI,
    "mod", fitModFunI,
    "n", n,
    "samp", substr(propVarCase, 1, 3),
    "pref", prefPar,
    "repP", repelAreaProp,
    "repE", repEffect,
    "nugV", nuggetVar,
    sep = "_"
  )
  
  fitModFun = getFitModFuns()[fitModFunI]
  
  # get seismic + well data and truth
  # list(xgrid=xgrid, ygrid=ygrid, surfMat=surfMat, 
  #      nx=nx, ny=ny, xstart=xstart, ystart=ystart, xend=xend, yend=yend)
  
  # seismic data
  out = readSurfaceRMS(paste0("/nr/sand/user/jpaige/synthetic_model/RegularizedPred_", parI, ".txt"))
  seismicDat = out$surfFrame
  
  # truth
  out = readSurfaceRMS(paste0("/nr/sand/user/jpaige/synthetic_model/RegularizedSand_", repI, ".txt"))
  truth = out$surfFrame
  
  # well data
  out = load(paste0("savedOutput/simStudy/wellDat_par", parI, "_rep", repI, ".RData"))
  
  # interpolate truth to well points
  truthWells = bilinearInterp(seismicDat[,1:2], truth[,3], transform=logit, invTransform=expit)
  
  # Fit model if need be
  if(!file.exists() || rerunModel) {
    out = fitModFun(wellDat, seismicDat)
    predMat = out$predMat # doesn't include nugget
    predAggMat = out$predAggMat # doesn't include nugget?
    obsMat = out$obsMat # doesn't include nuggget
    
    save(predMat, predAggMat, truthWells, obsMat, file=paste0("savedOutput/simStudy/modeRes_", resI, ".RData"))
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











