# functions for the simulation study

# construct a rectangular prediction grid over a domain
# # xLims: a 2-vector with the min and max x value
# # yLims: a 2-vector with the min and max y value
# upScaleFac: upscale factor (new number of pred pts = normal / 2^upScaleFac)
getSimStudyPredGrid = function(upScaleFac=2, getGoodCoords=FALSE) {
  
  predGridFile = paste0("savedOutput/global/simPredGrid_", upScaleFac, ".RData")
  
  if(!file.exists(predGridFile)) {
    out = readSurfaceRMS("data/seisTruthReplicates/RegularizedPred_1.txt")
    regPred = out$surfFrame
    # xlims = c(-12.5, 15012.5)
    # ylims = c(-8.3335, 5008.4336)
    
    allXs = sort(unique(regPred[,1]))
    allYs = sort(unique(regPred[,2]))
    
    newXs = allXs[seq(from=1, to=length(allXs), by=upScaleFac)]
    newYs = allYs[seq(from=1, to=length(allYs), by=upScaleFac)]
    
    goodXs = regPred[,1] %in% newXs
    goodYs = regPred[,2] %in% newYs
    goodCoords = goodXs & goodYs
    
    predGrid = as.matrix(regPred[goodCoords, 1:2])
    
    save(predGrid, goodCoords, file=predGridFile)
  } else {
    out = load(predGridFile)
  }
  
  if(getGoodCoords) {
    list(predGrid=predGrid, goodCoords=goodCoords)
  } else {
    predGrid
  }
}

subsampleSimStudyGrid = function(gridDat, upScaleFac=3) {
  
  allXs = sort(unique(gridDat[,1]))
  allYs = sort(unique(gridDat[,2]))
  
  newXs = allXs[seq(from=1, to=length(allXs), by=upScaleFac)]
  newYs = allYs[seq(from=1, to=length(allYs), by=upScaleFac)]
  
  goodXs = gridDat[,1] %in% newXs
  goodYs = gridDat[,2] %in% newYs
  goodCoords = goodXs & goodYs
  
  goodCoords
  
}

getFitModFuns = function() {
  funs = list()
  
  # basic SPDE (no seismic)
  # funs = c(funs, list(function(...) {fitSPDEsimDat(addKDE=TRUE, ...)}))
  
  # SPDE with seismic
  funs = c(funs, list(fitSPDEsimDat))
  
  # SPDE with seismic + kde covariate
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
  xlims = c(-12.5, 15012.5)
  ylims = c(-8.3335, 5008.4336)
  domainArea = diff(xlims) * diff(ylims)
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
# 
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
setupSimStudy = function(adaptScen=c("batch", "adaptPref", "adaptVar")) {
  adaptScen = match.arg(adaptScen)
  
  if(adaptScen == "batch") {
    seed = 1
  } else if(adaptScen == "adaptPref") {
    seed = 2
  } else {
    seed = 3
  }
  set.seed(seed)
  
  adaptScenCap = str_to_title(adaptScen)
  inputListFile = paste0("savedOutput/simStudy/simParList", adaptScenCap, ".RData")
  
  if(adaptScen == "batch") {
    n = c(20, 40, 60)
    propVarCase = c("realistic", "uniform")
    prefPar = c(1.5, 3)
    repelAreaProp = c(0, 0.001, 0.01)
  } else {
    n = c(10, 20, 30)
    propVarCase = c("spde", "self")
    prefPar = c(3)
    repelAreaProp = c(0, 0.01)
  }
  fitModFunI = 1:length(getFitModFuns())
  
  
  repEffect = Inf
  nuggetVar = 0.1^2
  nsim = 100
  sigmaSq = 1
  
  # Generate all combinations of varying parameters
  
  # parameters used to sample the wells (ID corresponds to well sampling scenario)
  sampleParCombs = expand.grid(
    repelAreaProp = repelAreaProp,
    propVarCase = propVarCase,
    sampModFunI = fitModFunI, 
    repEffect = repEffect, 
    nuggetVar = nuggetVar, 
    sigmaSq = sigmaSq, 
    prefPar = prefPar, 
    stringsAsFactors = FALSE
  )
  
  # well dataset IDs, including 100 replicate IDs + corresponding parameters (+ eventually max n, seed)
  wellDatCombs = expand.grid(
    repelAreaProp = repelAreaProp,
    propVarCase = propVarCase,
    sampModFunI = fitModFunI, 
    repI = 1:nsim, 
    repEffect = repEffect, 
    nuggetVar = nuggetVar, 
    sigmaSq = sigmaSq, 
    prefPar = prefPar, 
    stringsAsFactors = FALSE
  )
  
  # model fit IDs, including full well dataset info + n (well dataset subset size) + model
  modelFitCombs = expand.grid(
    repelAreaProp = repelAreaProp,
    propVarCase = propVarCase,
    sampModFunI = fitModFunI, 
    fitModFunI = fitModFunI, 
    repI = 1:nsim, 
    n = n, 
    repEffect = repEffect, 
    nuggetVar = nuggetVar, 
    sigmaSq = sigmaSq, 
    prefPar = prefPar, 
    stringsAsFactors = FALSE
  )
  
  # remove cases we don't consider:
  # - only use 1 beta value for uniform case
  # - don't use Diggle model for variance adaptive sampling?
  # - don't use repelAreaProp == .01 for n in 100,500 or repelAreaProp==.001 for n == 500
  # - don't generate data based on propVarCase in adaptive settings
  # - don't use sampModFunI in batch settings
  # - either sampModFunI == fitModFunI or sampModFunI == 1 in adaptive case (for modelTab)
  
  removeBadRows = function(tab, modelTab=FALSE) {
    # first remove extra beta value for uniform case
    badBetas = (tab$propVarCase == "uniform") & (tab$prefPar == 3)
    
    # then remove bad values of repelAreaProp
    if(!is.null(tab$n)) {
      badRepArea = tab$repelAreaProp * tab$n > 0.5
    } else {
      badRepArea = rep(FALSE, nrow(tab))
    }
    
    # remove rows (and column) based on sampModFunI in batch settings
    if(adaptScen == "batch") {
      badSampMod = tab$sampModFunI != fitModFunI[1]
      tab$sampModFunI = NULL
    } else {
      badSampMod = rep(FALSE, nrow(tab))
    }
    
    if(!modelTab) {
      # just for the data table (everything but modelFitCombs)
      
      # remove rows (and column) based on propVarCase in adaptive settings
      if(adaptScen == "batch") {
        badPropVar = rep(FALSE, nrow(tab))
      } else {
        badPropVar = tab$propVarCase != propVarCase[1]
        tab$propVarCase = NULL
      }
      
      badRows = badPropVar
    } else {
      # just for modelFitCombs
      
      # either sampModFunI == fitModFunI or sampModFunI == 1 in adaptive case 
      # depending on propVarCase (for modelTab)
      
      goodSampMod = rep(TRUE, nrow(tab))
      
      if(adaptScen != "batch") {
        goodSampMod[tab$propVarCase == "spde"] = 
          tab$sampModFunI[tab$propVarCase == "spde"] == 1
        goodSampMod[tab$propVarCase == "self"] = 
          tab$sampModFunI[tab$propVarCase == "self"] == tab$fitModFunI[tab$propVarCase == "self"]
      }
      
      badRows = !goodSampMod
    }
    
    badRows = badBetas | badRepArea | badSampMod | badRows
    tab[!badRows,]
  }
  
  sampleParCombs = removeBadRows(sampleParCombs)
  wellDatCombs = removeBadRows(wellDatCombs)
  modelFitCombs = removeBadRows(modelFitCombs, modelTab=TRUE)
  
  # add in IDs into the combination lists
  sampleParCombs$sampleParI = 1:nrow(sampleParCombs)
  wellDatCombs$wellDatI = 1:nrow(wellDatCombs)
  modelFitCombs$modelFitI = 1:nrow(modelFitCombs)
  
  # add in IDs from corresponding combination lists
  wellDatCombs = merge(wellDatCombs, sampleParCombs)
  wellDatCombs = wellDatCombs[order(wellDatCombs$wellDatI),]
  modelFitCombs = merge(modelFitCombs, wellDatCombs)
  modelFitCombs = modelFitCombs[order(modelFitCombs$modelFitI),]
  wellDatCombs$n = max(n)
  wellDatCombs$n[wellDatCombs$repelAreaProp * wellDatCombs$n > 0.5] = n[length(n)-1]
  
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

simStudyWellSamplerPar = function(i=1, adaptScen=c("batch", "adaptPref", "adaptVar"), 
                                  regenData=FALSE, verbose=FALSE) {
  
  tryCatch(simStudyWellSampler(i, adaptScen, regenData, verbose), 
           error = function(e) {
             logfile <- paste0("savedOutput/simStudy/well_", adaptScen, "_", i, "_err.txt")
             sink(logfile)
             cat("Error at i =", i, ":\n")
             cat(paste("Call stack:\n", paste(deparse(sys.calls()), collapse = "\n")), "\n")
             cat(conditionMessage(e), "\n")
             sink()
             stop("error")
           })
  
}

# generates sequentially-sampled well data for the simulation study
simStudyWellSampler = function(i=1, adaptScen=c("batch", "adaptPref", "adaptVar"), 
                               regenData=FALSE, verbose=FALSE) {
  adaptScen = match.arg(adaptScen)
  
  if(verbose) {
    print(paste0("generating well data for i: ", i, ", adapt scenario: ", adaptScen))
  }
  
  adaptScenCap = str_to_title(adaptScen)
  inputListFile = paste0("savedOutput/simStudy/simParList", adaptScenCap, ".RData")
  
  # parameters
  out = load(inputListFile)
  thisPar = wellDatCombsList[[i]]
  propVarCase = thisPar$propVarCase
  prefPar = thisPar$prefPar
  sampleParI = thisPar$sampleParI
  wellDatI = thisPar$wellDatI
  repI = thisPar$repI
  sigmaSqErr = thisPar$nuggetVar
  sigmaSq = thisPar$sigmaSq
  repelAmount = thisPar$repEffect
  nWells = thisPar$n
  if(adaptScen != "batch") {
    modelFitter = getFitModFuns()[[thisPar$sampModFunI]]
  }
  repelDist = repAreaToDist(thisPar$repelAreaProp)
  seed = thisPar$seed
  
  # if the well data already exists and we don't want to regenerate it, don't
  wellDatFile = paste0("savedOutput/simStudy/wellDat/wellDat_", adaptScen, "_par", sampleParI, "_rep", repI, ".RData")
  if(file.exists(wellDatFile) && !regenData) {
    return(invisible(NULL))
  }
  
  # get truth and data
  
  # seismic data
  out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedPred_", repI, ".txt"), force01=TRUE)
  seismicDat = out$surfFrame
  
  # truth
  out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", repI, ".txt"), force01=TRUE)
  truthDat = out$surfFrame
  
  # subsample
  goodCoords = subsampleSimStudyGrid(seismicDat)
  seismicDat = seismicDat[goodCoords,]
  truth = truth[goodCoords,]
  
  # set repulsion parameters
  repelType = ifelse(repelDist == 0, "none", "rbf")
  if(repelDist != 0) {
    bwRepel = repelDist
  } else {
    bwRepel = NULL
  }
  
  if(adaptScen != "batch") {
    # sample the wells
    # wellDat = wellSampler(truthDat, seismicDat, modelFitter, nWells=nWells, minN=4,
    #                       predGrid=cbind(east=seismicDat$east, north=seismicDat$north),
    #                       transform=logit, invTransform=expit, prefPar=prefPar,
    #                       samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr,
    #                       repelType=repelType, bwRepel=bwRepel,
    #                       repelAmount=repelAmount, seed=seed,
    #                       int.strategy="eb", strategy="gaussian")
    
    # for testing purposes:
    wellDat = wellSampler(truthDat=truthDat, seismicDat=seismicDat, 
                          modelFitter=modelFitter, nWells=5, minN=4,
                          predGrid=cbind(east=seismicDat$east, north=seismicDat$north),
                          transform=logit, invTransform=expit, prefPar=prefPar,
                          samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr,
                          repelType=repelType, bwRepel=bwRepel,
                          repelAmount=repelAmount, seed=seed,
                          int.strategy="eb", strategy="gaussian")
    
    # profvis(wellDat <- wellSampler(truthDat, seismicDat, modelFitter, nWells=5, minN=4,
    #                                predGrid=cbind(east=seismicDat$east, north=seismicDat$north),
    #                                transform=logit, invTransform=expit, prefPar=prefPar,
    #                                samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr,
    #                                repelType=repelType, bwRepel=bwRepel,
    #                                repelAmount=repelAmount, seed=seed,
    #                                int.strategy="eb", strategy="gaussian"))
  } else {
    
    if(propVarCase == "realistic") {
      # indep
      otherRepI = ((repI + 1) %% 100) + 1
      out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", otherRepI, ".txt"), force01=TRUE)
      indepDat = out$surfFrame
      
      # subsample
      indepDat = indepDat[goodCoords,]
      
      # standardize seismic, truth, and indep data on logit scale
      seismicDatStd = seismicDat
      truthDatStd = truthDat
      indepDatStd = indepDat
      seismicDatStd[,3] = logit(seismicDatStd[,3])
      truthDatStd[,3] = logit(truthDatStd[,3])
      indepDatStd[,3] = logit(indepDatStd[,3])
      seismicDatStd[,3] = (seismicDat[,3] - mean(seismicDat[,3]))/sd(seismicDat[,3])
      truthDatStd[,3] = (truthDat[,3] - mean(truthDat[,3]))/sd(truthDat[,3])
      indepDatStd[,3] = (indepDat[,3] - mean(indepDat[,3]))/sd(indepDat[,3])
      
      # combine them into a realistic mix and convert back to [0,1] scale
      sampleDat = seismicDat
      sampleDat[,3] = 0.5 * seismicDatStd[,3] + 0.25 * truthDatStd[,3] + 0.25 * indepDatStd[,3]
    } else if(propVarCase == "uniform") {
      # in this case, seismic data doesn't matter that much
      sampleDat = seismicDat
      sampleDat[,3] = rep(0, nrow(sampleDat))
    } else {
      stop("unrecognized propVarCase")
    }
    
    sampleDat[,3] = expit(sqrt(sigmaSq) * sampleDat[,3])
    
    # do batch sampling
    wellDat = wellSampler(nWells=nWells, wellDat=NULL, seismicDat=seismicDat, 
                     predGrid=as.matrix(sampleDat[,1:2]), truthDat=truthDat, 
                     preds=sampleDat[,3], 
                     transform=logit, invTransform=expit, prefPar=prefPar, 
                     samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr, 
                     repelType=repelType, bwRepel=bwRepel, 
                     rbf="uniform", repelAmount=Inf, 
                     seed=seed, int.strategy="eb", strategy="gaussian")$wellDat
    
    
  }
  
  
  save(wellDat, simPar=thisPar,
       file=wellDatFile)
  invisible(NULL)
}

runSimStudyIPar = function(i, significance=c(.8, .95), 
                           adaptScen=c("batch", "adaptPref", "adaptVar"), 
                           regenData=FALSE, verbose=FALSE) {
  
  tmplogfile <- paste0("savedOutput/simStudy/scores_", adaptScen, "_", i, "_tmp.txt")
  sink(tmplogfile)
  cat("Started scores job i =", i, ":\n")
  sink()
  
  tryCatch(runSimStudyI(i, significance, adaptScen, regenData, verbose), 
           error = function(e) {
             logfile <- paste0("savedOutput/simStudy/scores_", adaptScen, "_", i, "_err.txt")
             sink(logfile)
             cat("Error at i =", i, ":\n")
             cat(paste("Call stack:\n", paste(deparse(sys.calls()), collapse = "\n")), "\n")
             cat(conditionMessage(e), "\n")
             sink()
             stop("error")
           })
  
  # remove the tmp log file to indicate the job is complete
  system(paste0("rm ", tmplogfile))
  
}

runSimStudyI = function(i, significance=c(.8, .95), 
                        adaptScen=c("batch", "adaptPref", "adaptVar"), 
                        regenData=FALSE, verbose=FALSE) {
  startT = proc.time()[3]
  adaptScen = match.arg(adaptScen)
  
  if(verbose) {
    print(paste0("generating model predictions for i: ", i, ", adapt scenario: ", adaptScen))
  }
  
  adaptScenCap = str_to_title(adaptScen)
  inputListFile = paste0("savedOutput/simStudy/simParList", adaptScenCap, ".RData")
  out = load(inputListFile)
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
  fitModFun = getFitModFuns()[[fitModFunI]]
  
  # get seismic + well data and truth
  # list(xgrid=xgrid, ygrid=ygrid, surfMat=surfMat, 
  #      nx=nx, ny=ny, xstart=xstart, ystart=ystart, xend=xend, yend=yend)
  
  # seismic data
  out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedPred_", repI, ".txt"), force01=TRUE)
  seismicDat = out$surfFrame
  
  # truth
  out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", repI, ".txt"), force01=TRUE)
  truth = out$surfFrame
  
  # subsample
  goodCoords = subsampleSimStudyGrid(seismicDat)
  seismicDat = seismicDat[goodCoords,]
  truth = truth[goodCoords,]
  
  # well data
  wellDatFile = paste0("savedOutput/simStudy/wellDat/wellDat_", adaptScen, "_par", sampleParI, "_rep", repI, ".RData")
  out = load(wellDatFile)
  
  # subset well data based on n for this run
  wellDat = wellDat[1:n,]
  
  # interpolate truth to well points
  truthWells = bilinearInterp(wellDat[,1:2], truth, transform=logit, invTransform=expit)
  
  # Fit model and calculate scores if need be
  scoresFile = paste0("savedOutput/simStudy/scores/scores_", adaptScen, "_", i, ".RData")
  if(!file.exists(scoresFile) || regenData) {
    if(i != 4) {
      inputList = list(wellDat, seismicDat)
    } else {
      repDist = repAreaToDist(repelAreaProp)
      inputList = list(wellDat, seismicDat, repelDist=repDist)
    }
    
    out = do.call("fitModFun", list(wellDat, seismicDat))
    
    if(fitModFunI == 4) {
      predMat = out$predMat.y # doesn't include nugget
      predAggMat = out$pred.yAggMat # doesn't include nugget?
      obsMat = out$obsMat.y # doesn't include nugget
    } else {
      predMat = out$predMat # doesn't include nugget
      predAggMat = out$predAggMat # doesn't include nugget?
      obsMat = out$obsMat # doesn't include nugget
    }
    
    
    # calculate scoring rules and metrics based on predictions
    pwScoresMean = getScores(truth[,3], estMat=predMat, significance=significance)
    pwScoresWorst = getScores(truth[,3], estMat=predMat, significance=significance, aggFun=getWorst)
    aggScores = getScores(mean(truth[,3]), estMat=matrix(predAggMat, nrow=1), significance=significance)
    pwScoresMax = getScores(max(truth[,3]), estMat=matrix(apply(predMat, 2, max), nrow=1), significance=significance)
    pwScoresMin = getScores(min(truth[,3]), estMat=matrix(apply(predMat, 2, min), nrow=1), significance=significance)
    
    # calculate informative summary statistics, first wrt wells, then over grid
    corSeisTruthWells = cor(truthWells, wellDat[,3])
    corSeisTruthTrue = cor(seismicDat[,3], truth[,3])
    varTruth = var(truth[,3])
    varSeis = var(seismicDat[,3])
    
    ests = rowMeans(predMat)
    estsWells = rowMeans(obsMat)
    varEst = var(ests)
    corEstTruthWells = cor(estsWells, truthWells)
    corEstTruthTrue = cor(ests, truth[,3])
    
    endT = proc.time()[3]
    totT = endT - startT
    
    # save results
    save(pwScoresMean, pwScoresWorst, aggScores, pwScoresMax, pwScoresMin, 
         corSeisTruthWells, corSeisTruthTrue, varTruth, varSeis, 
         varEst, corEstTruthWells, corEstTruthTrue, totT=totT, file=scoresFile)
  }
  
  invisible(NULL)
}

# get scores from seismic data
getSeismicEsts = function(i, regenData=FALSE, significance=c(.8, .95)) {
  
  startT = proc.time()[3]
  
  # Fit model and calculate scores if need be
  scoresFile = paste0("savedOutput/simStudy/scores/scores_seismic_rep", i, ".RData")
  if(!file.exists(scoresFile) || regenData) {
    
    # get seismic data
    out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedPred_", i, ".txt"), force01=TRUE)
    seismicDat = out$surfFrame
    
    # truth
    out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", i, ".txt"), force01=TRUE)
    truth = out$surfFrame
    
    # subsample
    goodCoords = subsampleSimStudyGrid(seismicDat)
    seismicDat = seismicDat[goodCoords,]
    truth = truth[goodCoords,]
    
    # generate prediction "distribution"
    predMat = cbind(seismicDat[,3], seismicDat[,3])
    meanSeis = mean(seismicDat[,3])
    predAggMat = matrix(c(meanSeis, meanSeis), nrow=1)
    
    # calculate scoring rules and metrics based on predictions
    system.time(pwScoresMean <- getScores(truth[,3], estMat=predMat, significance=significance))[3] # 33 seconds?!
    pwScoresMean = getScores(truth[,3], estMat=predMat, significance=significance)
    pwScoresWorst = getScores(truth[,3], estMat=predMat, significance=significance, aggFun=getWorst)
    aggScores = getScores(mean(truth[,3]), estMat=matrix(predAggMat, nrow=1), significance=significance)
    pwScoresMax = getScores(max(truth[,3]), estMat=matrix(apply(predMat, 2, max), nrow=1), significance=significance)
    pwScoresMin = getScores(min(truth[,3]), estMat=matrix(apply(predMat, 2, min), nrow=1), significance=significance)
    
    # calculate informative summary statistics, first wrt wells, then over grid
    corSeisTruthWells = NA
    corSeisTruthTrue = cor(seismicDat[,3], truth[,3])
    varTruth = var(truth[,3])
    varSeis = var(seismicDat[,3])
    
    ests = seismicDat[,3]
    estsWells = NA
    varEst = var(ests)
    corEstTruthWells = NA
    corEstTruthTrue = cor(ests, truth[,3])
    
    endT = proc.time()[3]
    totT = endT - startT
    
    # save results
    save(pwScoresMean, pwScoresWorst, aggScores, pwScoresMax, pwScoresMin, 
         corSeisTruthWells, corSeisTruthTrue, varTruth, varSeis, 
         varEst, corEstTruthWells, corEstTruthTrue, totT=totT, file=scoresFile)
  }
}

# Final sim study funs ----

# generate simulation parameter tables (simParList*.RData files)
setupParSimStudy = function() {
  setupSimStudy("batch")
  setupSimStudy("adaptPref")
  setupSimStudy("adaptVar")
}

# generates well data for the simulation study
getWellDatSimStudy = function(nCores=8, adaptScen=c("batch", "adaptPref", "adaptVar"), 
                              doPar=TRUE, regenData=FALSE) {
  adaptScen = match.arg(adaptScen)
  
  # load simulation parameters
  adaptScenCap = str_to_title(adaptScen)
  inputListFile = paste0("savedOutput/simStudy/simParList", adaptScenCap, ".RData")
  out = load(inputListFile)
  
  # determine how many scenarios to simulate for:
  n = nrow(wellDatCombs)
  
  if(doPar) {
    # start parallel cluster
    cl = makeCluster(nCores)
    clusterEvalQ(cl, source("R/setup.R"))
    
    # generate well data in parallel
    tmp = parLapply(cl, 1:n, simStudyWellSamplerPar, adaptScen=adaptScen, regenData=regenData, verbose=TRUE)
    
    # remember to stop the cluster
    stopCluster(cl)
  } else {
    for(i in 1:n) {
      simStudyWellSampler(i, adaptScen=adaptScen, regenData=regenData, verbose=TRUE)
    }
  }
  
  invisible(NULL)
}

# fits models based on generated well data for the simulation study
fitModsSimStudy = function(nCores=8, adaptScen=c("batch", "adaptPref", "adaptVar"), maxRepI=100, 
                           doPar=TRUE, regenData=FALSE) {
  adaptScen = match.arg(adaptScen)
  
  # load simulation parameters
  adaptScenCap = str_to_title(adaptScen)
  inputListFile = paste0("savedOutput/simStudy/simParList", adaptScenCap, ".RData")
  out = load(inputListFile)
  
  # figure out which parameter sets have repI <= maxRepI
  is = 1:nrow(modelFitCombs)
  is = is[modelFitCombs$repI <= maxRepI]
  
  if(doPar) {
    # start parallel cluster
    cl = makeCluster(nCores)
    clusterEvalQ(cl, source("R/setup.R"))
    
    # generate well data in parallel
    tmp = parLapply(cl, is, runSimStudyIPar, adaptScen=adaptScen, regenData=regenData, verbose=TRUE)
    
    # remember to stop the cluster
    stopCluster(cl)
  } else {
    for(i in is) {
      runSimStudyI(i, adaptScen=adaptScen, regenData=regenData, verbose=TRUE)
    }
  }
}

getAllSeisEstsSimStudy = function(regenData=FALSE, significance=c(.8, .95)) {
  
  for(i in 1:100) {
    print(paste0("generating seismic estimates for iteration ", i, "/100"))
    getSeismicEsts(i, regenData=regenData, significance=significance)
  }
  
  invisible(NULL)
}

showSimStudyRes = function(adaptScen=c("batch", "adaptPref", "adaptVar"), maxRepI=100) {
  
  adaptScen = match.arg(adaptScen)
  
  # load simulation parameters and module run parameters
  adaptScenCap = str_to_title(adaptScen)
  inputListFile = paste0("savedOutput/simStudy/simParList", adaptScenCap, ".RData")
  out = load(inputListFile)
  
  # figure out which parameter sets have repI <= maxRepI
  is = 1:nrow(modelFitCombs)
  is = is[modelFitCombs$repI <= maxRepI]
  subModelCombs = modelFitCombs[modelFitCombs$repI <= maxRepI,]
  
  # 
  # for(i in is) {
  #   # check the parameters for this run
  #   runInfo = modelFitCombs[i,]
  #   parI = runInfo$sampleParI
  #   
  #   
  #   scoresFile = paste0("savedOutput/simStudy/scores/scores_", adaptScen, "_", i, ".RData")
  #   # save(pwScoresMean, pwScoresWorst, aggScores, pwScoresMax, pwScoresMin, 
  #   #      corSeisTruthWells, corSeisTruthTrue, varTruth, varSeis, 
  #   #      varEst, corEstTruthWells, corEstTruthTrue, totT=totT, file=scoresFile)
  #   
  #   out = load(scoresFile)
  # }
  
  # modI: the type of model being fit (0 = seismic only)
  # thisParI: the parameter scenario we will aggregate results for
  getModScores = function(modI, thisParI) {
    
    # indices of the scores files we'll need to load and the corresponding 
    # rows of modelFitCombs
    if(modI != 0) {
      thisIs = which((subModelCombs$sampleParI == thisParI) & (subModelCombs$fitModFunI == modI))
      
      # make cases go from small n to large n
      thisNs = subModelCombs$n[thisIs]
      ordI = order(thisNs)
      thisIs = thisIs[ordI]
    } else {
      thisIs = 1:maxRepI
    }
    
    pwScoresMeanAll = c()
    pwScoresWorstAll = c()
    aggScoresAll = c()
    pwScoresMaxAll = c()
    pwScoresMinAll = c()
    corSeisTruthWellsAll = c()
    corSeisTruthTrueAll = c()
    varTruthAll = c()
    varSeisAll = c()
    varEstAll = c()
    corEstTruthWellsAll = c()
    corEstTruthTrueAll = c()
    totTAll = c()
    repIAll = c()
    nAll = c()
    for(i in thisIs) {
      thesePar = subModelCombs[i,]
      thisModelI = thesePar$modelFitI
      
      if(modI == 0) {
        # get the seismic only estimate. This doesn't depend on well data at all
        # get the maximum number of reps
        scoresFile = paste0("savedOutput/simStudy/scores/scores_seismic_rep", i, ".RData")
        
      } else {
        scoresFile = paste0("savedOutput/simStudy/scores/scores_", adaptScen, "_", thisModelI, ".RData")
      }
      # save(pwScoresMean, pwScoresWorst, aggScores, pwScoresMax, pwScoresMin, 
      #      corSeisTruthWells, corSeisTruthTrue, varTruth, varSeis, 
      #      varEst, corEstTruthWells, corEstTruthTrue, totT=totT, file=scoresFile)
      
      out = load(scoresFile)
      
      pwScoresMeanAll = rbind(pwScoresMeanAll, pwScoresMean)
      pwScoresWorstAll = rbind(pwScoresWorstAll, pwScoresWorst)
      aggScoresAll = rbind(aggScoresAll, aggScores)
      pwScoresMaxAll = rbind(pwScoresMaxAll, pwScoresMax)
      pwScoresMinAll = rbind(pwScoresMinAll, pwScoresMin)
      corSeisTruthWellsAll = c(corSeisTruthWellsAll, corSeisTruthWells)
      corSeisTruthTrueAll = c(corSeisTruthTrueAll, corSeisTruthTrue)
      varTruthAll = c(varTruthAll, varTruth)
      varSeisAll = c(varSeisAll, varSeis)
      varEstAll = c(varEstAll, varEst)
      corEstTruthWellsAll = c(corEstTruthWellsAll, corEstTruthWells)
      corEstTruthTrueAll = c(corEstTruthTrueAll, corEstTruthTrue)
      totTAll = c(totTAll, totT)
      repIAll = c(repIAll, thesePar$repI)
      nAll = c(nAll, thesePar$n)
    }
    
    list(pwScoresMeanAll=pwScoresMeanAll, pwScoresWorstAll=pwScoresWorstAll, 
         aggScoresAll=aggScoresAll, pwScoresMaxAll=pwScoresMaxAll, 
         pwScoresMinAll=pwScoresMinAll, corSeisTruthWellsAll=corSeisTruthWellsAll, 
         corSeisTruthTrueAll=corSeisTruthTrueAll, varTruthAll=varTruthAll, 
         varSeisAll=varSeisAll, varEstAll=varEstAll, 
         corEstTruthWellsAll=corEstTruthWellsAll, corEstTruthTrueAll=corEstTruthTrueAll, 
         totTAll=totTAll, repIAll=repIAll, nAll=nAll)
  }
  
  # for each set of sampling parameters, show results
  for(i in 1:nrow(sampleParCombs)) {
    
    thesePar = sampleParCombs[i,]
    
    dirName = paste0("figures/simStudy/")
    
    sampleParI = thesePar$sampleParI # same as i
    fitModFunI = thesePar$fitModFunI
    repelAreaProp = thesePar$repelAreaProp
    propVarCase = thesePar$propVarCase
    prefPar = thesePar$prefPar
    repEffect = thesePar$repEffect
    fileRoot = paste0("_i", sampleParI, 
                      "_case", propVarCase, 
                      "_modI", fitModFunI, 
                      "_repA", repelAreaProp, 
                      "_pref", prefPar)
    
    # get scores for each model
    seisScores = getModScores(0, i)
    spdeScores = getModScores(1, i)
    spdeKernScores = getModScores(2, i)
    diggleScores = getModScores(3, i)
    watsonScores = getModScores(4, i)
    # score names:
    # pwScoresMeanAll, pwScoresWorstAll, aggScoresAll, pwScoresMaxAll, 
    # pwScoresMinAll, corSeisTruthWellsAll, corSeisTruthTrueAll, varTruthAll, 
    # varSeisAll, varEstAll, corEstTruthWellsAll, corEstTruthTrueAll, 
    # totTAll, repIAll, nAll
    
    modCols = c("black", "lightblue", "blue", "purple", "orange", "green")
    pch = c(5, 15:19)
    
    browser()
    seisScores$aggScoresAll
    
    
    
  }
  
}





