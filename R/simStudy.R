# functions for the simulation study

# THIS FUNCTION NO LONGER USED!!!
# construct a rectangular prediction grid over a domain
# # xLims: a 2-vector with the min and max x value
# # yLims: a 2-vector with the min and max y value
# upScaleFac: upscale factor (new number of pred pts = normal / 2^upScaleFac)
getSimStudyPredGrid = function(upScaleFac=2, getGoodCoords=FALSE) {
  stop("deprecated")
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

# low: lower end of CI on variance
# high: upper end of CI on variance
getPrecPrior = function(low=.05^2, high=0.5^2, alpha=.05) {
  
  optFun = function(par) {
    a = exp(par[1])
    b = exp(par[2])
    mean((log(c(1/qgamma(1-alpha/2, a, b), 1/qgamma(alpha/2, a, b))) - 
            log(c(low, high)))^2)
  }
  
  opt = optim(c(0, 0), fn=optFun)
  
  optPar = exp(opt$par)
  a = optPar[1]
  b = optPar[2]
  print(paste0("(a=", round(a, 4), ", b=", round(b,4), ") yields ", 1-alpha, "x100% CI for variance of (", 
               round(1/qgamma(1-alpha/2, a, b), 4), ", ", round(1/qgamma(0.5, a, b), 4), ", ", round(1/qgamma(alpha/2, a, b), 4), ")"))
  
  
  c(a=a, b=b, low=1/qgamma(1-alpha/2, a, b), high=1/qgamma(alpha/2, a, b))
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
  
  # SPDE with design sampling probabilities accounted for
  funs = c(funs, list(function(...) {fitSPDEsimDat(addLogProbs=TRUE, ...)}))
  
  funs
}

getFitModName = function(ind) {
  c("Seismic", "SPDE", "SPDEK", "Diggle", "Watson", "SPDED")[ind+1]
}


# proportion of domain area -> repel radius bandwidth
repAreaToDist = function(repArea=.01, anisFac=1) {
  # A = pi R^2
  # R = sqrt(A/pi)
  xlims = c(-12.5, 15012.5)/anisFac
  ylims = c(-8.3335, 5008.4336)
  domainArea = diff(xlims) * diff(ylims)
  sqrt(repArea * domainArea / pi)
}


oldToNewRepA = function(repArea=.01) {
  xlims = c(-12.5, 15012.5)
  ylims = c(-8.3335, 5008.4336)
  domainArea = diff(xlims) * diff(ylims)
  rad = repAreaToDist(repArea)
  pi * (0.5 * rad)^2 / domainArea # turns out to be previous area / 4
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

# Returns the normalization factor of the truth used in the batch case for the 
# given simulation replicate, i.e. if normLgtTruth = fac * (Lgttruth - meanLgtTruth), 
# returns fac, meanLgtTruth, and normLgtTruth
getNormFac = function(repI=NULL, seismicDat=NULL, truthDat=NULL, indepDat=NULL, 
                      goodCoords=NULL, subsampled=TRUE, truthFacOnly=FALSE, takeLogit=TRUE, 
                      getResidSD=FALSE) {
  
  if(is.null(seismicDat)) {
    # seismic data
    out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedPred_", repI, ".txt"), force01=TRUE)
    seismicDat = out$surfFrame
  }
  
  if(is.null(truthDat)) {
    # truth
    out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", repI, ".txt"), force01=TRUE)
    truthDat = out$surfFrame
  }
  
  if(is.null(indepDat)) {
    # indep
    otherRepI = ((repI + 1) %% 100) + 1
    out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", otherRepI, ".txt"), force01=TRUE)
    indepDat = out$surfFrame
  }
  
  if(!subsampled) {
    # subsample
    
    if(is.null(goodCoords)) {
      goodCoords = subsampleSimStudyGrid(seismicDat)
    }
    seismicDat = seismicDat[goodCoords,]
    truthDat = truthDat[goodCoords,]
    indepDat = indepDat[goodCoords,]
  }
  
  if(truthFacOnly) {
    # standardize truth on logit scale
    truthDatStd = truthDat
    if(takeLogit) {
      truthDatStd[,3] = logit(truthDatStd[,3])
    }
    truthSD = sd(truthDatStd[,3])
    
    1/truthSD
  } else {
    # standardize seismic, truth, and indep data on logit scale
    seismicDatStd = seismicDat
    truthDatStd = truthDat
    indepDatStd = indepDat
    if(takeLogit) {
      seismicDatStd[,3] = logit(seismicDatStd[,3])
      truthDatStd[,3] = logit(truthDatStd[,3])
      indepDatStd[,3] = logit(indepDatStd[,3])
    }
    if(getResidSD) {
      residSD = sd(resid(lm(truthDatStd[,3] ~ seismicDatStd[,3])))
    } else {
      residSD = NULL
    }
    seismicMean = mean(seismicDatStd[,3])
    truthMean = mean(truthDatStd[,3])
    indepMean = mean(indepDatStd[,3])
    seismicSD = sd(seismicDatStd[,3])
    truthSD = sd(truthDatStd[,3])
    indepSD = sd(indepDatStd[,3])
    seismicDatStd[,3] = (seismicDatStd[,3] - seismicMean)/seismicSD
    truthDatStd[,3] = (truthDatStd[,3] - truthMean)/truthSD
    indepDatStd[,3] = (indepDatStd[,3] - indepMean)/indepSD
    
    list(seismicDatStd=seismicDatStd, truthDatStd=truthDatStd, indepDatStd=indepDatStd, 
         seismicMean=seismicMean, truthMean=truthMean, indepMean=indepMean, 
         seismicSD=seismicSD, truthSD=truthSD, indepSD=indepSD, residSD=residSD)
  }
}

# returns all normalizing factors for the truth (1/sigma_eta)
getAllNormFacs = function(takeLogit=TRUE, val=NULL, ...) {
  
  allFacs = numeric(100)
  for(i in 1:100) {
    print(paste0("i = ", i, "/100"))
    if(is.null(val)) {
      allFacs[i] = getNormFac(repI=i, seismicDat=NULL, truthDat=NULL, indepDat=NULL, 
                              goodCoords=NULL, subsampled=FALSE, truthFacOnly=TRUE, 
                              takeLogit=takeLogit, ...)
    } else {
      tmp = getNormFac(repI=i, seismicDat=NULL, truthDat=NULL, indepDat=NULL, 
                              goodCoords=NULL, subsampled=FALSE, truthFacOnly=FALSE, 
                              takeLogit=takeLogit, ...)
      allFacs[i] = tmp[[val]]
    }
    
  }
  
  allFacs
}
# mean(getAllNormFacs())
# [1] 2.318985

getAllLMbetas = function() {
  
  allBetas = numeric(100)
  allResidVars = numeric(100)
  for(i in 1:100) {
    print(paste0("i = ", i, "/100"))
    repI = i
    
    # seismic data
    out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedPred_", repI, ".txt"), force01=TRUE)
    seismicDat = out$surfFrame
    
    # truth
    out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", repI, ".txt"), force01=TRUE)
    truthDat = out$surfFrame
    
    
    goodCoords = subsampleSimStudyGrid(seismicDat)
    seismicDat = seismicDat[goodCoords,]
    truthDat = truthDat[goodCoords,]
    
    mod = lm(I(logit(truthDat[,3])) ~ I(logit(seismicDat[,3])))
    
    allBetas[i] = coef(mod)[2]
    allResidVars[i] = var(resid(mod))
  }
  
  list(allBetas=allBetas, allResidVars=allResidVars)
}

# Main simulation study functions -----

# NOTE:
# repI: ID of truth replicate. 100 total.
# sampleParI: ID of parameters related to well sampling. 4 total. Depends on:
#   repelAreaProp (1-2)
#   propVarCase (1-2)
# wellDatI: ID of full well dataset. 250 total. Depends on:
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
    n = c(20, 40, 60, 250)
    propVarCase = c("realistic", "diggle", "cluster", "realNoClust", "seismic", "uniform")
    prefPar = c(1, 2)
    repelAreaProp = c(0, 0.004, 0.02)
  } else {
    # n = c(10, 20, 30)
    n = c(20, 40, 60)
    propVarCase = c("spde", "self")
    prefPar = c(2)
    repelAreaProp = c(0, 0.02)
  }
  fitModFunI = 1:length(getFitModFuns())
  
  if(adaptScen != "batch") {
    fitModFunI = fitModFunI[-5]
  }
  
  
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
  # - don't consider SPDED model in uniform case, since sample probs are uniform
  
  removeBadRows = function(tab, modelTab=FALSE) {
    # first remove extra beta value for uniform case
    badBetas = (tab$propVarCase == "uniform") & (tab$prefPar == 2)
    if(!is.null(tab$n) && (250 %in% n) && ("n" %in% names(tab))) {
      badBetas = badBetas | ((tab$n == 250) & ((tab$propVarCase != "uniform") & (tab$prefPar != 2)))
    }
    
    # then remove bad values of repelAreaProp
    if(!is.null(tab$n)) {
      # badRepArea = (tab$repelAreaProp * tab$n > 0.5) | ((tab$repelAreaProp != 0) & tab$n > 200)
      badRepArea = (tab$repelAreaProp * tab$n > 1.2) # under manuscript parameterization, equivalent to > 0.3
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
      } else {
        # in uniform case, logit probabilities are all identical, so don't include SPDED model then
        goodSampMod[tab$propVarCase == "uniform"] = 
          tab$fitModFunI[tab$propVarCase == "uniform"] != 5
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
  
  # for well sampling, sample the max n necessary. We can subset later
  wellDatCombs$n = max(n)
  getMaxN = function(i, tab) {
    repAreaProp = tab$repelAreaProp[i]
    prefPar = tab$prefPar[i]
    propVarCase = tab$propVarCase[i]
    
    # sample max number of wells below repulsion limit. (previously: Also in big n case use no repulsion)
    out = max(n[n * repAreaProp <= 1.2])
    # if((out == 250) && (repAreaProp != 0)) {
    #   out = n[length(n)-1]
    # }
    
    # sample max number of wells below prefPar limit
    if((out == 250) && (prefPar != 2) && (propVarCase != "uniform")) {
      out = n[length(n)-1]
    }
    out
  }
  # wellDatCombs$n[wellDatCombs$repelAreaProp * wellDatCombs$n > 0.5] = n[length(n)-1]
  wellDatCombs$n = sapply(1:nrow(wellDatCombs), getMaxN, tab=wellDatCombs)
  
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
                                  regenData=FALSE, verbose=FALSE, batchSize=5) {
  
  tmplogfile <- paste0("savedOutput/simStudy/well_", adaptScen, "_", i, "_tmp.txt")
  sink(tmplogfile)
  cat("Started wells job i =", i, ":\n")
  sink()
  
  tryCatch(simStudyWellSampler(i, adaptScen, regenData, verbose, batchSize), 
           error = function(e) {
             logfile <- paste0("savedOutput/simStudy/well_", adaptScen, "_", i, "_err.txt")
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

# generates sequentially-sampled well data for the simulation study
# anisFac: only used in the adaptive case as model parameter
simStudyWellSampler = function(i=1, adaptScen=c("batch", "adaptPref", "adaptVar"), 
                               regenData=FALSE, verbose=FALSE, batchSize=5, anisFac=3, 
                               doPlot=FALSE, doDebug=FALSE, extraVerbose=doDebug) {
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
  sampModFunI = thisPar$sampModFunI
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
  truthDat = truthDat[goodCoords,]
  
  # set repulsion parameters
  repelType = ifelse(repelDist == 0, "none", "rbf")
  if(repelDist != 0) {
    bwRepel = repelDist
  } else {
    bwRepel = NULL
  }
  
  if(adaptScen != "batch") {
    
    truthFac = getNormFac(seismicDat=1, truthDat=truthDat, indepDat=1, 
                          subsampled=TRUE, goodCoords=goodCoords, truthFacOnly=TRUE)
    
    # use the same preferentiality prior as the realistic non-adaptive case
    prefPar = sqrt(0.25) * prefPar * truthFac
    
    mesh = getSPDEmeshSimStudy(anisFac=anisFac)
    if(sampModFunI < 3) {
      inputList = list(mesh=mesh)
    } else if (sampModFunI == 3) {
      inputList = list(prefMean=prefPar, mesh=mesh)
    } else if (sampModFunI == 4) {
      inputList = list(repelDist=repelDist, prefMean=prefPar, mesh=mesh, anisFac=anisFac)
    }
    
    if(doPlot && (sampModFunI %in% c(3, 4))) {
      # for testing purposes, get the point process results also
      inputList = c(inputList, list(getPPres=TRUE))
    }
    
    modelFitter = getFitModFuns()[[sampModFunI]]
    
    
    # sample the wells
    wellDat = wellSampler(truthDat=truthDat, seismicDat=seismicDat, 
                          modelFitter=modelFitter, nWells=nWells, minN=4,
                          predGrid=cbind(east=seismicDat$east, north=seismicDat$north),
                          transform=logit, invTransform=expit, prefPar=prefPar,
                          samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr,
                          repelType=repelType, bwRepel=bwRepel, batchSize=batchSize, 
                          repelAmount=repelAmount, seed=seed, 
                          int.strategy="eb", strategy="gaussian", 
                          getProbsNoRepOnly=TRUE, fitInputs=inputList, 
                          verbose=extraVerbose, doDebug=doDebug, anisFac=anisFac)
    
    # profvis(wellDat <- wellSampler(truthDat, seismicDat, modelFitter, nWells=5, minN=4,
    #                                predGrid=cbind(east=seismicDat$east, north=seismicDat$north),
    #                                transform=logit, invTransform=expit, prefPar=prefPar,
    #                                samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr,
    #                                repelType=repelType, bwRepel=bwRepel,
    #                                repelAmount=repelAmount, seed=seed,
    #                                int.strategy="eb", strategy="gaussian"))
  } else {
    
    if(propVarCase %in% c("realistic", "diggle", "cluster", "realNoClust", "seismic")) {
      # indep
      otherRepI = ((repI + 1) %% 100) + 1
      out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", otherRepI, ".txt"), force01=TRUE)
      indepDat = out$surfFrame
      
      # subsample
      indepDat = indepDat[goodCoords,]
      
      # standardize (center + normalize)
      normDat = getNormFac(seismicDat=seismicDat, truthDat=truthDat, indepDat=indepDat, 
                           subsampled=TRUE, goodCoords=goodCoords)
      seismicDatStd = normDat$seismicDatStd
      truthDatStd = normDat$truthDatStd
      indepDatStd = normDat$indepDatStd
      
      # combine them into a realistic mix and convert back to [0,1] scale
      sampleDat = seismicDat
      
      if(propVarCase == "realistic") {
        seisProp = 0.5
        truthProp = 0.25
        indProp = 0.25
      } else if(propVarCase == "diggle") {
        seisProp = 0
        truthProp = 1
        indProp = 0
      } else if(propVarCase == "cluster") {
        seisProp = 0
        truthProp = 0
        indProp = 1
      } else if(propVarCase == "realNoClust") {
        seisProp = 0.5
        truthProp = 0.5
        indProp = 0
      } else if(propVarCase == "seismic") {
        seisProp = 1
        truthProp = 0
        indProp = 0
      }
      sampleDat[,3] = sqrt(seisProp) * seismicDatStd[,3] + sqrt(truthProp) * truthDatStd[,3] + sqrt(indProp) * indepDatStd[,3]
      sampleDat[,3] = sampleDat[,3] * (1/sd(sampleDat[,3]))
    } else if(propVarCase == "uniform") {
      # in this case, sampleDat doesn't matter
      sampleDat = seismicDat
      sampleDat[,3] = rep(0, nrow(sampleDat))
    } else {
      stop("unrecognized propVarCase")
    }
    
    # must expit here, since preds argument assumes 0-1 scale
    sampleDat[,3] = expit(sampleDat[,3])
    
    # do batch sampling
    wellDat = wellSampler(nWells=nWells, wellDat=NULL, seismicDat=seismicDat, 
                          predGrid=as.matrix(sampleDat[,1:2]), truthDat=truthDat, 
                          preds=sampleDat[,3], 
                          transform=logit, invTransform=expit, prefPar=prefPar, 
                          samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr, 
                          repelType=repelType, bwRepel=bwRepel, 
                          rbf="uniform", repelAmount=Inf, batchSize=batchSize, 
                          seed=seed, int.strategy="eb", strategy="gaussian", 
                          getProbsNoRepOnly=TRUE, verbose=extraVerbose, 
                          doDebug=doDebug)
    
  }
  
  if(doPlot) {
    
    browser()
    
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


# head(modelFitCombs[(modelFitCombs$fitModFunI == 3) & (modelFitCombs$n == 250) & (modelFitCombs$prefPar == 3) & (modelFitCombs$propVarCase == "diggle") & (modelFitCombs$repelAreaProp == 0),])
#      repelAreaProp propVarCase repI repEffect nuggetVar sigmaSq prefPar fitModFunI   n modelFitI wellDatI sampleParI       seed
# 3520             0      diggle    1       Inf      0.01       1       3          3 250     44023     1804         22 1868918159
# 3950             0      diggle    2       Inf      0.01       1       3          3 250     44073     1819         22 1215976688
# 4329             0      diggle    3       Inf      0.01       1       3          3 250     44123     1834         22 1311296880
# 4706             0      diggle    4       Inf      0.01       1       3          3 250     44173     1849         22 1099855696
# 5107             0      diggle    5       Inf      0.01       1       3          3 250     44223     1864         22 1179314858
# 5494             0      diggle    6       Inf      0.01       1       3          3 250     44273     1879         22 1360835856
# head(modelFitCombs[(modelFitCombs$fitModFunI == 4) & (modelFitCombs$n == 60) & (modelFitCombs$prefPar == 3) & (modelFitCombs$propVarCase == "diggle") & (modelFitCombs$repelAreaProp == 0),])
#      repelAreaProp propVarCase repI repEffect nuggetVar sigmaSq prefPar fitModFunI  n modelFitI wellDatI sampleParI       seed
# 3535             0      diggle    1       Inf      0.01       1       3          4 60     41949     1804         22  329769516
# 3938             0      diggle    2       Inf      0.01       1       3          4 60     42024     1819         22  407592108
# 4325             0      diggle    3       Inf      0.01       1       3          4 60     42099     1834         22 1588549675
# 4710             0      diggle    4       Inf      0.01       1       3          4 60     42174     1849         22   44923799
# 5099             0      diggle    5       Inf      0.01       1       3          4 60     42249     1864         22  239254633
# 5481             0      diggle    6       Inf      0.01       1       3          4 60     42324     1879         22  638595106
# head(modelFitCombs[(modelFitCombs$fitModFunI == 4) & (modelFitCombs$n == 40) & (modelFitCombs$prefPar == 3) & (modelFitCombs$propVarCase == "diggle") & (modelFitCombs$repelAreaProp == 0.01),])
#       repelAreaProp propVarCase repI repEffect nuggetVar sigmaSq prefPar fitModFunI  n modelFitI wellDatI sampleParI       seed
# 41229          0.02      diggle    1       Inf      0.01       1       3          4 40     34451     1806         24  199193277
# 41585          0.02      diggle    2       Inf      0.01       1       3          4 40     34526     1821         24  520229110
# 41906          0.02      diggle    3       Inf      0.01       1       3          4 40     34601     1836         24 1235426228
# 42247          0.02      diggle    4       Inf      0.01       1       3          4 40     34676     1851         24 1101297221
# 42572          0.02      diggle    5       Inf      0.01       1       3          4 40     34751     1866         24   66030971
# 42902          0.02      diggle    6       Inf      0.01       1       3          4 40     34826     1881         24 1176183191
# head(modelFitCombs[(modelFitCombs$fitModFunI == 5) & (modelFitCombs$n == 60) & (modelFitCombs$prefPar == 3) & (modelFitCombs$propVarCase == "diggle") & (modelFitCombs$repelAreaProp == 0),])
#      repelAreaProp propVarCase repI repEffect nuggetVar sigmaSq prefPar fitModFunI  n modelFitI wellDatI sampleParI       seed
# 3530             0      diggle    1       Inf      0.01       1       3          5 60     41964     1804         22 1467909580
# 3948             0      diggle    2       Inf      0.01       1       3          5 60     42039     1819         22  884959045
# 4331             0      diggle    3       Inf      0.01       1       3          5 60     42114     1834         22 1463736114
# 4708             0      diggle    4       Inf      0.01       1       3          5 60     42189     1849         22 1274854046
# 5096             0      diggle    5       Inf      0.01       1       3          5 60     42264     1864         22 1251339718
# 5494             0      diggle    6       Inf      0.01       1       3          5 60     42339     1879         22 2144331088
runSimStudyI = function(i, significance=c(.8, .95), 
                        adaptScen=c("batch", "adaptPref", "adaptVar"), 
                        regenData=FALSE, verbose=FALSE, doPlot=FALSE, anisFac=3) {
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
  seismicDat[,1] = seismicDat[,1]/anisFac
  
  # truth
  out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", repI, ".txt"), force01=TRUE)
  truth = out$surfFrame
  truth[,1] = truth[,1]/anisFac
  
  # subsample
  goodCoords = subsampleSimStudyGrid(seismicDat)
  seismicDat = seismicDat[goodCoords,]
  truth = truth[goodCoords,]
  
  # well data
  wellDatFile = paste0("savedOutput/simStudy/wellDat/wellDat_", adaptScen, "_par", sampleParI, "_rep", repI, ".RData")
  out = load(wellDatFile)
  
  if((fitModFunI == 5) || doPlot) {
    if("logProbsNoRep" %in% names(wellDat)) {
      logProbsNoRep = wellDat$logProbsNoRep
    } else {
      logProbsNoRep = wellDat$allLogProbsNoRep[,1]
    }
    
  }
  
  # subset well data based on n for this run
  if("wellDat" %in% names(wellDat)) {
    wellDat = wellDat$wellDat[1:n,]
  } else {
    wellDat = wellDat[1:n,]
  }
  wellDat[,1] = wellDat[,1]/anisFac
  
  # interpolate truth to well points
  truthWells = bilinearInterp(wellDat[,1:2], truth, transform=logit, invTransform=expit)
  
  # Fit model and calculate scores if need be
  scoresFile = paste0("savedOutput/simStudy/scores/scores_", adaptScen, "_", i, ".RData")
  if(!file.exists(scoresFile) || regenData) {
    
    truthFac = getNormFac(seismicDat=1, truthDat=truth, indepDat=1, 
                          subsampled=TRUE, goodCoords=goodCoords, truthFacOnly=TRUE)
    
    if(propVarCase %in% c("uniform", "cluster", "seismic")) {
      # make sure priors are centered around the truth in this case
      prefPar = 0
    } else if(propVarCase == "realNoClust") {
      prefPar = sqrt(0.5) * prefPar * truthFac
    } else if(propVarCase == "realistic") {
      prefPar = sqrt(0.25) * prefPar * truthFac
    } else if(propVarCase == "diggle") {
      prefPar = prefPar * truthFac
    }
    
    mesh = getSPDEmeshSimStudy(anisFac=anisFac)
    if(fitModFunI < 3) {
      inputList = list(wellDat, seismicDat, mesh=mesh)
    } else if (fitModFunI == 3) {
      inputList = list(wellDat, seismicDat, prefMean=prefPar, mesh=mesh)
    } else if (fitModFunI == 4) {
      repDist = repAreaToDist(repelAreaProp)
      inputList = list(wellDat, seismicDat, repelDist=repDist, prefMean=prefPar, mesh=mesh, anisFac=anisFac)
    }
    else if (fitModFunI == 5) {
      inputList = list(wellDat, seismicDat, logProbsNoRep=logProbsNoRep, mesh=mesh)
    }
    
    if(doPlot && (fitModFunI %in% c(3, 4))) {
      # for testing purposes, get the point process results also
      inputList = c(inputList, list(getPPres=TRUE))
    }
    
    out = do.call("fitModFun", inputList)
    
    if(fitModFunI == 4) {
      predMat = out$predMat.y # doesn't include nugget
      predAggMat = out$pred.yAggMat # doesn't include nugget
      obsMat = out$obsMat.y # doesn't include nugget
    } else {
      predMat = out$predMat # doesn't include nugget
      predAggMat = out$predAggMat # doesn't include nugget?
      obsMat = out$obsMat # doesn't include nugget
    }
    fixedEffectSummary = out$fixedEffectSummary
    parameterSummaryTable = out$parameterSummaryTable
    
    
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
    
    if(doPlot) {
      
      browser()
      
      plotSeisGrid(cbind(out$predCoords, out$xPred[,3]))
      points(out$obsCoords[,1], out$obsCoords[,2])
      
      # plot some figures for testing purposes
      
      gEast = seismicDat$east
      gNorth = seismicDat$north
      gSeismic = seismicDat$seismicEst
      gTruth = truth[,3]
      pEast = wellDat$east
      pNorth = wellDat$north
      pVolFrac = wellDat$volFrac
      
      pSeismic = bilinearInterp(cbind(pEast, pNorth), seismicDat, 
                                transform=logit, invTransform = expit)
      
      if("logProbsNoRep" %in% names(wellDat)) {
        logProbsNoRepWells = wellDat$logProbsNoRep
      }
      
      
      
      preds = ests
      sds = out$predSDs
      predQuants = sapply(1:length(preds), function(i) {
        ecdf(predMat[i,])(gTruth[i])
      })
      predQuants[predQuants > .975] = 1
      predQuants[predQuants < .025] = 0
      
      eastGrid = sort(unique(gEast))
      northGrid = sort(unique(gNorth))
      
      # basic model testing
      summary(lm(I(logit(truth[,3])) ~ I(logit(seismicDat[,3]))))
      # Coefficients:
      #   Estimate Std. Error t value Pr(>|t|)    
      # (Intercept)                1.15260    0.01583   72.79   <2e-16 ***
      #   I(logit(seismicDat[, 3]))  1.73034    0.01493  115.87   <2e-16 ***
      #   ---
      #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # Residual standard error: 0.3249 on 20299 degrees of freedom
      # Multiple R-squared:  0.3981,	Adjusted R-squared:  0.3981 
      # F-statistic: 1.343e+04 on 1 and 20299 DF,  p-value: < 2.2e-16
      summary(lm(logProbsNoRepWells ~ I(logit(wellDat$volFrac)) + I(logit(pSeismic))))
      summary(lm(logProbsNoRep ~ I(logit(gTruth)) + I(logit(gSeismic))))
      
      summary(out$mod)
      
      mean((truth[,3] - seismicDat[,3])^2)
      mean((truth[,3] - preds)^2)
      
      # plotting setup
      ticks = seq(0, 1, by=.2)
      tickLabs = as.character(ticks)
      
      seqCols = function(n) {purpleYellowSeqCols(n)}
      quantCols = function(n) {
        tmp = seqCols(n)
        tmp[n] = "red"
        tmp[1] = "blue"
        tmp
      }
      
      pdf(file=paste0("figures/testSimStudy/testPreds_simStudy_", adaptScen, "_par", 
                      sampleParI, "_rep", repI, ".pdf"), width=8, height=ifelse(anisFac==1,5, 7.5))
      par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
      squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
             xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
             zlim=range(c(preds, gTruth)), asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
      points(pEast, pNorth, cex=.5)
      
      squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
             xlab="Easting", ylab="Northing", main="Seismic Estimate", colScale=seqCols, 
             asp=1, smallplot=c(.83,.87,.25,.8))
      points(pEast, pNorth, cex=.5)
      
      splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
            resetGraphics=FALSE, 
            zlim=range(c(preds, gTruth, pVolFrac)), xlab="Easting", ylab="Northing", main="Well Data", 
            asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
      
      squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
             colScale=seqCols, zlim=range(c(preds, gTruth)), 
             xlab="Easting", ylab="Northing", main="Estimate", 
             asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
      points(pEast, pNorth, cex=.5)
      dev.off()
      
      pdf(file=paste0("figures/testSimStudy/testQuants_simStudy_", adaptScen, "_par", 
                      sampleParI, "_rep", repI, ".pdf"), width=8, height=5)
      par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
      squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
             xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
             asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
      points(pEast, pNorth, cex=.5)
      
      squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
             xlab="Easting", ylab="Northing", main="Seismic Estimate", 
             asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
      points(pEast, pNorth, cex=.5)
      
      splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
            resetGraphics=FALSE, 
            xlab="Easting", ylab="Northing", main="Well Data", 
            asp=1, legend.args=list(smallplot=c(.83,.87,.25,.8)), ticks=ticks, tickLabels=tickLabs)
      
      squilt(gEast, gNorth, predQuants, grid=list(x=eastGrid, y=northGrid), 
             colScale=quantCols, 
             xlab="Easting", ylab="Northing", main="Truth Quantiles", 
             asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
      points(pEast, pNorth, cex=.5)
      dev.off()
      
      pdf(file=paste0("figures/testSimStudy/testLogitSDs_simStudy_", adaptScen, "_par", 
                      sampleParI, "_rep", repI, ".pdf"), width=8, height=5)
      par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
      squilt(gEast, gNorth, logit(gTruth), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
             xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
             asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
      points(pEast, pNorth, cex=.5)
      
      squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
             xlab="Easting", ylab="Northing", main="Seismic Estimate", 
             asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
      points(pEast, pNorth, cex=.5)
      
      splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
            resetGraphics=FALSE, 
            xlab="Easting", ylab="Northing", main="Well Data", 
            asp=1, legend.args=list(smallplot=c(.83,.87,.25,.8)), ticks=ticks, tickLabels=tickLabs)
      
      squilt(gEast, gNorth, predQuants, grid=list(x=eastGrid, y=northGrid), 
             colScale=seqCols, 
             xlab="Easting", ylab="Northing", main="Truth Quantiles", 
             asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
      points(pEast, pNorth, cex=.5)
      dev.off()
      
      if(fitModFunI %in% c(3:4)) {
        lambdas = exp(rowMeans(out$predMat.pp))
        lambdas = lambdas * (1/sum(lambdas))
        probsNoRep = exp(logProbsNoRep)
        probsNoRep = probsNoRep * (1/sum(probsNoRep))
        probLims = range(c(lambdas, probsNoRep))
        
        if(fitModFunI == 3) {
          estV = rowMeans(out$spatialPredMat)
          estFixedPt = rowMeans(out$fixedPredMat)
        } else {
          estV = rowMeans(out$spatialPredMat.y)
          estFixedPt = rowMeans(out$fixedPredMat.y)
        }
        
        lmod = lm(I(logit(truth[,3])) ~ I(logit(seismicDat[,3])))
        summary(lmod)
        trueV = resid(lmod)
        rangeV = range(c(estV, trueV))
        
        trueFixedPt = fitted(lmod)
        
        rangeFixedPt = range(c(trueFixedPt, estFixedPt))
        
        pdf(file=paste0("figures/testSimStudy/testLambdas_simStudy_", adaptScen, "_par", 
                        sampleParI, "_rep", repI, ".pdf"), width=16, height=7.5)
        par(mfrow=c(2, 4), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
        
        # top row
        squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
               xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
               zlim=range(c(preds, gTruth, pVolFrac)), 
               asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
        points(pEast, pNorth, cex=.5)
        
        squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
               colScale=seqCols, zlim=range(c(preds, gTruth, pVolFrac)), 
               xlab="Easting", ylab="Northing", main="Estimated Sand Volume Frac", 
               asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
        points(pEast, pNorth, cex=.5)
        
        squilt(gEast, gNorth, estV, grid=list(x=eastGrid, y=northGrid), 
               colScale=seqCols, zlim=rangeV, 
               xlab="Easting", ylab="Northing", main="Estimated V", 
               asp=1, smallplot=c(.83,.87,.25,.8))
        points(pEast, pNorth, cex=.5)
        
        squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
               xlab="Easting", ylab="Northing", main="Seismic Estimate", colScale=seqCols, 
               asp=1, smallplot=c(.83,.87,.25,.8))
        points(pEast, pNorth, cex=.5)
        
        # bottom row
        squilt(gEast, gNorth, probsNoRep, grid=list(x=eastGrid, y=northGrid), 
               colScale=seqCols, 
               xlab="Easting", ylab="Northing", main="True probs no rep", 
               asp=1, smallplot=c(.83,.87,.25,.8), zlim=probLims)
        points(pEast, pNorth, cex=.5)
        
        squilt(gEast, gNorth, lambdas, grid=list(x=eastGrid, y=northGrid), 
               xlab="Easting", ylab="Northing", main="Sampling Intensity Estimate", colScale=seqCols, 
               asp=1, smallplot=c(.83,.87,.25,.8), zlim=probLims)
        points(pEast, pNorth, cex=.5)
        
        squilt(gEast, gNorth, trueV, grid=list(x=eastGrid, y=northGrid), 
               colScale=seqCols, zlim=rangeV, 
               xlab="Easting", ylab="Northing", main="V", 
               asp=1, smallplot=c(.83,.87,.25,.8))
        points(pEast, pNorth, cex=.5)
        
        splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
              resetGraphics=FALSE, 
              zlim=range(c(preds, gTruth, pVolFrac)), xlab="Easting", ylab="Northing", main="Well Data", 
              asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
        
        dev.off()
        
        
        
        
        pdf(file=paste0("figures/testSimStudy/testLambdas2_simStudy_", adaptScen, "_par", 
                        sampleParI, "_rep", repI, ".pdf"), width=16, height=7.5)
        par(mfrow=c(2, 4), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
        
        # top row
        squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
               xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
               zlim=range(c(preds, gTruth, pVolFrac)), 
               asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
        points(pEast, pNorth, cex=.5)
        
        squilt(gEast, gNorth, probsNoRep, grid=list(x=eastGrid, y=northGrid), 
               colScale=seqCols, 
               xlab="Easting", ylab="Northing", main="True probs no rep", 
               asp=1, smallplot=c(.83,.87,.25,.8), zlim=probLims)
        points(pEast, pNorth, cex=.5)
        
        squilt(gEast, gNorth, trueV, grid=list(x=eastGrid, y=northGrid), 
               colScale=seqCols, zlim=rangeV, 
               xlab="Easting", ylab="Northing", main="V", 
               asp=1, smallplot=c(.83,.87,.25,.8))
        points(pEast, pNorth, cex=.5)
        
        squilt(gEast, gNorth, trueFixedPt, grid=list(x=eastGrid, y=northGrid), 
               colScale=seqCols, zlim=rangeFixedPt, 
               xlab="Easting", ylab="Northing", main="Fixed part", 
               asp=1, smallplot=c(.83,.87,.25,.8))
        points(pEast, pNorth, cex=.5)
        
        # bottom row
        squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
               colScale=seqCols, zlim=range(c(preds, gTruth, pVolFrac)), 
               xlab="Easting", ylab="Northing", main="Estimated Sand Volume Frac", 
               asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
        points(pEast, pNorth, cex=.5)
        
        squilt(gEast, gNorth, lambdas, grid=list(x=eastGrid, y=northGrid), 
               xlab="Easting", ylab="Northing", main="Sampling Intensity Estimate", colScale=seqCols, 
               asp=1, smallplot=c(.83,.87,.25,.8), zlim=probLims)
        points(pEast, pNorth, cex=.5)
        
        squilt(gEast, gNorth, estV, grid=list(x=eastGrid, y=northGrid), 
               colScale=seqCols, zlim=rangeV, 
               xlab="Easting", ylab="Northing", main="Estimated V", 
               asp=1, smallplot=c(.83,.87,.25,.8))
        points(pEast, pNorth, cex=.5)
        
        squilt(gEast, gNorth, estFixedPt, grid=list(x=eastGrid, y=northGrid), 
               xlab="Easting", ylab="Northing", main="Estimated fixed part", colScale=seqCols, 
               asp=1, smallplot=c(.83,.87,.25,.8), zlim=rangeFixedPt)
        points(pEast, pNorth, cex=.5)
        
        dev.off()
      }
      
      browser()
      
    }
    
    # save results
    save(pwScoresMean, pwScoresWorst, aggScores, pwScoresMax, pwScoresMin, 
         corSeisTruthWells, corSeisTruthTrue, varTruth, varSeis, 
         varEst, corEstTruthWells, corEstTruthTrue, totT=totT, 
         fixedEffectSummary=fixedEffectSummary, parameterSummaryTable=parameterSummaryTable, 
         file=scoresFile)
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
                           doPar=TRUE, regenData=FALSE, fitModFunIs=1:5) {
  adaptScen = match.arg(adaptScen)
  
  startT = proc.time()[3]
  
  # load simulation parameters
  adaptScenCap = str_to_title(adaptScen)
  inputListFile = paste0("savedOutput/simStudy/simParList", adaptScenCap, ".RData")
  out = load(inputListFile)
  
  # figure out which parameter sets have repI <= maxRepI
  is = 1:nrow(modelFitCombs)
  is = is[(modelFitCombs$repI <= maxRepI) & (modelFitCombs$fitModFunI %in% fitModFunIs)]
  
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
  
  endT = proc.time()[3]
  totHours = (endT - startT)/(60 * 60)
  
  print(paste0("took ", totHours, " hours, or ", totHours/maxRepI, " hours per replicate"))
}

getAllSeisEstsSimStudy = function(regenData=FALSE, significance=c(.8, .95)) {
  
  for(i in 1:100) {
    print(paste0("generating seismic estimates for iteration ", i, "/100"))
    getSeismicEsts(i, regenData=regenData, significance=significance)
  }
  
  invisible(NULL)
}

# NOTE: currently only handles batch case
showSimStudyRes = function(adaptScen=c("batch", "adaptPref", "adaptVar"), maxRepI=100) {
  
  # set arg defaults ----
  adaptScen = match.arg(adaptScen)
  
  # load simulation parameters and model run parameters
  adaptScenCap = str_to_title(adaptScen)
  inputListFile = paste0("savedOutput/simStudy/simParList", adaptScenCap, ".RData")
  out = load(inputListFile)
  
  # helper functions ----
  
  # modI: the type of model being fit (0 = seismic only)
  # thisParI: the parameter scenario we will aggregate results for
  getModScores = function(modI, thisParI, sampModI=NULL) {
    
    # indices of the scores files we'll need to load and the corresponding 
    # rows of modelFitCombs
    if(modI != 0) {
      if(!is.null(sampModI)) {
        thisLs = (subModelCombs$sampleParI %in% thisParI) & (subModelCombs$fitModFunI == modI) & 
          (subModelCombs$sampModFunI == sampModI)
        if(sampModI == modI) {
          thisLs = thisLs & (subModelCombs$propVarCase == "self") # can choose either spde or self for spde, but only self for others
        }
        thisIs = which(thisLs)
      } else {
        thisIs = which((subModelCombs$sampleParI %in% thisParI) & (subModelCombs$fitModFunI == modI))
      }
      
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
    parEstsAll = c()
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
      #      varEst, corEstTruthWells, corEstTruthTrue, totT=totT, 
      #      fixedEffectSummary=fixedEffectSummary, parameterSummaryTable=parameterSummaryTable, 
      #      file=scoresFile)
      
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
      if(modI != 0) {
        parEstsAll = rbind(parEstsAll, c(fixedEffectSummary[,1], parameterSummaryTable[,1]))
      }
    }
    
    if(modI != 0) {
      colnames(parEstsAll) = c(row.names(fixedEffectSummary), row.names(parameterSummaryTable))
    }
    
    list(pwScoresMeanAll=pwScoresMeanAll, pwScoresWorstAll=pwScoresWorstAll, 
         aggScoresAll=aggScoresAll, pwScoresMaxAll=pwScoresMaxAll, 
         pwScoresMinAll=pwScoresMinAll, corSeisTruthWellsAll=corSeisTruthWellsAll, 
         corSeisTruthTrueAll=corSeisTruthTrueAll, varTruthAll=varTruthAll, 
         varSeisAll=varSeisAll, varEstAll=varEstAll, 
         corEstTruthWellsAll=corEstTruthWellsAll, corEstTruthTrueAll=corEstTruthTrueAll, 
         totTAll=totTAll, repIAll=repIAll, nAll=nAll, parEstsAll=parEstsAll)
  }
  
  # combine scores output from getModScores
  bindModScores = function(modScores1, modScores2) {
    pwScoresMeanAll = rbind(modScores1$pwScoresMeanAll, modScores2$pwScoresMeanAll)
    pwScoresWorstAll = rbind(modScores1$pwScoresWorstAll, modScores2$pwScoresWorstAll)
    aggScoresAll = rbind(modScores1$aggScoresAll, modScores2$aggScoresAll)
    pwScoresMaxAll = rbind(modScores1$pwScoresMaxAll, modScores2$pwScoresMaxAll)
    pwScoresMinAll = rbind(modScores1$pwScoresMinAll, modScores2$pwScoresMinAll)
    corSeisTruthWellsAll = c(modScores1$corSeisTruthWellsAll, modScores2$corSeisTruthWellsAll)
    corSeisTruthTrueAll = c(modScores1$corSeisTruthTrueAll, modScores2$corSeisTruthTrueAll)
    varTruthAll = c(modScores1$varTruthAll, modScores2$varTruthAll)
    varSeisAll = c(modScores1$varSeisAll, modScores2$varSeisAll)
    varEstAll = c(modScores1$varEstAll, modScores2$varEstAll)
    corEstTruthWellsAll = c(modScores1$corEstTruthWellsAll, modScores2$corEstTruthWellsAll)
    corEstTruthTrueAll = c(modScores1$corEstTruthTrueAll, modScores2$corEstTruthTrueAll)
    totTAll = c(modScores1$totTAll, modScores2$totTAll)
    repIAll = c(modScores1$repIAll, modScores2$repIAll)
    nAll = c(modScores1$nAll, modScores2$nAll)
    parEstsAll = rbind(modScores1$parEstsAll, modScores2$parEstsAll)
    
    list(pwScoresMeanAll=pwScoresMeanAll, pwScoresWorstAll=pwScoresWorstAll, 
         aggScoresAll=aggScoresAll, pwScoresMaxAll=pwScoresMaxAll, 
         pwScoresMinAll=pwScoresMinAll, corSeisTruthWellsAll=corSeisTruthWellsAll, 
         corSeisTruthTrueAll=corSeisTruthTrueAll, varTruthAll=varTruthAll, 
         varSeisAll=varSeisAll, varEstAll=varEstAll, 
         corEstTruthWellsAll=corEstTruthWellsAll, corEstTruthTrueAll=corEstTruthTrueAll, 
         totTAll=totTAll, repIAll=repIAll, nAll=nAll, parEstsAll=parEstsAll)
  }
  
  # modI: the type of model being fit (0 = seismic only)
  # thisParI: the parameter scenario we will aggregate results for
  # otherParVal: fixed value of the other parameter, other than parName
  getModScoresAcrossPar = function(modI, thisN, parName=c("prefPar", "repelAreaProp"), 
                                   fixedParVal, propVarCase=NULL) {
    
    parName = match.arg(parName)
    fixedParName = ifelse(parName == "prefPar", "repelAreaProp", "prefPar")
    
    # get indices of the corresponding parI
    thisSampleParI = sampleParCombs$sampleParI[sampleParCombs[[fixedParName]] == fixedParVal]
    
    # indices of the scores files we'll need to load and the corresponding 
    # rows of modelFitCombs
    if(modI != 0) {
      thisLs = (subModelCombs$sampleParI %in% thisSampleParI) & 
        (subModelCombs$fitModFunI == modI) & 
        (subModelCombs$n == thisN)
      if(!is.null(propVarCase)) {
        
        if(parName == "prefPar") {
          # uniform is equivalent to prefPar == 0 regardless of the prepVarCase
          thisLs = thisLs & (subModelCombs$propVarCase %in% c("uniform", propVarCase))
        } else {
          thisLs = thisLs & (subModelCombs$propVarCase == propVarCase)
        }
      }
      thisIs = which(thisLs)
      
      # make cases go from small to large value of parName
      thisParVals = subModelCombs[[parName]][thisIs]
      ordI = order(thisParVals)
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
    parEstsAll = c()
    otherParAll = c()
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
      #      varEst, corEstTruthWells, corEstTruthTrue, totT=totT, 
      #      fixedEffectSummary=fixedEffectSummary, parameterSummaryTable=parameterSummaryTable, 
      #      file=scoresFile)
      
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
      parEstsAll = rbind(parEstsAll, c(rowMeans(fixedEffectSummary), c(parameterSummaryTable[,1])))
      otherParAll = c(otherParAll, thesePar[[parName]])
    }
    
    out = list(pwScoresMeanAll=pwScoresMeanAll, pwScoresWorstAll=pwScoresWorstAll, 
               aggScoresAll=aggScoresAll, pwScoresMaxAll=pwScoresMaxAll, 
               pwScoresMinAll=pwScoresMinAll, corSeisTruthWellsAll=corSeisTruthWellsAll, 
               corSeisTruthTrueAll=corSeisTruthTrueAll, varTruthAll=varTruthAll, 
               varSeisAll=varSeisAll, varEstAll=varEstAll, 
               corEstTruthWellsAll=corEstTruthWellsAll, corEstTruthTrueAll=corEstTruthTrueAll, 
               totTAll=totTAll, repIAll=repIAll, parEstsAll=parEstsAll, otherParAll=otherParAll)
  }
  
  collectScoreTab = function(seisScores, spdeScores, spdeKernScores, diggleScores, watsonScores, designScores=NULL, 
                             type=c("agg", "max", "min", "mean", "worst", "par"), 
                             varyN=TRUE, varyParName=NULL, adaptType=c("spde", "self")) {
    
    type = match.arg(type)
    adaptType = match.arg(adaptType)
    keepMod = rep(TRUE, 6) #  in the case type is a parameter, not all models have the parameter
    if(type == "agg") {
      typeName = "aggScoresAll"
    } else if(type == "max") {
      typeName = "pwScoresMaxAll"
    } else if(type == "min") {
      typeName = "pwScoresMinAll"
    } else if(type == "mean") {
      typeName = "pwScoresMeanAll"
    } else if(type == "worst") {
      typeName = "pwScoresWorstAll"
    } else if(type == "par") {
      typeName = "parEstsAll"
    } else {
      stop("bad type")
    }
    
    # collect aggregate scores for each model (repeat seismic estimates for each  
    # value of n we are considering)
    if(type != "par") {
      thisSeis = seisScores[[typeName]]
      thisSeis = matrix(rep(as.matrix(t(thisSeis)), length(spdeScores$totTAll)/length(seisScores$totTAll)), ncol=ncol(thisSeis), byrow = TRUE)
      thisSeis = as.data.frame(thisSeis)
      names(thisSeis) = names(seisScores[[typeName]])
    } else {
      thisSeis = NULL
    }
    
    thisSPDE = spdeScores[[typeName]]
    thisKern = spdeKernScores[[typeName]]
    thisDiggle = diggleScores[[typeName]]
    thisWatson = watsonScores[[typeName]]
    if(!is.null(designScores)) {
      thisDesign = designScores[[typeName]]
    }
    
    # if type == "par", make sure the "prefPar" variable names don't collide
    if(type == "par") {
      colnames(thisDiggle)[colnames(thisDiggle) == "prefPar"] = "pref"
      colnames(thisWatson)[colnames(thisWatson) == "prefPar"] = "pref"
    }
    
    # add in info on n (or the other parameter), model
    if(varyN) {
      addedVar = spdeScores[["nAll"]] # doesn't matter which model, just pick one
      varyParName = "n"
      
      if(!is.null(designScores)) {
        addedVarDes = designScores[["nAll"]]
      }
    } else {
      
      if(is.null(varyParName)) {
        stop("if !varyN, must provide varyParName")
      }
      
      addedVar = spdeScores[["otherParAll"]]
      if(!is.null(designScores)) {
        addedVarDes = designScores[["otherParAll"]]
      }
    }
    
    if(adaptScen == "batch") {
      thisSPDE = cbind(Model="SPDE", n=addedVar, thisSPDE)
      thisKern = cbind(Model="SPDEK", n=addedVar, thisKern)
      thisDiggle = cbind(Model="Diggle", n=addedVar, thisDiggle)
      thisWatson = cbind(Model="Watson", n=addedVar, thisWatson)
      if(!is.null(designScores)) {
        thisDesign = cbind(Model="SPDED", n=addedVarDes, thisDesign)
      }
    } else {
      
      if(adaptType == "spde") {
        thisSPDE = cbind(Model="SPDE->SPDE", n=addedVar, thisSPDE)
        thisKern = cbind(Model="SPDE->SPDEK", n=addedVar, thisKern)
        thisDiggle = cbind(Model="SPDE->Diggle", n=addedVar, thisDiggle)
        thisWatson = cbind(Model="SPDE->Watson", n=addedVar, thisWatson)
      } else if(adaptType == "self") {
        thisSPDE = cbind(Model="SPDE->SPDE", n=addedVar, thisSPDE)
        thisKern = cbind(Model="SPDEK->SPDEK", n=addedVar, thisKern)
        thisDiggle = cbind(Model="Diggle->Diggle", n=addedVar, thisDiggle)
        thisWatson = cbind(Model="Watson->Watson", n=addedVar, thisWatson)
      } else {
        stop("not possible")
      }
      
    }
    
    if(!is.null(thisSeis)) {
      if(nrow(thisSeis) == 0) {
        thisSeis = cbind(Model=character(0), n=numeric(0), thisSeis)
      } else {
        thisSeis = cbind(Model="Seismic", n=addedVar, thisSeis)
      }
      
      colnames(thisSeis)[2] = varyParName
    }
    colnames(thisSPDE)[2] = varyParName
    colnames(thisKern)[2] = varyParName
    colnames(thisDiggle)[2] = varyParName
    colnames(thisWatson)[2] = varyParName
    if(!is.null(designScores)) {
      # thisDesign = cbind(Model="SPDE + design", n=addedVar, thisDesign)
      colnames(thisDesign)[2] = varyParName
    } else {
      thisDesign = NULL
    }
    
    # combine results from all models into a single table
    if(type == "par") {
      
      # set names of parameters to be interprettable, particularly fixed effects
      colnames(thisSPDE)[colnames(thisSPDE) == "X2"] = "seismic_y"
      colnames(thisKern)[colnames(thisKern) == "X2"] = "seismic_y"
      colnames(thisKern)[colnames(thisKern) == "X3"] = "design"
      colnames(thisDiggle)[colnames(thisDiggle) == "X_y"] = "seismic_y"
      colnames(thisDiggle)[colnames(thisDiggle) == "X_pp"] = "seismic_p"
      colnames(thisWatson)[colnames(thisWatson) == "X.y2"] = "seismic_y"
      colnames(thisWatson)[colnames(thisWatson) == "X.y1"] = "X1"
      colnames(thisWatson)[colnames(thisWatson) == "X.pp2"] = "seismic_p"
      colnames(thisWatson)[colnames(thisWatson) == "X.pp1"] = "X1.pp"
      if(!is.null(designScores)) {
        colnames(thisDesign)[colnames(thisDesign) == "X2"] = "seismic_y"
        colnames(thisDesign)[colnames(thisDesign) == "X3"] = "design"
      }
      
      
      # in this case buffer tables with NAs if models don't have given parameters
      if(is.null(thisDesign) ) {
        thisDesign = thisSPDE[numeric(0),]
      }
      if(is.null(thisSeis)) {
        thisSeis = thisSPDE[numeric(0),]
      }
      # Reduce(function(x, y) merge(x, y, all = TRUE), 
      #        list(thisSeis, thisSPDE, thisKern, thisDiggle, thisWatson, thisDesign))
      dplyr::bind_rows(as.data.frame(thisSeis), as.data.frame(thisSPDE), as.data.frame(thisKern), as.data.frame(thisDiggle), as.data.frame(thisWatson), as.data.frame(thisDesign))
    } else {
      rbind(thisSeis, thisSPDE, thisKern, thisDiggle, thisWatson, thisDesign)
    }
  }
  
  
  mean_se <- function(x) {
    x <- na.omit(x)
    if (length(x) == 0) return(c(y = NA, ymin = NA, ymax = NA))
    m <- mean(x)
    se <- sd(x) / sqrt(length(x))
    return(c(y = m, ymin = m - qnorm(.975)*se, ymax = m + qnorm(.975)*se))
  }
  
  
  makeBoxplotsVsN = function(type=c("agg", "max", "min", "mean", "worst", "par"), 
                             seisScores, spdeScores, spdeKernScores, diggleScores, 
                             watsonScores, designScores=NULL, 
                             spdeKernScoresSelf=NULL, diggleScoresSelf=NULL, 
                             watsonScoresSelf=NULL) {
    type = match.arg(type)
    if(type == "agg") {
      typeName = "Aggregate"
    } else if(type == "max") {
      typeName = "Max"
    } else if(type == "min") {
      typeName = "Min"
    } else if(type == "mean") {
      typeName = "Mean"
    } else if(type == "worst") {
      typeName = "Worst"
    } else if(type == "par") {
      typeName = "Estimated"
    } else {
      stop("bad type")
    }
    
    if(adaptScen != "batch") {
      # in this case, get scores with SPDE model sampling, self sampling, and all together
      # remember to remove overlap between the two cases for the SPDE and Seismic models
      tab = collectScoreTab(seisScores, spdeScores, spdeKernScores, diggleScores, watsonScores, 
                            type=type, adaptType="spde")
      tabSelf = collectScoreTab(seisScores, spdeScoresSelf, spdeKernScoresSelf, diggleScoresSelf, watsonScoresSelf, 
                                type=type, adaptType="self")
      tempTab = tabSelf
      tempTab = tempTab[!(tempTab$Model %in% c("SPDE->SPDE", "Seismic")),]
      tabComb = rbind(tab, tempTab)
    } else {
      tab = collectScoreTab(seisScores, spdeScores, spdeKernScores, diggleScores, watsonScores, designScores, 
                            type=type)
    }
    
    
    # plots for each score
    makeThisScorePlots = function(thisTab, adaptType=c("spde", "self", "comb")) {
      adaptType = match.arg(adaptType)
      
      # Ensure the number of colors matches the number of unique models
      # if(adaptScen == "batch") {
      #   modCols = c("grey", "turquoise1", "blue", "purple", "maroon2", "seagreen")
      # } else {
      #   modCols = c("grey", "turquoise1", "blue", "purple", "maroon2")
      #   modColsSelf = c("grey", "turquoise1", "steelblue1", "violet", "palevioletred1")
      #   modColsComb = c("grey", "turquoise1", "blue", "purple", "maroon2", "steelblue1", "violet", "palevioletred1")
      #   # modColsSelf = c("skyblue", "mediumorchid1", "palevioletred1")
      # }
      unique_models = unique(thisTab$Model)
      if (adaptScen == "batch") {
        baseCols = modCols
      } else {
        if (adaptType == "spde") {
          baseCols = modCols
        } else if (adaptType == "self") {
          baseCols = modColsSelf
        } else {
          baseCols = modColsComb
        }
      }
      
      # restrict to models actually present in thisTab, 
      # but preserve the order from baseCols
      presentModels = intersect(names(baseCols), as.character(thisTab$Model))
      thisModCols = baseCols[presentModels]
      
      # set factor levels of Model in the same order as colors
      thisTab$Model = factor(thisTab$Model, levels = names(thisModCols))
      
      # pdf(paste0("figures/simStudy/", fileRoot, "/", type, thisVar, "_", fileRoot, ".pdf"), width=5, height=5)
      adaptFRoot = ""
      if(adaptScen != "batch") {
        adaptFRoot = paste0("_", adaptType)
      }
      
      thisTab = thisTab[!is.na(thisTab[[thisVar]]),]
      thisTab[[thisVar]] = as.numeric(thisTab[[thisVar]])
      levelsN = sort(unique(as.numeric(thisTab$n)))
      thisTab$n = factor(thisTab$n, levels = as.character(levelsN))
      
      # fname = paste0("figures/simStudy/", thisDirRoot, "/", fileSubRoot, "/", fileRoot, "/", type, "_", fileRoot, adaptFRoot, "_", thisVar, ".pdf")
      fname = paste0("figures/simStudy/", thisDirRoot, "/", fileRoot, "/", type, "_", fileRoot, adaptFRoot, "_", thisVar, ".pdf")
      pdf(fname, width=5, height=5)
      
      # Create the boxplot
      p = ggplot(thisTab, aes(x = factor(n), y = .data[[thisVar]], fill = Model))
      if(!(thisVar %in% c("Coverage80", "Coverage95"))) {
        p = p + geom_boxplot()
      }
      p = p +
        # stat_summary(fun = mean, geom = "point", shape = 20, size = 1, color = "black", 
        #              position = position_dodge(width = 0.75)) +
        stat_summary(fun = mean, geom = "point", shape = 21, size = 2, 
                     color = "black", aes(fill = Model),
                     position = position_dodge(width = 0.75)) + 
        stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black", 
                     position = position_dodge(width = 0.75)) +
        # scale_fill_manual(values = setNames(thisModCols[1:length(unique_models)], unique_models)) +
        scale_fill_manual(values = thisModCols) + 
        labs(
          # title = paste0(myTitleCase(thisVar), " vs. n (", typeName, ")"),
          title = paste0(myTitleCase(thisVar), " vs. n"),
          x = "n",
          y = myTitleCase(thisVar),
          fill = "Model"
        ) +
        theme_minimal()
      if(!(thisVar %in% c("Bias", "Var", "Width80", "Width95", "Coverage80", "Coverage95", 
                          "pref", "seismic_y", "seismic_p", "design"))) {
        # seismic estimates have 0 variance
        p = p + scale_y_log10()
      }
      if(grepl("Coverage", thisVar)) {
        cvg = as.numeric(substr(thisVar, nchar(thisVar)-1, nchar(thisVar))) / 100
        p = p + geom_hline(yintercept = cvg, color = "darkgrey", linetype = "dashed") 
      }
      if(thisVar %in% c()) {
        cvg = as.numeric(substr(thisVar, nchar(thisVar)-1, nchar(thisVar))) / 100
        p = p + geom_hline(yintercept = cvg, color = "darkgrey", linetype = "dashed") 
      }
      
      print(p)
      
      dev.off()
      
    }
    
    if(type != "par") {
      varNames = names(spdeScores$aggScoresAll)
    } else {
      # varNames = names(spdeScores$parEstsAll)
      varNames = c("pref", "spatialRange", "spatialVar", "errorVar", "seismic_y", "seismic_p", "design")
    }
    
    for(j in 1:length(varNames)) {
      thisVar = varNames[j]
      
      # keepMod = rep(TRUE, 6)
      # if(thisVar == "prefPar") {
      #   keepMod[c(1, 2, 3, 6)] = FALSE
      # } else if(thisVar == "spatialRange") {
      #   keepMod[1] = FALSE
      # } else if(thisVar == "sptialVar") {
      #   keepMod[1] = FALSE
      # } else if(thisVar == "errorVar") {
      #   keepMod[1] = FALSE
      # } else if(thisVar == "seismic") {
      #   keepMod[1] = FALSE
      # } else if(thisVar == "design") {
      #   keepMod[c(1, 2, 4, 5)] = FALSE
      # } 
      # allMods = c("Seismic", "SPDE", "SPDEK", "Diggle", "Watson", "SPDED")
      # keptsMods = allMods[keepMod]
      
      # if(!dir.exists(paste0("figures/simStudy/", thisDirRoot, "/", fileSubRoot))) {
      #   dir.create(paste0("figures/simStudy/", thisDirRoot, "/", fileSubRoot))
      # }
      # if(!dir.exists(paste0("figures/simStudy/", thisDirRoot, "/", fileSubRoot, "/", fileRoot))) {
      #   dir.create(paste0("figures/simStudy/", thisDirRoot, "/", fileSubRoot, "/", fileRoot))
      # }
      if(!dir.exists(paste0("figures/simStudy/", thisDirRoot, "/", fileRoot))) {
        dir.create(paste0("figures/simStudy/", thisDirRoot, "/", fileRoot))
      }
      
      
      if(adaptScen == "batch") {
        
        # tab = tab[tab$Model %in% keptsMods]
        
        makeThisScorePlots(tab)
      } else {
        
        # # 'Model' variable in adaptive tables is of form [sampleMod->fitMod], 
        # # and we keep only rows from specific fitted models
        # afterArrow <- sub(".*->", "", tab$Model)
        # tab = tab[afterArrow %in% keptsMods,]
        # 
        # afterArrow <- sub(".*->", "", tabSelf$Model)
        # tabSelf = tabSelf[afterArrow %in% keptsMods,]
        # 
        # afterArrow <- sub(".*->", "", tabComb$Model)
        # tabComb = tabComb[afterArrow %in% keptsMods,]
        
        makeThisScorePlots(tab)
        makeThisScorePlots(tabSelf, adaptType="self")
        makeThisScorePlots(tabComb, adaptType="comb")
      }
      
    }
  }
  
  makeBoxplotsVsPar = function(type=c("agg", "max", "min", "mean", "worst", "par"), 
                               thisN, parName=c("prefPar", "repelAreaProp"), 
                               fixedParVal, propVarCase) {
    type = match.arg(type)
    if(type == "agg") {
      typeName = "Aggregate"
    } else if(type == "max") {
      typeName = "Max"
    } else if(type == "min") {
      typeName = "Min"
    } else if(type == "mean") {
      typeName = "Mean"
    } else if(type == "worst") {
      typeName = "Worst"
    } else if(type == "par") {
      typeName = "Estimated"
    } else {
      stop("bad type")
    }
    
    fixedParName = ifelse(parName == "prefPar", "repelAreaProp", "prefPar")
    
    
    if(parName == "prefPar") {
      thisFileRoot = paste0(propVarCase, "_prefParAll_repelAreaProp", fixedParVal, "_n", thisN, "_", adaptScen)
    } else {
      thisFileRoot = paste0(propVarCase, "_repelAreaPropAll_prefPar", fixedParVal, "_n", thisN, "_", adaptScen)
    }
    
    
    # if(!dir.exists(paste0("figures/simStudy/", thisDirRoot, "/", fileSubRoot))) {
    #   dir.create(paste0("figures/simStudy/", thisDirRoot, "/", fileSubRoot))
    # }
    # if(!dir.exists(paste0("figures/simStudy/", thisDirRoot, "/", fileSubRoot, "/", fileRoot))) {
    #   dir.create(paste0("figures/simStudy/", thisDirRoot, "/", fileSubRoot, "/", fileRoot))
    # }
    
    # unlike makeBoxPlotsVsN, need to make most of the model score tables here. 
    # Don't need to mess with seisScores, since they don't depend on anything
    
    spdeScores = getModScoresAcrossPar(1, thisN=thisN, parName=parName, 
                                       fixedParVal=fixedParVal, propVarCase=propVarCase)
    spdeKernScores = getModScoresAcrossPar(2, thisN=thisN, parName=parName, 
                                           fixedParVal=fixedParVal, propVarCase=propVarCase)
    diggleScores = getModScoresAcrossPar(3, thisN=thisN, parName=parName, 
                                         fixedParVal=fixedParVal, propVarCase=propVarCase)
    watsonScores = getModScoresAcrossPar(4, thisN=thisN, parName=parName, 
                                         fixedParVal=fixedParVal, propVarCase=propVarCase)
    if((adaptScen == "batch") && (propVarCase != "uniform")) {
      designScores = getModScoresAcrossPar(5, thisN=thisN, parName=parName, 
                                           fixedParVal=fixedParVal, propVarCase=propVarCase)
    } else {
      designScores = NULL
    }
    
    tab = collectScoreTab(seisScores, spdeScores, spdeKernScores, diggleScores, watsonScores, designScores, 
                          type=type, varyN=FALSE, varyParName=parName)
    
    # Ensure the number of colors matches the number of unique models
    unique_models <- unique(tab$Model)
    if (adaptScen == "batch") {
      baseCols = modCols
    } else {
      if (adaptType == "spde") {
        baseCols = modCols
      } else if (adaptType == "self") {
        baseCols = modColsSelf
      } else {
        baseCols = modColsComb
      }
    }
    
    # restrict to models actually present in thisTab, 
    # but preserve the order from baseCols
    presentModels = intersect(names(baseCols), as.character(tab$Model))
    thisModCols = baseCols[presentModels]
    
    # set factor levels of Model in the same order as colors
    tab$Model = factor(tab$Model, levels = names(thisModCols))
    
    # plot each score
    if(type != "par") {
      varNames = names(spdeScores$aggScoresAll)
    } else {
      varNames = c("pref", "spatialRange", "spatialVar", "errorVar", "seismic_y", "seismic_p", "design")
    }
    
    # reparameterize repelAreaProp (divide by 4)
    if(parName == "repelAreaProp") {
      tab$repelAreaProp = as.character(as.numeric(tab$repelAreaProp) * 0.25) # faster than call to oldToNewRepA()
    }
    
    for(j in 1:length(varNames)) {
      thisVar = varNames[j]
      
      if(!dir.exists(paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot))) {
        dir.create(paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot))
      }
      
      thisTab = tab
      thisTab = thisTab[!is.na(thisTab[[thisVar]]),]
      thisTab[[thisVar]] = as.numeric(thisTab[[thisVar]])
      
      pdf(paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot, "/", type, "_", thisFileRoot, "_", thisVar, ".pdf"), width=5, height=5)
      
      # Create the plot
      if(parName == "prefPar") {
        p = ggplot(thisTab, aes(x = factor(prefPar), y = .data[[thisVar]], fill = Model))
      } else {
        p = ggplot(thisTab, aes(x = factor(repelAreaProp), y = .data[[thisVar]], fill = Model))
      }
      if(!(thisVar %in% c("Coverage80", "Coverage95"))) {
        p = p + geom_boxplot()
      }
      p = p +
        # stat_summary(fun = mean, geom = "point", shape = 20, size = 1, color = "black", 
        #              position = position_dodge(width = 0.75)) +
        stat_summary(fun = mean, geom = "point", shape = 21, size = 2, 
                     color = "black", aes(fill = Model),
                     position = position_dodge(width = 0.75)) + 
        stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black", 
                     position = position_dodge(width = 0.75)) +
        # scale_fill_manual(values = setNames(modCols[1:length(unique_models)], unique_models)) +
        scale_fill_manual(values = thisModCols) + 
        labs(
          # title = paste0(myTitleCase(thisVar), " vs. ", parName, 
          #                " (", adaptScenCap, ", ", typeName, ", ", 
          #                fixedParName, "=", fixedParVal, ", n=", thisN, ")"),
          title = paste0(myTitleCase(thisVar), " vs. ", parName, 
                         " (", fixedParName, "=", fixedParVal, ", n=", thisN, ")"),
          x = parName,
          y = myTitleCase(thisVar),
          fill = "Model"
        ) +
        theme_minimal()
      if(!(thisVar %in% c("Bias", "Var", "Width80", "Width95", "Coverage80", "Coverage95", 
                          "pref", "seismic_y", "seismic_p", "design"))) {
        # seismic estimates have 0 variance
        p = p + scale_y_log10()
      }
      if(grepl("Coverage", thisVar)) {
        cvg <- as.numeric(substr(thisVar, nchar(thisVar)-1, nchar(thisVar))) / 100
        p = p + geom_hline(yintercept = cvg, color = "darkgrey", linetype = "dashed") 
      }
      
      print(p)
      
      dev.off()
    }
  }
  
  # plots and tables setup ----
  print("plotting...")
  
  # make sure computer knows uniform case has no preferentiality
  modelFitCombs$prefPar[modelFitCombs$propVarCase == "uniform"] = 0
  wellDatCombs$prefPar[wellDatCombs$propVarCase == "uniform"] = 0
  sampleParCombs$prefPar[sampleParCombs$propVarCase == "uniform"] = 0
  
  # figure out which parameter sets have repI <= maxRepI
  is = 1:nrow(modelFitCombs)
  is = is[modelFitCombs$repI <= maxRepI]
  
  if(adaptScen != "batch") {
    # TODO: subset modelFitCombs below (via subModelCombs) based on adaptScen
    # JP: I don't think the above TODO is necessary?
  }
  subModelCombs = modelFitCombs[modelFitCombs$repI <= maxRepI,]
  
  thisDirRoot = paste0(adaptScen, "/")
  if(adaptScen == "batch") {
    modCols = c(Seismic="grey", SPDE="turquoise1", SPDEK="blue", Diggle="purple", Watson="maroon2", SPDED="seagreen")
  } else {
    modCols = c(Seismic="grey", SPDE="turquoise1", SPDEK="blue", Diggle="purple", Watson="maroon2")
    modColsSelf = c("grey", "turquoise1", "steelblue1", "violet", "palevioletred1")
    modColsComb = c("grey", "turquoise1", "blue", "purple", "maroon2", "steelblue1", "violet", "palevioletred1")
    # modColsSelf = c("skyblue", "mediumorchid1", "palevioletred1")
    names(modColsSelf) = c("SPDE->SPDE", "SPDE->SPDEK", "SPDE->Diggle", "SPDE->Watson")
    names(modColsComb) = c("SPDE->SPDE", "SPDEK->SPDEK", "Diggle->Diggle", "Watson->Watson")
  }
  
  pch = c(5, 15:19)
  allTypes = c("agg", "max", "min", "mean", "worst", "par")
  
  if(adaptScen == "batch") {
    nUnique = c(20, 40, 60, 250)
  } else {
    # nUnique = c(10, 20, 30)
    nUnique = c(20, 40, 60)
  }
  
  # if(adaptScen == "batch") {
  #   propVarCases = c("uniform", "realistic")
  # } else {
  #   stop("adaptive scenarios not yet supported")
  # }
  
  seisScores = getModScores(0, thisParI=1) # seismic scores don't depend on thisParI
  
  # boxplots vs n ----
  print("boxplots vs n...")
  
  # for each set of sampling parameters, show results
  if(adaptScen == "batch") {
    for(i in 1:nrow(sampleParCombs)) {
      
      thesePar = sampleParCombs[i,]
      
      sampleParI = thesePar$sampleParI # same as i
      repelAreaProp = thesePar$repelAreaProp
      propVarCase = thesePar$propVarCase
      prefPar = thesePar$prefPar
      repEffect = thesePar$repEffect
      
      fileSubRoot = paste0("i", sampleParI, 
                           "_", propVarCase)
      fileRoot = paste0("i", sampleParI, 
                        "_", propVarCase, 
                        "_repA", repelAreaProp, 
                        "_pref", prefPar, 
                        "_", adaptScen)
      
      # get scores for each model
      spdeScores = getModScores(1, i)
      spdeKernScores = getModScores(2, i)
      diggleScores = getModScores(3, i)
      watsonScores = getModScores(4, i)
      if(propVarCase != "uniform") {
        designScores = getModScores(5, i)
      } else {
        designScores = spdeScores
      }
      
      # score names:
      # pwScoresMeanAll, pwScoresWorstAll, aggScoresAll, pwScoresMaxAll, 
      # pwScoresMinAll, corSeisTruthWellsAll, corSeisTruthTrueAll, varTruthAll, 
      # varSeisAll, varEstAll, corEstTruthWellsAll, corEstTruthTrueAll, 
      # totTAll, repIAll, nAll
      
      for(j in 1:length(allTypes)) {
        
        makeBoxplotsVsN(allTypes[j], seisScores, spdeScores, spdeKernScores, diggleScores, 
                        watsonScores, designScores)
      }
      
    }
  } else {
    # adaptive case
    
    repPars = sort(unique(sampleParCombs$repelAreaProp))
    for(i in 1:length(repPars)) {
      thisRepPar = repPars[i]
      
      fileSubRoot = ""
      fileRoot = paste0(adaptScen, "_repA", thisRepPar)
      
      # SPDE sampling case
      
      # get sampleParI's corresponding to thisRepPar
      thisIs = sampleParCombs$sampleParI[(sampleParCombs$sampModFunI == 1) & 
                                           (sampleParCombs$repelAreaProp == thisRepPar)]
      
      # get scores for each model
      spdeScores = getModScores(1, thisIs, 1)
      spdeKernScores = getModScores(2, thisIs, 1)
      diggleScores = getModScores(3, thisIs, 1)
      watsonScores = getModScores(4, thisIs, 1)
      
      thisIs = sampleParCombs$sampleParI[sampleParCombs$repelAreaProp == thisRepPar]
      
      spdeScoresSelf = spdeScores
      spdeKernScoresSelf = getModScores(2, thisIs, 2)
      diggleScoresSelf = getModScores(3, thisIs, 3)
      watsonScoresSelf = getModScores(4, thisIs, 4)
      
      # spdeScoresComb = bindModScores(spdeScores, spdeScoresSelf)
      # spdeKernScoresComb = bindModScores(spdeKernScores, spdeKernScoresSelf)
      # diggleScoresComb = bindModScores(diggleScores, diggleScoresSelf)
      # watsonScoresComb = bindModScores(watsonScores, watsonScoresSelf)
      
      for(j in 1:length(allTypes)) {
        
        makeBoxplotsVsN(allTypes[j], seisScores, spdeScores, spdeKernScores, diggleScores, 
                        watsonScores, spdeKernScoresSelf=spdeKernScoresSelf, diggleScoresSelf=diggleScoresSelf, 
                        watsonScoresSelf=watsonScoresSelf)
        
      }
      
    }
  }
  
  
  # boxplots vs parameters ----
  print("boxplots vs parameters...")
  
  propVarCases = c("uniform", "diggle", "seismic", "cluster", "realNoClust", "realistic")
  
  if(adaptScen == "batch") {
    
    
    repelAreaPropUnique = sort(unique(sampleParCombs$repelAreaProp))
    
    for(caseI in 1:length(propVarCases)) {
      
      thisPropVarCase = propVarCases[caseI]
      
      
      # wrt prefPar (fix repelAreaProp)
      if(thisPropVarCase != "uniform") {
        
        for(i in 1:length(repelAreaPropUnique)) {
          thisRepelAreaProp = repelAreaPropUnique[i]
          
          thisUniqueNs = nUnique[(thisRepelAreaProp * nUnique) <= 1.2]
          # if(thisRepelAreaProp != 0) {
          #   thisUniqueNs = thisUniqueNs[thisUniqueNs < 200]
          # }
          
          for(j in 1:length(allTypes)) {
            
            for(k in 1:length(thisUniqueNs)) {
              thisN = thisUniqueNs[k]
              
              if(adaptScen == "batch") {
                makeBoxplotsVsPar(allTypes[j], parName="prefPar", thisN=thisN, 
                                  fixedParVal=thisRepelAreaProp, propVarCase=thisPropVarCase)
              } else {
                stop("still need to work out adaptive case")
              }
            }
          }
        }
        
      }
      
      # wrt repelAreaProp (fix prefPar)
      prefParUnique = sort(unique(sampleParCombs$prefPar[sampleParCombs$propVarCase == thisPropVarCase]))
      for(i in 1:length(prefParUnique)) {
        thisPrefPar = prefParUnique[i]
        
        for(j in 1:length(allTypes)) {
          
          for(k in 1:length(nUnique)) {
            thisN = nUnique[k]
            
            if(adaptScen == "batch") {
              if((thisN > 200) && (thisPrefPar < 2)) {
                next
              }
              makeBoxplotsVsPar(allTypes[j], parName="repelAreaProp", thisN=thisN, 
                                fixedParVal=thisPrefPar, propVarCase=thisPropVarCase)
            } else {
              stop("still need to work out adaptive case")
            }
          }
        }
      }
      
    }
    
  }
  
  print("copying figures to manuscript fig dir...")
  # copy only the figures we will actually keep for the manuscript
  copyDirFiltered(srcDir=paste0("figures/simStudy/", thisDirRoot), 
                  dstDir=paste0("~/fig/simStudy/", adaptScen), 
                  includeSubstr = c("agg", "par", "mean.*Coverage95"), excludeSubstr = c("80", "_Var", "RMSE"))
  
  browser()
  
}


showSimStudyRes2 = function(adaptScen = c("batch", "adaptPref", "adaptVar"),
                            maxRepI = 100,
                            regenData = FALSE, 
                            adaptType=c("spde", "self", "comb"), 
                            doBoxN=TRUE, doBoxPhi=TRUE, doBoxRep=TRUE, doScatterGamma=TRUE) {
  adaptScen = match.arg(adaptScen)
  adaptType = match.arg(adaptType)
  adaptScenCap = str_to_title(adaptScen)
  adaptTypeCap = str_to_title(adaptType)
  inputListFile = paste0("savedOutput/simStudy/simParList", adaptScenCap, ".RData")
  load(inputListFile)
  
  # make score/par table ----
  if(adaptScen == "batch") {
    adaptTypeCap = ""
  }
  mergedFile = paste0("savedOutput/simStudy/mergedScores_", adaptScen, adaptTypeCap, ".RData")
  
  if (!file.exists(mergedFile) || regenData) {
    message("Regenerating merged score table...")
    subModelCombs = modelFitCombs[modelFitCombs$repI <= maxRepI, ]
    nTotal = nrow(subModelCombs) + maxRepI
    startTime = Sys.time()
    
    mergedTab = lapply(seq_len(nTotal), function(i) {
      if (i %% 500 == 0 || i == nTotal) {
        elapsed = Sys.time() - startTime
        estTotal = as.numeric(elapsed) / i * nTotal
        estRemaining = estTotal - as.numeric(elapsed)
        message(sprintf("Progress: %d/%d (%.1f%%), Estimated time left: %.1fs", i, nTotal, i/nTotal*100, estRemaining))
      }
      
      if(i <= nrow(subModelCombs)) {
        row = subModelCombs[i, ]
      } else {
        # seismic model cases. Doesn't matter what the main variables are since we will overwrite them later
        i = i - nrow(subModelCombs)
        row = subModelCombs[1, ]
        row$repI = i
        row$fitModFunI = 0
      }
      fitModFunI = row$fitModFunI
      modelFitI = row$modelFitI
      
      
      scoresFile = if (fitModFunI == 0) {
        paste0("savedOutput/simStudy/scores/scores_seismic_rep", i, ".RData")
      } else {
        paste0("savedOutput/simStudy/scores/scores_", adaptScen, "_", modelFitI, ".RData")
      }
      
      if (!file.exists(scoresFile)) return(NULL)
      
      out = load(scoresFile)
      
      if (row$propVarCase == "uniform") row$prefPar = 0
      row$repelAreaProp = row$repelAreaProp / 4 # reparameterize to the new one
      
      if(fitModFunI != 0) {
        fixEffs = setNames(data.frame(t(fixedEffectSummary[, 1])), paste0(rownames(fixedEffectSummary), "_param"))
        parInfo = setNames(data.frame(t(parameterSummaryTable[, 1])), paste0(rownames(parameterSummaryTable), "_param"))
        scoreRow = cbind(
          row,
          setNames(as.data.frame(pwScoresMean), paste0(colnames(pwScoresMean), "_pwMean")),
          setNames(as.data.frame(pwScoresWorst), paste0(colnames(pwScoresWorst), "_pwWorst")),
          setNames(as.data.frame(aggScores), paste0(colnames(aggScores), "_agg")),
          setNames(as.data.frame(pwScoresMax), paste0(colnames(pwScoresMax), "_pwMax")),
          setNames(as.data.frame(pwScoresMin), paste0(colnames(pwScoresMin), "_pwMin")),
          corSeisTruthWells = corSeisTruthWells,
          corSeisTruthTrue = corSeisTruthTrue,
          varTruth = varTruth,
          varSeis = varSeis,
          varEst = varEst,
          corEstTruthWells = corEstTruthWells,
          corEstTruthTrue = corEstTruthTrue,
          totT = totT,
          fixEffs, 
          parInfo
        )
      } else {
        scoreRow = cbind(
          row,
          setNames(as.data.frame(pwScoresMean), paste0(colnames(pwScoresMean), "_pwMean")),
          setNames(as.data.frame(pwScoresWorst), paste0(colnames(pwScoresWorst), "_pwWorst")),
          setNames(as.data.frame(aggScores), paste0(colnames(aggScores), "_agg")),
          setNames(as.data.frame(pwScoresMax), paste0(colnames(pwScoresMax), "_pwMax")),
          setNames(as.data.frame(pwScoresMin), paste0(colnames(pwScoresMin), "_pwMin")),
          corSeisTruthWells = corSeisTruthWells,
          corSeisTruthTrue = corSeisTruthTrue,
          varTruth = varTruth,
          varSeis = varSeis,
          varEst = varEst,
          corEstTruthWells = corEstTruthWells,
          corEstTruthTrue = corEstTruthTrue,
          totT = totT
        )
      }
      
      return(scoreRow)
    })
    
    # mergedTab = do.call(rbind, Filter(Negate(is.null), mergedTab))
    mergedTab = do.call(dplyr::bind_rows, Filter(Negate(is.null), mergedTab))
    
    # rename from old names to names in writeup
    names(mergedTab)[names(mergedTab) == "prefPar_param"] = "gamma_param"
    names(mergedTab)[names(mergedTab) == "prefPar"] = "phi"
    
    save(mergedTab, subModelCombs, file = mergedFile)
  } else {
    message("Loading existing merged score table...")
    load(mergedFile)
  }
  
  # process table ----
  mergedTab$Model = getFitModName(mergedTab$fitModFunI)
  
  # plots and tables setup ----
  print("plotting...")
  
  mean_se <- function(x) {
    x <- na.omit(x)
    if (length(x) == 0) return(c(y = NA, ymin = NA, ymax = NA))
    m <- mean(x)
    se <- sd(x) / sqrt(length(x))
    return(c(y = m, ymin = m - qnorm(.975)*se, ymax = m + qnorm(.975)*se))
  }
  
  # make sure computer knows uniform case has no preferentiality
  modelFitCombs$prefPar[modelFitCombs$propVarCase == "uniform"] = 0
  wellDatCombs$prefPar[wellDatCombs$propVarCase == "uniform"] = 0
  sampleParCombs$prefPar[sampleParCombs$propVarCase == "uniform"] = 0
  
  # figure out which parameter sets have repI <= maxRepI
  is = 1:nrow(modelFitCombs)
  is = is[modelFitCombs$repI <= maxRepI]
  
  if(adaptScen != "batch") {
    # TODO: subset modelFitCombs below (via subModelCombs) based on adaptScen
    # JP: I don't think the above TODO is necessary?
  }
  subModelCombs = modelFitCombs[modelFitCombs$repI <= maxRepI,]
  
  thisDirRoot = paste0(adaptScen, "/")
  if(adaptScen == "batch") {
    modCols = c(Seismic="grey", SPDE="turquoise1", SPDEK="blue", Diggle="purple", Watson="maroon2", SPDED="seagreen")
  } else {
    modCols = c(Seismic="grey", SPDE="turquoise1", SPDEK="blue", Diggle="purple", Watson="maroon2")
    modColsSelf = c("grey", "turquoise1", "steelblue1", "violet", "palevioletred1")
    modColsComb = c("grey", "turquoise1", "blue", "purple", "maroon2", "steelblue1", "violet", "palevioletred1")
    # modColsSelf = c("skyblue", "mediumorchid1", "palevioletred1")
    names(modColsSelf) = c("SPDE->SPDE", "SPDE->SPDEK", "SPDE->Diggle", "SPDE->Watson")
    names(modColsComb) = c("SPDE->SPDE", "SPDEK->SPDEK", "Diggle->Diggle", "Watson->Watson")
  }
  
  pch = c(5, 15:19)
  allTypes = c("agg", "max", "min", "mean", "worst", "par")
  
  if(adaptScen == "batch") {
    nUnique = c(20, 40, 60, 250)
  } else {
    # nUnique = c(10, 20, 30)
    nUnique = c(20, 40, 60)
  }
  
  library(ggplot2)
  library(dplyr)
  
  scoreTypes = c("_agg", "_pwMax", "_pwMin", "_pwMean", "_pwWorst", "_param")
  
  scoreTypesNamed = list(
    agg = "Aggregate",
    max = "Max",
    min = "Min",
    mean = "Mean",
    worst = "Worst", 
    par = "Param"
  )
  scoreTypesNameRoot = list(
    agg = "agg",
    max = "max",
    min = "min",
    mean = "mean",
    worst = "worst", 
    par = "par"
  )
  
  propVarCases = unique(mergedTab$propVarCase)
  
  # Automatically make boxplots of the format we want consistently
  # Inputs:
  # thisTab: data.frame with the information to boxplot
  # parName: name of variable on the horizontal axis
  # scoreCol: column name of variable to plot
  # fixedParNames: names of the fixed variables (with one value in thisTab)
  # fname: name of the file to save the plot to
  # adaptType: if adaptScen isn't batch, this represents which sampling/fitting 
  #            model combinations to plot
  makeBoxplot = function(thisTab, parName, scoreCol, fixedParNames, fname) {
    scoreColName = sub("(_pwMean|_pwWorst|_agg|_pwMax|_pwMin|_param)$", "", scoreCol)
    parTitle = sub("(_pwMean|_pwWorst|_agg|_pwMax|_pwMin|_param)$", "", parName)
    
    # add in seismic data:
    # copy tabSeismic, once per unique combination of n, repelAreaProp, and phi
    nVals = sort(unique(thisTab$n))
    repelVals = sort(unique(thisTab$repelAreaProp))
    phiVals = sort(unique(thisTab$phi))
    
    parCombs = expand.grid(n = nVals, repelAreaProp = repelVals, phi = phiVals)
    
    tabList = lapply(1:nrow(parCombs), function(i) {
      parRow = parCombs[i, ]
      tab = tabSeismic
      tab$n = parRow$n
      tab$repelAreaProp = parRow$repelAreaProp
      tab$phi = parRow$phi
      tab
    })
    tabSeisFull = do.call(rbind, tabList)
    thisTab = rbind(thisTab, tabSeisFull)
    
    # filter out NAs
    thisTab = thisTab[!is.na(thisTab[[scoreCol]]), ]
    thisTab[[scoreCol]] = as.numeric(thisTab[[scoreCol]])
    
    # don't include seismic model for uncertainty related scores/metrics
    if(greplAny(c("Coverage", "IntervalScore", "Width"), scoreColName)) {
      thisTab = thisTab[thisTab$Model != "Seismic",]
    }
    
    unique_models = unique(thisTab$Model)
    if (adaptScen == "batch") {
      baseCols = modCols
    } else {
      if (adaptType == "spde") {
        baseCols = modCols
      } else if (adaptType == "self") {
        baseCols = modColsSelf
      } else {
        baseCols = modColsComb
      }
    }
    presentModels = intersect(names(baseCols), as.character(thisTab$Model))
    thisModCols = baseCols[presentModels]
    thisTab$Model = factor(thisTab$Model, levels = names(thisModCols))
    
    pdf(fname, width = 5, height = 5)
    p = ggplot(thisTab, aes(x = factor(.data[[parName]]), y = .data[[scoreCol]], fill = Model))
    if(!(scoreColName %in% c("Coverage80", "Coverage95"))) {
      p = p + geom_boxplot()
    }
    p = p + stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black",
                         position = position_dodge(width = 0.75)) +
      stat_summary(fun = mean, geom = "point", shape = 21, size = 2,
                   color = "black", aes(fill = Model),
                   position = position_dodge(width = 0.75)) +
      scale_fill_manual(values = thisModCols) +
      # labs(title = paste0(myTitleCase(scoreColName), " vs. ", parName,
      #                     " (", fixedParName, "=", unique(thisTab[[fixedParName]]), ")"),
      #      x = parName, y = myTitleCase(scoreColName), fill = "Model") +
      labs(
        title = paste0(
          myTitleCase(scoreColName), " vs. ", parTitle, " (",
          paste(paste0(fixedParNames, "=", sapply(fixedParNames, function(nm) unique(thisTab[[nm]]))), collapse = ", "),
          ")"),
        x = parTitle, y = myTitleCase(scoreColName), fill = "Model") + 
      theme_minimal()
    
    if (!scoreColName %in% c("Bias", "Var", "Width80", "Width95", "Coverage80", "Coverage95",
                             "gamma", "seismic_y", "seismic_p", "design")) {
      p = p + scale_y_log10()
    }
    if (grepl("Coverage", scoreColName)) {
      cvg = as.numeric(substr(scoreColName, nchar(scoreColName)-1, nchar(scoreColName))) / 100
      p = p + geom_hline(yintercept = cvg, color = "darkgrey", linetype = "dashed")
    }
    print(p)
    dev.off()
    
    invisible(NULL)
  }
  
  # Automatically make pairplots of the format we want consistently
  # Inputs:
  # thisTab: data.frame with the information to boxplot
  # parName: name of variable on the horizontal axis
  # scoreCol: column name of variable to plot
  # fixedParNames: names of the fixed variables (with one value in thisTab)
  # fname: name of the file to save the plot to
  # adaptType: if adaptScen isn't batch, this represents which sampling/fitting 
  #            model combinations to plot
  makeScatterplot = function(thisTab, parName, scoreCol, fixedParNames, fname) {
    scoreColName = sub("(_pwMean|_pwWorst|_agg|_pwMax|_pwMin|_param)$", "", scoreCol)
    parTitle = sub("(_pwMean|_pwWorst|_agg|_pwMax|_pwMin|_param)$", "", parName)
    thisTab = thisTab[!is.na(thisTab[[scoreCol]]), ]
    thisTab = thisTab[!is.na(thisTab[[parName]]), ]
    thisTab[[scoreCol]] = as.numeric(thisTab[[scoreCol]])
    
    unique_models = unique(thisTab$Model)
    if (adaptScen == "batch") {
      baseCols = modCols
    } else {
      if (adaptType == "spde") {
        baseCols = modCols
      } else if (adaptType == "self") {
        baseCols = modColsSelf
      } else {
        baseCols = modColsComb
      }
    }
    presentModels = intersect(names(baseCols), as.character(thisTab$Model))
    thisModCols = baseCols[presentModels]
    thisTab$Model = factor(thisTab$Model, levels = names(thisModCols))
    
    
    pdf(fname, width = 5, height = 5)
    p = ggplot(thisTab, aes(x = .data[[parName]], y = .data[[scoreCol]], color = Model, shape=Model)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = thisModCols) +
      labs(
        title = paste0(
          myTitleCase(scoreColName), " vs. ", parTitle, " (",
          paste(paste0(fixedParNames, "=", sapply(fixedParNames, 
                                                  function(nm) unique(thisTab[[nm]]))), collapse = ", "),
          ")"),
        x = parTitle, y = myTitleCase(scoreColName), color = "Model") +
      theme_minimal()
    
    if (!scoreColName %in% c("Bias", "Var", "Width80", "Width95", "Coverage80", "Coverage95",
                             "gamma", "seismic_y", "seismic_p", "design")) {
      p = p + scale_y_log10()
    }
    if (grepl("Coverage", scoreColName)) {
      cvg = as.numeric(substr(scoreColName, nchar(scoreColName)-1, nchar(scoreColName))) / 100
      p = p + geom_hline(yintercept = cvg, color = "darkgrey", linetype = "dashed")
    }
    print(p)
    dev.off()
    
    
    # pdf(fname, width = 5, height = 5)
    # p = ggplot(thisTab, aes(x = .data[[parName]], y = .data[[scoreCol]], fill = Model)) +
    #   geom_boxplot() +
    #   stat_summary(fun = mean, geom = "point", shape = 21, size = 2,
    #                color = "black", aes(fill = Model),
    #                position = position_dodge(width = 0.75)) +
    #   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black",
    #                position = position_dodge(width = 0.75)) +
    #   scale_fill_manual(values = thisModCols) +
    #   # labs(title = paste0(myTitleCase(scoreColName), " vs. ", parName,
    #   #                     " (", fixedParName, "=", unique(thisTab[[fixedParName]]), ")"),
    #   #      x = parName, y = myTitleCase(scoreColName), fill = "Model") +
    #   labs(
    #     title = paste0(
    #       myTitleCase(scoreColName), " vs. ", parName, " (",
    #       paste(paste0(fixedParNames, "=", sapply(fixedParNames, function(nm) unique(thisTab[[nm]]))), collapse = ", "),
    #       ")"),
    #     x = parName, y = myTitleCase(scoreColName), fill = "Model") + 
    #   theme_minimal()
    # 
    # if (!scoreColName %in% c("Bias", "Var", "Width80", "Width95", "Coverage80", "Coverage95",
    #                          "pref", "seismic_y", "seismic_p", "design")) {
    #   p = p + scale_y_log10()
    # }
    # if (grepl("Coverage", scoreColName)) {
    #   cvg = as.numeric(substr(scoreColName, nchar(scoreColName)-1, nchar(scoreColName))) / 100
    #   p = p + geom_hline(yintercept = cvg, color = "darkgrey", linetype = "dashed")
    # }
    # print(p)
    # dev.off()
    
    invisible(NULL)
  }
  
  getLatexLongtable = function(tab) {
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install the 'dplyr' package.")
    
    require(dplyr)
    
    # get names of scores to put in table
    thisScoreNames = names(tab)[grepl("_agg", names(tab))]
    
    # exclude scores we don't care about
    thisScoreNames = thisScoreNames[!greplAny(c("80", "Var", "MSE", "RMSE"), thisScoreNames)]
    
    thisScoreTitles = sub("(_pwMean|_pwWorst|_agg|_pwMax|_pwMin|_param)$", "", thisScoreNames)
    thisScoreTitles = gsub("IntervalScore", "IS", thisScoreTitles)
    thisScoreTitles = gsub("Coverage", "Cvg", thisScoreTitles)
    
    
    # Define model order
    modelOrder = c("SPDE", "SPDEK", "Diggle", "Watson", "SPDED")
    
    
    # Compute mean and SE for each group
    summaryTab = tab %>%
      group_by(phi, repelAreaProp, n, Model) %>%
      summarise(across(all_of(thisScoreNames), list(
        mean = ~mean(.),
        se = ~sd(.) / sqrt(length(.))
      ), .names = "{.col}{.fn}"), .groups = "drop")
    
    # Start LaTeX longtable
    latex = "\\begin{longtable}{rrrl" 
    latex = paste0(latex, strrep("c", length(thisScoreNames)), "}\n")
    latex = paste0(latex, "\\caption{Mean ± 2SE of aggregate NTG scores and metrics for the ", 
                   tab$propVarCase[1], " sampling scenario.} \\label{tab:", tab$propVarCase[1], "}\\\\\n")
    latex = paste0(latex, "\\toprule\n")
    latex = paste0(latex, "$\\phi$ & $p_{\\tiny \\mbox{rep}}$ & $n$ & Model & ", paste(thisScoreTitles, collapse = " & "), " \\\\\n")
    latex = paste0(latex, "\\midrule\n")
    latex = paste0(latex, "\\endfirsthead\n")
    latex = paste0(latex, "\\toprule\n")
    latex = paste0(latex, "$\\phi$ & $p_{\\tiny \\mbox{rep}}$ & $n$ & Model & ", paste(thisScoreTitles, collapse = " & "), " \\\\\n")
    latex = paste0(latex, "\\midrule\n")
    latex = paste0(latex, "\\endhead\n")
    
    
    # # Add rows
    # for (i in seq_len(nrow(summaryTab))) {
    #   row = summaryTab[i, ]
    #   # setting = paste0("phi=", row$phi, ", prep=", row$prep, ", n=", row$n)
    #   setting = paste0("$\\phi$ & $p_{\\tiny \\mbox{rep}}$ & $n$ & ", paste(thisScoreTitles, collapse = " & "), " \\\\\n")
    #   values = sapply(thisScoreNames, function(score) {
    #     meanVal = row[[paste0(score, "mean")]]
    #     seVal = row[[paste0(score, "se")]]
    #     if (grepl("Coverage", score)) {
    #       meanVal = round(meanVal * 100)
    #       seVal = round(qnorm(.975) * seVal * 100)
    #       sprintf("%d (%d)", meanVal, seVal)
    #     } else {
    #       sprintf("%.3f (%.3f)", meanVal, qnorm(.975) * seVal)
    #     }
    #   })
    #   latex = paste0(latex, setting, " & ", paste(values, collapse = " & "), " \\\\\n")
    # }
    
    
    # Initialize previous values
    prevPhi = prevPrep = prevN = NA
    
    for (i in seq_len(nrow(summaryTab))) {
      row = summaryTab[i, ]
      phiVal = row$phi
      prepVal = row$repelAreaProp
      nVal = row$n
      
      # Determine whether to print setting values
      phiPrint = if (!identical(phiVal, prevPhi)) phiVal else ""
      prepPrint = if (!identical(prepVal, prevPrep)) prepVal else ""
      nPrint = if (!identical(nVal, prevN)) nVal else ""
      
      # Update previous values
      prevPhi = phiVal
      prevPrep = prepVal
      prevN = nVal
      
      # Format score values
      values = sapply(thisScoreNames, function(score) {
        meanVal = row[[paste0(score, "mean")]]
        seVal = row[[paste0(score, "se")]]
        if (grepl("Coverage", score)) {
          meanVal = round(meanVal * 100)
          seVal = round(qnorm(.975) * seVal * 100)
          sprintf("%d (%d)", meanVal, seVal)
        } else {
          sprintf("%.3f (%.3f)", meanVal, qnorm(.975) * seVal)
        }
      })
      
      # Print row with conditional setting values
      latex = paste0(latex, phiPrint, " & ", prepPrint, " & ", nPrint, " & ", row$Model, " & ", paste(values, collapse = " & "), " \\\\\n")
    }
    
    
    latex = paste0(latex, "\\bottomrule\n\\end{longtable}\n")
    return(latex)
  }
  
  # precompute special cases: seismic model, uniform sampling
  tabSeismic = mergedTab %>% filter(fitModFunI == 0)
  tabUniform = mergedTab %>% filter(propVarCase == "uniform")
  
  
  # adjust mergedTab to not include seismic results (they will be added in later)
  mergedTab = mergedTab %>% filter(fitModFunI != 0)
  
  # Make tables: ----
  for (case in propVarCases) {
    tabBase = mergedTab %>% filter(propVarCase %in% c(case, "uniform"))
    
    # filter out SPDED and SPDEK models from seismic case
    if(case == "seismic") {
      tabBase = tabBase %>% filter(!(Model %in% c("SPDED", "SPDEK")))
    }
    
    cat(getLatexLongtable(tabBase))
    browser()
  }
  
  # Make plots: ----
  for (case in propVarCases) {
    tabBase = mergedTab %>% filter(propVarCase %in% c(case, "uniform"))
    
    # filter out SPDED and SPDEK models from seismic case
    if(case == "seismic") {
      tabBase = tabBase %>% filter(!(Model %in% c("SPDED", "SPDEK")))
    }
    
    # Boxplots and tables vs n: ----
    if(doBoxN) {
      print("boxplots vs n...")
      for (sampI in sort(unique(tabBase$sampleParI))) {
        tab = tabBase %>% filter(sampleParI == sampI)
        # tab = bind_rows(tab, tabSeismic)
        
        thisDirRoot = adaptScen
        adaptFRoot = ""  # for now. Will be changed later depending on the type of plots
        
        fileRoot = paste0("i", sampI, "_", case, "_repA", unique(tab$repelAreaProp), "_pref", unique(tab$phi), "_", adaptScen)
        figDir = file.path("figures/simStudy", thisDirRoot, fileRoot)
        
        dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
        
        for (scoreTypeI in 1:length(scoreTypesNamed)) {
          # scoreCols = grep(scoreType, names(tab), value = TRUE)
          # scoreCols = grep(paste0("_", scoreType), names(tab), value = TRUE)
          scoreType = scoreTypes[scoreTypeI]
          scoreTypeName = scoreTypesNamed[[scoreTypeI]]
          scoreTypeNameRoot = scoreTypesNameRoot[[scoreTypeI]]
          scoreCols = grep(scoreType, names(tab), value = TRUE)
          
          if(scoreTypeNameRoot == "par") {
            # only include plots for these parameters
            scoreCols = scoreCols[greplAny(c("gamma", "spatialRange", "spatialVar", "errorVar", "seismic_y", "seismic_p", "design"), scoreCols)]
          }
          
          for (scoreCol in scoreCols) {
            
            thisTab = tab[!is.na(tab[[scoreCol]]),]
            
            scoreColName = sub(paste0(scoreType, "$"), "", scoreCol)
            adaptFRoot = ""
            if(adaptScen != "batch") {
              adaptFRoot = paste0("_", adaptType)
            }
            fname = paste0("figures/simStudy/", thisDirRoot, "/", fileRoot, "/", scoreTypeNameRoot, "_", fileRoot, adaptFRoot, "_", scoreColName, ".pdf")
            
            makeBoxplot(thisTab, parName="n", scoreCol=scoreCol, fixedParNames=c("phi", "repelAreaProp"), fname=fname)
            
            
            # 
            # 
            # unique_models = unique(tab$Model)
            # 
            # 
            # # restrict to models actually present in thisTab, 
            # # but preserve the order from baseCols
            # presentModels = intersect(names(baseCols), as.character(thisTab$Model))
            # thisModCols = baseCols[presentModels]
            # 
            # # set factor levels of Model in the same order as colors
            # thisTab$Model = factor(thisTab$Model, levels = names(thisModCols))
            # 
            # # pdf(paste0("figures/simStudy/", fileRoot, "/", type, thisVar, "_", fileRoot, ".pdf"), width=5, height=5)
            # adaptFRoot = ""
            # if(adaptScen != "batch") {
            #   adaptFRoot = paste0("_", adaptType)
            # }
            # 
            # thisTab = thisTab[!is.na(thisTab[[scoreCol]]),]
            # thisTab[[scoreCol]] = as.numeric(thisTab[[scoreCol]])
            # levelsN = sort(unique(as.numeric(thisTab$n)))
            # thisTab$n = factor(thisTab$n, levels = as.character(levelsN))
            # 
            # # fname = paste0("figures/simStudy/", thisDirRoot, "/", fileSubRoot, "/", fileRoot, "/", type, "_", fileRoot, adaptFRoot, "_", scoreColName, ".pdf")
            # fname = paste0("figures/simStudy/", thisDirRoot, "/", fileRoot, "/", scoreTypeNameRoot, "_", fileRoot, adaptFRoot, "_", scoreColName, ".pdf")
            # pdf(fname, width=5, height=5)
            # 
            # # Create the boxplot
            # p = ggplot(thisTab, aes(x = factor(n), y = .data[[scoreCol]], fill = Model))
            # if(!(scoreColName %in% c("Coverage80", "Coverage95"))) {
            #   p = p + geom_boxplot()
            # }
            # p = p +
            #   # stat_summary(fun = mean, geom = "point", shape = 20, size = 1, color = "black", 
            #   #              position = position_dodge(width = 0.75)) +
            #   stat_summary(fun = mean, geom = "point", shape = 21, size = 2, 
            #                color = "black", aes(fill = Model),
            #                position = position_dodge(width = 0.75)) + 
            #   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black", 
            #                position = position_dodge(width = 0.75)) +
            #   # scale_fill_manual(values = setNames(thisModCols[1:length(unique_models)], unique_models)) +
            #   scale_fill_manual(values = thisModCols) + 
            #   labs(
            #     # title = paste0(myTitleCase(scoreColName), " vs. n (", typeName, ")"),
            #     title = paste0(myTitleCase(scoreColName), " vs. n"),
            #     x = "n",
            #     y = myTitleCase(scoreColName),
            #     fill = "Model"
            #   ) +
            #   theme_minimal()
            # if(!(scoreColName %in% c("Bias", "Var", "Width80", "Width95", "Coverage80", "Coverage95", 
            #                          "pref", "seismic_y", "seismic_p", "design"))) {
            #   # seismic estimates have 0 variance
            #   p = p + scale_y_log10()
            # }
            # if(grepl("Coverage", scoreColName)) {
            #   cvg = as.numeric(substr(scoreColName, nchar(scoreColName)-1, nchar(scoreColName))) / 100
            #   p = p + geom_hline(yintercept = cvg, color = "darkgrey", linetype = "dashed") 
            # }
            # 
            # print(p)
            # 
            # dev.off()
          }
        }
      }
    }
    
    
    # Plot vs phi ----
    if(doBoxPhi) {
      print("boxplots vs phi...")
      
      fixedParName = "repelAreaProp"
      parName = "phi"
      
      for (repelVal in sort(unique(tabBase$repelAreaProp))) {
        validNs = sort(unique(tabBase$n[tabBase$repelAreaProp == repelVal & repelVal * tabBase$n <= 0.3]))
        
        for (nVal in validNs) {
          tab = tabBase %>% filter(n == nVal & repelAreaProp == repelVal)
          # tab = bind_rows(tab, tabSeismic)
          
          thisFileRoot = paste0(case, "_prefParAll_repelAreaProp", repelVal, "_n", nVal, "_", adaptScen)
          
          figDir = paste0("figures/simStudy/", adaptScen, "/", case, "_prefParAll_repelAreaProp", repelVal, "_n", nVal, "_", adaptScen)
          dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
          
          for (scoreTypeI in 1:length(scoreTypesNamed)) {
            # scoreCols = grep(scoreType, names(tab), value = TRUE)
            # scoreCols = grep(paste0("_", scoreType), names(tab), value = TRUE)
            scoreType = scoreTypes[scoreTypeI]
            scoreTypeName = scoreTypesNamed[[scoreTypeI]]
            scoreTypeNameRoot = scoreTypesNameRoot[[scoreTypeI]]
            scoreCols = grep(scoreType, names(tab), value = TRUE)
            
            if(scoreTypeNameRoot == "par") {
              # only include plots for these parameters
              scoreCols = scoreCols[greplAny(c("gamma", "spatialRange", "spatialVar", "errorVar", "seismic_y", "seismic_p", "design"), scoreCols)]
            }
            
            for (scoreCol in scoreCols) {
              scoreColName = sub(paste0(scoreType, "$"), "", scoreCol)
              
              thisTab = tab[!is.na(tab[[scoreCol]]),]
              
              adaptFRoot = ""
              if(adaptScen != "batch") {
                adaptFRoot = paste0("_", adaptType)
              }
              # fname = paste0("figures/simStudy/", thisDirRoot, "/", fileRoot, "/", scoreTypeNameRoot, "_", fileRoot, adaptFRoot, "_", scoreColName, ".pdf")
              fname = paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot, "/", scoreTypeNameRoot, "_", thisFileRoot, "_", scoreColName, ".pdf")
              
              makeBoxplot(thisTab, parName=parName, scoreCol=scoreCol, fixedParNames=c("n", fixedParName), fname=fname)
              
              
              # # Ensure the number of colors matches the number of unique models
              # unique_models <- unique(tab$Model)
              # if (adaptScen == "batch") {
              #   baseCols = modCols
              # } else {
              #   if (adaptType == "spde") {
              #     baseCols = modCols
              #   } else if (adaptType == "self") {
              #     baseCols = modColsSelf
              #   } else {
              #     baseCols = modColsComb
              #   }
              # }
              # 
              # # restrict to models actually present in thisTab, 
              # # but preserve the order from baseCols
              # presentModels = intersect(names(baseCols), as.character(thisTab$Model))
              # thisModCols = baseCols[presentModels]
              # 
              # # set factor levels of Model in the same order as colors
              # thisTab$Model = factor(thisTab$Model, levels = names(thisModCols))
              # 
              # 
              # if(!dir.exists(paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot))) {
              #   dir.create(paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot))
              # }
              # 
              # thisTab = thisTab[!is.na(thisTab[[scoreCol]]),]
              # thisTab[[scoreCol]] = as.numeric(thisTab[[scoreCol]])
              # 
              # pdf(paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot, "/", scoreTypeNameRoot, "_", thisFileRoot, "_", scoreColName, ".pdf"), width=5, height=5)
              # 
              # # Create the plot
              # if(parName == "prefPar") {
              #   p = ggplot(thisTab, aes(x = factor(prefPar), y = .data[[scoreCol]], fill = Model))
              # } else {
              #   p = ggplot(thisTab, aes(x = factor(repelAreaProp), y = .data[[scoreCol]], fill = Model))
              # }
              # if(!(scoreColName %in% c("Coverage80", "Coverage95"))) {
              #   p = p + geom_boxplot()
              # }
              # p = p +
              #   # stat_summary(fun = mean, geom = "point", shape = 20, size = 1, color = "black", 
              #   #              position = position_dodge(width = 0.75)) +
              #   stat_summary(fun = mean, geom = "point", shape = 21, size = 2, 
              #                color = "black", aes(fill = Model),
              #                position = position_dodge(width = 0.75)) + 
              #   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black", 
              #                position = position_dodge(width = 0.75)) +
              #   # scale_fill_manual(values = setNames(modCols[1:length(unique_models)], unique_models)) +
              #   scale_fill_manual(values = thisModCols) + 
              #   labs(
              #     # title = paste0(myTitleCase(scoreColName), " vs. ", parName, 
              #     #                " (", adaptScenCap, ", ", typeName, ", ", 
              #     #                fixedParName, "=", fixedParVal, ", n=", thisN, ")"),
              #     title = paste0(myTitleCase(scoreColName), " vs. ", parName, 
              #                    " (", fixedParName, "=", fixedParVal, ", n=", thisN, ")"),
              #     x = parName,
              #     y = myTitleCase(scoreColName),
              #     fill = "Model"
              #   ) +
              #   theme_minimal()
              # if(!(scoreColName %in% c("Bias", "Var", "Width80", "Width95", "Coverage80", "Coverage95", 
              #                     "pref", "seismic_y", "seismic_p", "design"))) {
              #   # seismic estimates have 0 variance
              #   p = p + scale_y_log10()
              # }
              # if(grepl("Coverage", scoreColName)) {
              #   cvg <- as.numeric(substr(scoreColName, nchar(scoreColName)-1, nchar(scoreColName))) / 100
              #   p = p + geom_hline(yintercept = cvg, color = "darkgrey", linetype = "dashed") 
              # }
              # 
              # print(p)
              # 
              # dev.off()
            }
          }
        }
      }
    }
    
    
    # Plot vs repelAreaProp ----
    if(doBoxRep) {
      print("boxplots vs repelAreaProp...")
      
      fixedParName = "phi"
      parName = "repelAreaProp"
      
      for (nVal in sort(unique(tabBase$n))) {
        for (prefVal in sort(unique(tabBase$phi))) {
          tab = tabBase %>% filter(n == nVal & phi == prefVal)
          # tab = bind_rows(tab, tabSeismic)
          
          if(nrow(tab) == 0) {
            print(paste0("no values to plot for n=", nVal, ", phi=", prefVal))
            next
          }
          
          
          figDir = paste0("figures/simStudy/", adaptScen, "/", case, "_repelAreaPropAll_prefPar", prefVal, "_n", nVal, "_", adaptScen)
          dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
          
          thisFileRoot = paste0(case, "_repelAreaPropAll_prefPar", prefVal, "_n", nVal, "_", adaptScen)
          
          for (scoreTypeI in 1:length(scoreTypesNamed)) {
            # scoreCols = grep(scoreType, names(tab), value = TRUE)
            # scoreCols = grep(paste0("_", scoreType), names(tab), value = TRUE)
            scoreType = scoreTypes[scoreTypeI]
            scoreTypeName = scoreTypesNamed[[scoreTypeI]]
            scoreTypeNameRoot = scoreTypesNameRoot[[scoreTypeI]]
            scoreCols = grep(scoreType, names(tab), value = TRUE)
            
            if(scoreTypeNameRoot == "par") {
              # only include plots for these parameters
              scoreCols = scoreCols[greplAny(c("gamma", "spatialRange", "spatialVar", "errorVar", "seismic_y", "seismic_p", "design"), scoreCols)]
            }
            
            for (scoreCol in scoreCols) {
              scoreColName = sub(paste0(scoreType, "$"), "", scoreCol)
              
              thisTab = tab[!is.na(tab[[scoreCol]]),]
              
              adaptFRoot = ""
              if(adaptScen != "batch") {
                adaptFRoot = paste0("_", adaptType)
              }
              
              fname = paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot, "/", scoreTypeNameRoot, "_", thisFileRoot, "_", scoreColName, ".pdf")
              
              makeBoxplot(thisTab, parName=parName, scoreCol=scoreCol, fixedParNames=c("n", fixedParName), fname=fname)
              
              # # Ensure the number of colors matches the number of unique models
              # unique_models <- unique(tab$Model)
              # if (adaptScen == "batch") {
              #   baseCols = modCols
              # } else {
              #   if (adaptType == "spde") {
              #     baseCols = modCols
              #   } else if (adaptType == "self") {
              #     baseCols = modColsSelf
              #   } else {
              #     baseCols = modColsComb
              #   }
              # }
              # 
              # # restrict to models actually present in thisTab, 
              # # but preserve the order from baseCols
              # presentModels = intersect(names(baseCols), as.character(thisTab$Model))
              # thisModCols = baseCols[presentModels]
              # 
              # # set factor levels of Model in the same order as colors
              # thisTab$Model = factor(thisTab$Model, levels = names(thisModCols))
              # 
              # 
              # if(!dir.exists(paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot))) {
              #   dir.create(paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot))
              # }
              # 
              # thisTab = thisTab[!is.na(thisTab[[scoreCol]]),]
              # thisTab[[scoreCol]] = as.numeric(thisTab[[scoreCol]])
              # 
              # pdf(paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot, "/", scoreTypeNameRoot, "_", thisFileRoot, "_", scoreColName, ".pdf"), width=5, height=5)
              # 
              # # Create the plot
              # if(parName == "prefPar") {
              #   p = ggplot(thisTab, aes(x = factor(prefPar), y = .data[[scoreCol]], fill = Model))
              # } else {
              #   p = ggplot(thisTab, aes(x = factor(repelAreaProp), y = .data[[scoreCol]], fill = Model))
              # }
              # if(!(scoreColName %in% c("Coverage80", "Coverage95"))) {
              #   p = p + geom_boxplot()
              # }
              # p = p +
              #   # stat_summary(fun = mean, geom = "point", shape = 20, size = 1, color = "black", 
              #   #              position = position_dodge(width = 0.75)) +
              #   stat_summary(fun = mean, geom = "point", shape = 21, size = 2, 
              #                color = "black", aes(fill = Model),
              #                position = position_dodge(width = 0.75)) + 
              #   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black", 
              #                position = position_dodge(width = 0.75)) +
              #   # scale_fill_manual(values = setNames(modCols[1:length(unique_models)], unique_models)) +
              #   scale_fill_manual(values = thisModCols) + 
              #   labs(
              #     # title = paste0(myTitleCase(scoreColName), " vs. ", parName, 
              #     #                " (", adaptScenCap, ", ", typeName, ", ", 
              #     #                fixedParName, "=", fixedParVal, ", n=", thisN, ")"),
              #     title = paste0(myTitleCase(scoreColName), " vs. ", parName, 
              #                    " (", fixedParName, "=", fixedParVal, ", n=", thisN, ")"),
              #     x = parName,
              #     y = myTitleCase(scoreColName),
              #     fill = "Model"
              #   ) +
              #   theme_minimal()
              # if(!(scoreColName %in% c("Bias", "Var", "Width80", "Width95", "Coverage80", "Coverage95", 
              #                          "pref", "seismic_y", "seismic_p", "design"))) {
              #   # seismic estimates have 0 variance
              #   p = p + scale_y_log10()
              # }
              # if(grepl("Coverage", scoreColName)) {
              #   cvg <- as.numeric(substr(scoreColName, nchar(scoreColName)-1, nchar(scoreColName))) / 100
              #   p = p + geom_hline(yintercept = cvg, color = "darkgrey", linetype = "dashed") 
              # }
              # 
              # print(p)
              # 
              # dev.off()
            }
          }
        }
      }
    }
    
    if(doScatterGamma) {
      print("scatterplots vs gamma...")
      
      # plot prefEst vs Bias ----
      for (repelVal in sort(unique(tabBase$repelAreaProp))) {
        validNs = sort(unique(tabBase$n[tabBase$repelAreaProp == repelVal & repelVal * tabBase$n <= 0.3]))
        
        for (nVal in validNs) {
          
          for (prefVal in sort(unique(tabBase$phi))) {
            
            tab = tabBase %>% filter(n == nVal & repelAreaProp == repelVal & phi == prefVal)
            
            if(nrow(tab) == 0) {
              print(paste0("no values to plot for n=", nVal, ", phi=", prefVal, ", prep=", repelVal))
              next
            }
            
            figDir = paste0("figures/simStudy/", adaptScen, "/", case, "_prefPairs", "_prefPar", prefVal, "_repelAreaProp", repelVal, "_n", nVal, "_", adaptScen)
            dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
            
            for (scoreTypeI in 1:length(scoreTypesNamed)) {
              # scoreCols = grep(scoreType, names(tab), value = TRUE)
              # scoreCols = grep(paste0("_", scoreType), names(tab), value = TRUE)
              scoreType = scoreTypes[scoreTypeI]
              scoreTypeName = scoreTypesNamed[[scoreTypeI]]
              scoreTypeNameRoot = scoreTypesNameRoot[[scoreTypeI]]
              scoreCols = grep(scoreType, names(tab), value = TRUE)
              
              if(scoreTypeNameRoot == "par") {
                # skip this case
                next
              }
              
              for (scoreCol in scoreCols) {
                scoreColName = sub(paste0(scoreType, "$"), "", scoreCol)
                
                thisTab = tab[!is.na(tab[[scoreCol]]),]
                
                adaptFRoot = ""
                if(adaptScen != "batch") {
                  adaptFRoot = paste0("_", adaptType)
                }
                
                # prefEst
                thisFileRoot = paste0(case, "_prefPairs", "_prefPar", prefVal, "_repelAreaProp", repelVal, "_n", nVal, "_", adaptScen)
                thisFileRoot2 = paste0(case, "_prefPairs", "_prefPar", prefVal, "_repelAreaProp", repelVal, "_n", nVal, "_", adaptScen)
                fname = paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot, "/", scoreTypeNameRoot, "_", thisFileRoot2, "_", scoreColName, ".pdf")
                
                makeScatterplot(thisTab, parName="gamma_param", scoreCol=scoreCol, fixedParNames=c("n", "phi", "repelAreaProp"), 
                                fname=fname)
                
                # corWells
                thisFileRoot2 = paste0(case, "_corWellsPairs", "_prefPar", prefVal, "_repelAreaProp", repelVal, "_n", nVal, "_", adaptScen)
                fname = paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot, "/", scoreTypeNameRoot, "_", thisFileRoot2, "_", scoreColName, ".pdf")
                
                # make the scatterplot
                makeScatterplot(thisTab, parName="corSeisTruthWells", scoreCol=scoreCol, fixedParNames=c("n", "phi", "repelAreaProp"), 
                                fname=fname)
                
                # corTruth
                thisFileRoot2 = paste0(case, "_corTruePairs", "_prefPar", prefVal, "_repelAreaProp", repelVal, "_n", nVal, "_", adaptScen)
                fname = paste0("figures/simStudy/", thisDirRoot, "/", thisFileRoot, "/", scoreTypeNameRoot, "_", thisFileRoot2, "_", scoreColName, ".pdf")
                
                # make the scatterplot
                makeScatterplot(thisTab, parName="corSeisTruthTrue", scoreCol=scoreCol, fixedParNames=c("n", "phi", "repelAreaProp"), 
                                fname=fname)
                
                
                
              }
            }
          }
        }
      }
    }
    
  }
  
  message("All plots saved.")
}



# call on the cluster to copy figures we need into the manuscript directory ~/fig
refreshManuscriptFigDir = function() {
  
  copyDirFiltered(srcDir=paste0("figures/simStudy/"), 
                  dstDir=paste0("~/fig/simStudy/"), 
                  includeSubstr = c("agg", "par", "mean.*Coverage95"), excludeSubstr = c("80", "_Var", "RMSE", "Pairs", "min", "max", "worst"))
  
  copyDirFiltered(srcDir=paste0("figures/simStudy/"), 
                  dstDir=paste0("~/fig/simStudy/"), 
                  includeSubstr = c("agg_cluster_prefPairs.*Bias"))
  
}





