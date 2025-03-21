# this script fits SPDE model to the data and generates predictions

# generate default priors for SPDE model
# from Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
# NOTE: by default, this constructs a spde prior with unit median marginal variance 
#       and median effective range equal to a fifth of the spatial range 
# 
# Inputs:
# mesh: SPDE mesh. Could be created by e.g. getSPDEmesh()
# U, alpha: the threshold and probability of crossing the threshold for the 
#              standard deviation in the PC prior. By default, p(marginalSD > 1) = 0.05
# medianRange: by default, one fifth of the domain diameter
# 
# Outputs: prior to give to INLA
getSPDEprior = function(mesh, U=1, alpha=0.05, medianRange=NULL) {
  size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
  
  # set default median range
  if(is.null(medianRange))
    range0 <- size/5
  else
    range0 = medianRange
  
  # set PC prior on SPDE model based on the mesh and hyperparameters
  spde = inla.spde2.pcmatern(mesh, prior.range=c(range0, 0.5), prior.sigma = c(U, alpha))
  
  # return results
  spde
}

# get a reasonable default mesh triangulation for the SPDE model for a general 
# set of locations
# 
# Inputs: 
# locs: general set of locations in 2d
# n: target number of triangles in mesh. See ?inla.mesh.2d
# max.n: maximum number of triangle in mesh. See ?inla.mesh.2d
# max.edge: 2-vector of max mesh edge lengths for SPDE model. See ?inla.mesh.2d
# offset: relates to extension of domain for mesh creation. See ?inla.mesh.2d
# cutoff: minimum distance between points. See ?inla.mesh.2d
# doPlot: whether to plot the resulting mesh for testing purposes
# ...: other arguments to inla.mesh.2d
# 
# Outputs:
# mesh to give to INLA
getSPDEmesh = function(locs=cbind(c(-1, -1, 1, 1), c(-1, 1, -1, 1)), n=3500, 
                       max.n=5000, doPlot=TRUE, max.edge=c(.01, .1), 
                       offset=-.08, cutoff=.005, ...) {
  
  
  # generate mesh on R2
  mesh = inla.mesh.2d(loc.domain=locs, n=n, max.n=max.n, offset=offset, 
                      cutoff=cutoff, max.edge=max.edge, ...)
  
  # plot the mesh if user wants
  if(doPlot) {
    plot(mesh)
  }
  
  mesh
}

# get a reasonable default mesh triangulation for the SPDE model for a 
# rectangular spatial domain
# 
# inputs: 
# lowerLeft: a 2-vector of form c(x, y) of lower left corner of the rectangular domain
# width: width of domain
# height: height of the domain
# n: target number of triangles in mesh. See ?inla.mesh.2d
# max.n: maximum number of triangle in mesh. See ?inla.mesh.2d
# scale: scaling parameter for setting defaults
# max.edge: 2-vector of max mesh edge lengths for SPDE model. See ?inla.mesh.2d
# offset: relates to extension of domain for mesh creation. See ?inla.mesh.2d
# cutoff: minimum distance between points. See ?inla.mesh.2d
# doPlot: whether to plot the resulting mesh for testing purposes
# ...: other arguments to inla.mesh.2d
# 
# Outputs:
# mesh to give to INLA
getSPDEmeshRect = function(lowerLeft=c(0,0), width=1, height=1, n=3500, max.n=5000, 
                           scale=max(c(width, height)), max.edge=c(.01, .1)*scale, 
                           offset=-.015, cutoff=.005*scale, doPlot=TRUE, ...) {
  
  locs = cbind(c(0, width, width, 0), c(0, 0, height, height))
  locs = sweep(locs, 2, lowerLeft, "+")
  
  # generate mesh on R2
  mesh = inla.mesh.2d(loc.domain=locs, n=n, max.n=max.n, offset=offset, 
                      cutoff=cutoff, max.edge=max.edge, ...)
  
  # plot the mesh if user wants
  if(doPlot) {
    plot(mesh)
    
    # add domain boundary
    domLocs = rbind(locs, locs[1,])
    lines(domLocs[,1], domLocs[,2], col="blue")
  }
  
  mesh
}

# get a reasonable default mesh triangulation for the SPDE model for the sim study data
# 
# Inputs:
# doPlot: whether or not to plot the mesh after generating it
# 
# Outputs:
# mesh to give to INLA
getSPDEmeshSimStudy = function(doPlot=FALSE) {
  
  getSPDEmeshRect(doPlot=doPlot)
}

# Fits SPDE model to well and seismic data from simulation study
# 
# Inputs:
# wellDat: A data.frame with columns:
#   east: easting
#   north: northing
#   volFrac: sand volume fraction (nett og gross)
# seismicDat: a data.frame containing a grid of of seismic estimates with columns:
#   east: easting
#   north: northing
#   est: central estimate
# predGrid: grid over which to make predictions
# transform: how to transform obsValues prior to modeling
# invTransform: inverse of transform. Used to backtransform predictions prior to aggregation
# mesh: SPDE mesh
# prior: SPDE prior
# significanceCI: the credible level of the CIs (e.g., .8 means 80% CI)
# int.strategy: inla integration strategy
# strategy: inla strategy
# nPostSamples: number posterior draws
# verbose: verbose argument to inla
# seed: random seed. Not set if NULL
# family: currently only normal is supported
# doModAssess: whether or not to calculate CPO, DIC, and WAIC
# previousFit: a previous INLA model fit used to initialize optimization
# improperCovariatePrior: if TRUE, N(0, infty) prior on covariates (aside from 
#                         intercept, which already has this prior)
# fixedParameters: A list of parameters to fix in the model rather than infer. 
#                  Contains some of all of the elements: spde$effRange, 
#                  spde$margVar, familyPrec, clusterPrec, beta (NOT TESTED)
# experimentalMode: Whether to use INLA variational inference tools (NOT TESTED)
# 
# Outputs:
# INLA model, predictions, summary statistics, input data, posterior draws, etc.
fitSPDEsimDat = function(wellDat, seismicDat, 
                         predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
                         transform=logit, invTransform=expit, 
                         mesh=getSPDEmeshSimStudy(), prior=getSPDEprior(mesh), 
                         significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                         nPostSamples=1000, verbose=TRUE, seed=123, 
                         family="normal", doModAssess=FALSE, previousFit=NULL, 
                         improperCovariatePrior=TRUE, fixedParameters=NULL, 
                         experimentalMode=FALSE) {
  
  # set defaults
  # family = match.arg(family)
  
  # construct prediction covariates
  xPred = cbind(1, seismicDat$seismicEst)
  
  # interpolate seismic data to the well points
  wellSeismicEsts = bilinearInterp(bilinearInterp[,1:2], seismicDat)
  
  # construct well data covariates
  xObs = cbind(1, wellSeismicEsts)
  
  # set observations
  obsValues = wellDat$volFrac
  obsCoords = cbind(wellDat$east, wellDat$north)
  
  fitSPDE(obsCoords, obsValues, xObs, 
          predPts, xPred, 
          mesh, prior, 
          significanceCI, int.strategy, strategy, 
          nPostSamples, verbose, link, seed, 
          family, doModAssess, previousFit, 
          improperCovariatePrior=improperCovariatePrior, 
          fixedParameters=fixedParameters, experimentalMode=experimentalMode)
}


# function for fitting the SPDE model to data
# 
# Inputs:
# obsCoords: data.frame with easting and northing columns for observations
# xObs: matrix with intercept and covariate information for observations
# predCoords: data.frame with easting and northing columns for prediction grid
# xPred: matrix with intercept and covariate information for prediction grid
# transform: how to transform obsValues prior to modeling
# invTransform: inverse of transform. Used to backtransform predictions prior to aggregation
# mesh: SPDE mesh
# prior: SPDE prior
# significanceCI: the credible level of the CIs (e.g., .8 means 80% CI)
# int.strategy: inla integration strategy
# strategy: inla strategy
# nPostSamples: number posterior draws
# verbose: verbose argument to inla
# link: link=1 is the canonical link in inla
# seed: random seed. Not set if NULL
# family: currently only normal is supported
# doModAssess: whether or not to calculate CPO, DIC, and WAIC
# previousFit: a previous INLA model fit used to initialize optimization
# improperCovariatePrior: if TRUE, N(0, infty) prior on covariates (aside from 
#                         intercept, which already has this prior)
# fixedParameters: A list of parameters to fix in the model rather than infer. 
#                  Contains some of all of the elements: spde$effRange, 
#                  spde$margVar, familyPrec, clusterPrec, beta (NOT TESTED)
# experimentalMode: Whether to use INLA variational inference tools (NOT TESTED)
# 
# Outputs:
# INLA model, predictions, summary statistics, input data, posterior draws, etc.
fitSPDE = function(obsCoords, obsValues, xObs=matrix(rep(1, length(obsValues)), ncol=1), 
                   predCoords, xPred = matrix(rep(1, nrow(predCoords)), ncol=1), 
                   transform=I, invTransform=I, 
                   mesh=getSPDEmesh(obsCoords), prior=getSPDEprior(mesh), 
                   significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                   nPostSamples=1000, verbose=TRUE, link=1, seed=NULL, 
                   family=c("normal", "binomial", "betabinomial"), 
                   doModAssess=FALSE, 
                   previousFit=NULL, improperCovariatePrior=TRUE, 
                   fixedParameters=NULL, experimentalMode=FALSE) {
  family = match.arg(family)
  startTime = proc.time()[3]
  if(!is.null(seed))
    set.seed(seed)
  
  if(experimentalMode) {
    if(strategy != "gaussian") {
      stop("only gaussian integration is possible if experimentalMode is TRUE")
    }
  }
  
  startTimeDefineModel = proc.time()[3]
  
  # fix SPDE model parameters if necessary
  if(!is.null(fixedParameters)) {
    prior = getSPDEModelFixedPar(mesh, effRange=fixedParameters$spde$effRange, 
                                   margVar=fixedParameters$spde$margVar)
  }
  
  # set family prior
  control.family = list(hyper = list(prec = list(prior="loggamma", param=c(0.1,0.1))))
  
  if(!is.null(fixedParameters$familyPrec)) {
    # fix the family precision parameter on INLA's latent scale
    control.family = list(initial=log(fixedParameters$familyPrec), fixed=TRUE)
  }
  
  # construct A matrix for observations
  m = nrow(obsCoords)
  AEst = inla.spde.make.A(mesh, loc = obsCoords)
  
  # construct A matrix for predictions
  APred = inla.spde.make.A(mesh, loc = predCoords)
  
  # make inla stack
  ys = transform(obsValues)
  n = ncol(AEst) # number of basis elements
  nObs = length(ys) # number of observations
  nPreds = nrow(predCoords)
  latticeInds = 1:n
  cluster = 1:length(ys)
  
  if(!is.null(fixedParameters$beta)) {
    offsetEst = xObs %*% fixedParameters$beta
    xObs = NULL
    offsetPred = xPred %*% fixedParameters$beta
    xPred = NULL
  } else {
    offsetEst = NULL
    offsetPred = NULL
  }
  
  # construct the observation stack 
  if(family == "normal") {
    if(!is.null(xObs)) {
      stack.est = inla.stack(A =list(AEst, 1),
                             effects =list(field=latticeInds, X=xObs),
                             data =list(y=ys, link=1),
                             tag ="est",
                             remove.unused=FALSE)
    } else {
      stack.est = inla.stack(A =list(AEst),
                             effects =list(field=latticeInds),
                             data =list(y=ys, link=1),
                             tag ="est",
                             remove.unused=FALSE)
    }
  } else {
    stop(paste0("family ", family, " not supported"))
  }
  
  endTimeDefineModel = proc.time()[3]
  totalTimeDefineModel = endTimeDefineModel - startTimeDefineModel
  
  # make mesh index
  mesh.index <- inla.spde.make.index(name = "field", n.spde = prior$n.spde)
  
  # fit model
  control.inla = list(strategy=strategy, int.strategy=int.strategy)
  modeControl = inla.set.control.mode.default()
  if(!is.null(previousFit)) {
    # initialize the fitting process based on a previous optimum
    
    # modeControl$result = previousFit
    modeControl$theta = previousFit$mode$theta
    modeControl$x = previousFit$mode$x
    modeControl$restart = TRUE
  }
  
  # Set distributional quantiles we're interested in: median and based on significanceCI
  allQuantiles = c(0.5, (1-significanceCI) / 2, 1 - (1-significanceCI) / 2)
  
  # fixed effect priors: are they improper or not?
  if(improperCovariatePrior) {
    controlFixed=list(quantiles=allQuantiles, mean=0, prec=0)
  } else {
    controlFixed=list(quantiles=allQuantiles)
  }
  
  # construct the stack
  stack.full = stack.est
  stackDat = inla.stack.data(stack.full, spde=prior)
  
  # see: inla.doc("loggamma")
  # shape=.1, scale=10 for unit mean, variance 100 prior
  
  # setup the formula
  thisFormula = "y ~ -1 + f(field, model=prior)"
  if(!is.null(xObs)) {
    thisFormula = paste0(thisFormula, " + X")
  }
  
  thisFormula = as.formula(thisFormula)
  
  startModelFitTime = proc.time()[3]
  if(family == "normal") {
    mod = inla(#y ~ - 1 + X + f(field, model=prior), 
      thisFormula, 
      data = stackDat, 
      control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles), 
      family=family, verbose=verbose, control.inla=control.inla, 
      control.compute=list(config=TRUE, cpo=doModAssess, dic=doModAssess, waic=doModAssess), 
      control.mode=modeControl, 
      offset=offsetEst, 
      control.fixed=controlFixed, 
      control.family=control.family)
  } else {
    stop(paste0("Unsupported family: ", family))
  }
  endModelFitTime = proc.time()[3]
  totalModelFitTime = endModelFitTime - startModelFitTime
  
  # get predictive surface, SD, and data
  n = nrow(obsCoords)
  # obsInds = 1:n
  obsInds = inla.stack.index(stack.full, "est")$data
  predInds = inla.stack.index(stack.full, "pred")$data
  index = inla.stack.index(stack.full, "pred")$data
  linpreds = mod[["summary.fitted.values"]]$mean
  linpred.sd = mod[["summary.fitted.values"]]$sd
  
  preds = linpreds
  predSDs = linpred.sd
  
  # generate samples from posterior
  startTimePosteriorSampling = proc.time()[3]
  postSamples = inla.posterior.sample(nPostSamples, mod)
  endTimePosteriorSampling = proc.time()[3]
  totalTimePosteriorSampling = endTimePosteriorSampling - startTimePosteriorSampling
  latentMat = sapply(postSamples, function(x) {x$latent})
  
  
  if(family == "normal") {
    browser() # get rid of hard coded [1]:
    nuggetVars = sapply(postSamples, function(x) {1 / x$hyperpar[1]})
  }
  
  latentVarNames = rownames(postSamples[[1]]$latent)
  fieldIndices = which(grepl("field", latentVarNames))
  fixedIndices = which(grepl("X", latentVarNames))
  
  # generate predictions: first without nugget then add the nugget in
  if(length(xPred) != 0)
    fixedPredMat = xPred  %*% latentMat[fixedIndices,]
  else
    fixedPredMat = 0
  
  spatialPredMat = APred %*% latentMat[fieldIndices,]
  
  predMat = fixedPredMat + spatialPredMat
  
  if(!is.null(offsetPred)) {
    predMat = sweep(predMat, 1, offsetPred, "+")
  }
  
  # do the same for the observations
  if(length(xObs) != 0)
    fixedObsMat = xObs  %*% latentMat[fixedIndices,]
  else
    fixedObsMat = 0
  
  spatialObsMat = AEst %*% latentMat[fieldIndices,]
  
  obsMat = fixedObsMat + spatialObsMat
  
  if(!is.null(offsetEst)) {
    obsMat = sweep(obsMat, 1, offsetEst, "+")
  }
  
  # add in nugget if necessary
  if(family == "normal") {
    # get cluster effect variance
    clusterVars = nuggetVars
    
    predMatNugget = predMat + matrix(rnorm(length(predMat), sd=rep(sqrt(clusterVars), each=nrow(predMat))), nrow=nrow(predMat))
    obsMatNugget = obsMat + matrix(rnorm(length(obsMat), sd=rep(sqrt(clusterVars), each=nrow(obsMat))), nrow=nrow(obsMat))
  } else {
    stop("family not supported")
  }
  
  # back transform predictions
  obsMat = invTransform(obsMat)
  obsMatNugget = invTransform(obsMatNugget)
  predMat = invTransform(predMat)
  predMatNugget = invTransform(predMatNugget)
  
  # get summary statistics
  obsEst = rowMeans(obsMat)
  obsSDs = apply(obsMat, 1, sd)
  obsLower = apply(obsMat, 1, quantile, probs=(1-significanceCI)/2)
  obsMedian = apply(obsMat, 1, median)
  obsUpper = apply(obsMat, 1, quantile, probs=1-(1-significanceCI)/2)
  
  obsNuggetEst = obsEst
  obsNuggetSDs = apply(obsMatNugget, 1, sd)
  obsNuggetLower = apply(obsMatNugget, 1, quantile, probs=(1-significanceCI)/2)
  obsNuggetMedian = apply(obsMatNugget, 1, median)
  obsNuggetUpper = apply(obsMatNugget, 1, quantile, probs=1-(1-significanceCI)/2)
  
  predEst = rowMeans(predMat)
  predSDs = apply(predMat, 1, sd)
  predLower = apply(predMat, 1, quantile, probs=(1-significanceCI)/2)
  predMedian = apply(predMat, 1, median)
  predUpper = apply(predMat, 1, quantile, probs=1-(1-significanceCI)/2)
  
  predNuggetEst = predEst
  predNuggetSDs = apply(predMatNugget, 1, sd)
  predNuggetLower = apply(predMatNugget, 1, quantile, probs=(1-significanceCI)/2)
  predNuggetMedian = apply(predMatNugget, 1, median)
  predNuggetUpper = apply(predMatNugget, 1, quantile, probs=1-(1-significanceCI)/2)
  
  # get aggregate summary statistics over prediction domain
  predAggEst = mean(predMat)
  predAggSDs = apply(colMeans(predMat), 1, sd)
  predAggLower = quantile(colMeans(predMat), probs=(1-significanceCI)/2)
  predAggMedian = quantile(colMeans(predMat), probs=.5)
  predAggUpper = quantile(colMeans(predMat), probs=1-(1-significanceCI)/2)
  
  if(length(xPred) != 0) {
    interceptSummary=mod$summary.fixed[1,1:5]
    fixedEffectSummary=mod$summary.fixed[,1:5]
  } 
  else {
    interceptSummary = matrix(rep(0, 5), nrow=1)
    fixedEffectSummary = mod$summary.fixed
  }
  rangeSummary=mod$summary.hyperpar[2,1:5]
  spatialSDSummary = mod$summary.hyperpar[3,1:5]
  
  # get posterior hyperparameter samples and transform them as necessary
  hyperMat = sapply(postSamples, function(x) {x$hyperpar})
  if(family == "normal") {
    clusterVarI = 1
    spatialRangeI = 2
    spatialSDI = 3
    if(!is.matrix(hyperMat)) {
      mat = NULL
    } else {
      mat = apply(hyperMat, 2, function(x) {c(totalVar=x[spatialSDI]^2+1/x[clusterVarI], spatialVar=x[spatialSDI]^2, errorVar=1/x[clusterVarI], 
                                              totalSD=sqrt(x[spatialSDI]^2+1/x[clusterVarI]), spatialSD=x[spatialSDI], errorSD=sqrt(1/x[clusterVarI]), 
                                              spatialRange=x[spatialRangeI])})
    }
  } else {
    stop("family not supported")
  }
  
  if(family == "normal")
    hyperNames = c("totalVar", "spatialVar", "errorVar", "totalSD", "spatialSD", "errorSD", "spatialRange")
  else 
    stop("family not supported")
  if(is.matrix(hyperMat)) {
    rownames(mat) = hyperNames
    
    getSummaryStatistics = function(draws) {
      c(Est=mean(draws, na.rm=TRUE), SD=sd(draws, na.rm=TRUE), 
        Qlower=quantile(probs=(1 - significanceCI) / 2, draws, na.rm=TRUE), 
        Q50=quantile(probs=0.5, draws, na.rm=TRUE), 
        Qupper=quantile(probs=1 - (1 - significanceCI) / 2, draws, na.rm=TRUE))
    }
    summaryNames = c("Est", "SD", "Qlower", "Q50", "Qupper")
    parameterSummaryTable = t(apply(mat, 1, getSummaryStatistics))
    colnames(parameterSummaryTable) = summaryNames
    
    # separate out default parameter summaries
    if(family == "normal") {
      sdSummary=parameterSummaryTable[6,]
      varSummary=parameterSummaryTable[3,]
      rangeSummary=parameterSummaryTable[7,]
    } else {
      stop("family not supported")
    }
  } else {
    parameterSummaryTable = NULL
    sdSummary = NULL
    varSummary = NULL
    rangeSummary = NULL
    overdispersionSummary = NULL
  }
  
  endTime = proc.time()[3]
  totalTime = endTime - startTime
  timings = data.frame(totalTime=totalTime, 
                       modelDefineTime=totalTimeDefineModel, 
                       modelFitTime=totalModelFitTime, 
                       posteriorSamplingTime=totalTimePosteriorSampling, 
                       otherTime=totalTime-(totalTimeDefineModel + totalModelFitTime + totalTimePosteriorSampling))
  timings$modelDefinePct = timings$modelDefineTime / timings$totalTime
  timings$modelFitTimePct = timings$modelFitTime / timings$totalTime
  timings$posteriorSamplingTimePct = timings$posteriorSamplingTime / timings$totalTime
  timings$otherTimePct = timings$otherTime / timings$totalTime
  
  list(mod=mod, 
       obsCoords=obsCoords, xObs=xObs, obsValues=obsValues, predCoords=predCoords, xPred=xPred, 
       obsEst=obsEst, obsSDs=obsSDs, obsLower=obsLower, obsMedian=obsMedian, obsUpper=obsUpper, 
       predEst=predEst, predSDs=predSDs, predLower=predLower, predMedian=predMedian, predUpper=predUpper, 
       predAggEst=predAggEst, predAggSDs=predAggSDs, predAggLower=predAggLower, predAggMedian=predAggMedian, predAggUpper=predAggUpper, 
       mesh=mesh, prior=prior, stack=stack.full, 
       interceptSummary=interceptSummary, fixedEffectSummary=fixedEffectSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, 
       parameterSummaryTable=parameterSummaryTable, 
       fixedEffectDraws=latentMat[fixedIndices,], 
       spatialPredMat=spatialPredMat, fixedPredMat=fixedPredMat, 
       spatialObsMat=spatialObsMat, fixedObsMat=fixedObsMat, 
       obsMat=obsMat, obsMatNugget=obsMatNugget, predMat=predMat, predMatNugget=predMatNugget, 
       hyperMat=hyperMat, timings=timings, sigmaEpsilonDraws=sqrt(clusterVars))
}







