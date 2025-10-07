# this script fits SPDE model to the data and generates predictions

# generate default priors for SPDE model
# from Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
# NOTE: by default, this constructs a spde prior with unit median marginal variance 
#       and median effective range equal to a fifth of the spatial range 
# 
# Inputs:
# mesh: SPDE mesh. Could be created by e.g. getSPDEmesh()
# U, alpha: the threshold and probability of crossing the threshold for the 
#              standard deviation in the PC prior. By default, p(marginalSD > 0.1) = 0.5
# medianRange: by default, one fifth of the domain diameter
# 
# Outputs: prior to give to INLA
getSPDEprior = function(mesh, U=0.1, alpha=0.5, medianRange=NULL) {
  size <- max(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
  
  # set default median range
  if(is.null(medianRange))
    range0 <- size/5 # 1483.947 for sim study mesh
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
  
  locs = cbind(c(0, width, width, 0, 0), c(0, 0, height, height, 0))
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
getSPDEmeshSimStudy = function(n=1500, max.n=2500, doPlot=FALSE, anisFac=1) {
  
  simStudyXlims = c(-12.5, 15012.5)/anisFac
  simStudyYlims = c(-8.3335, 5008.4336)
  
  getSPDEmeshRect(width=diff(simStudyXlims), 
                  height=diff(simStudyYlims), 
                  n=n, max.n=max.n, 
                  doPlot=doPlot)
}

getSPDEmeshTestDat = function(n=1500, max.n=2500, doPlot=FALSE) {
  
  simStudyXlims = c(0, 10000)
  simStudyYlims = c(0, 10000)
  
  getSPDEmeshRect(width=diff(simStudyXlims), 
                  height=diff(simStudyYlims), 
                  n=n, max.n=max.n, 
                  doPlot=doPlot)
}

getDomainSimStudy = function(returnClass=c("sf", "matrix")) {
  returnClass = match.arg(returnClass)
  
  simStudyXlims = c(-12.5, 15012.5)
  simStudyYlims = c(-8.3335, 5008.4336)
  
  locs = cbind(c(simStudyXlims[1], simStudyXlims[2], simStudyXlims[2], simStudyXlims[1], simStudyXlims[1]), 
               c(simStudyYlims[1], simStudyYlims[1], simStudyYlims[2], simStudyYlims[2], simStudyYlims[1]))
  
  if(returnClass == "matrix") {
    locs
  } else {
    require(sf)
    pol = st_polygon(
      list(
        locs
      )
    )
    pol
  }
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
# addKDE: Whether or not to add log of a kde density estimator as a covariate
# esthFromSeismic: whether to estimate h, kde's bandwidth parameter, from 
#                  seismic/grid data
# kde.args: arguments passed to kde function
# significanceCI: the credible level of the CIs (e.g., .8 means 80% CI)
# int.strategy: inla integration strategy
# strategy: inla strategy
# nPostSamples: number posterior draws
# verbose: verbose argument to inla
# seed: random seed. Not set if NULL
# family: currently only normal is supported
# doModAssess: whether or not to calculate CPO, DIC, and WAIC
# previousFit: a previous INLA model fit used to initialize optimization
# fixedParameters: A list of parameters to fix in the model rather than infer. 
#                  Contains some of all of the elements: spde$effRange, 
#                  spde$margVar, familyPrec, clusterPrec, beta (NOT TESTED)
# experimentalMode: Whether to use INLA variational inference tools (NOT TESTED)
# bwRepel: does nothing for this model
# 
# Outputs:
# INLA model, predictions, summary statistics, input data, posterior draws, etc.
fitSPDEsimDat = function(wellDat, seismicDat, 
                         predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
                         control.fixed = list(prec=list(default=0, X2=1/.5^2, X3=1), mean=list(default=0, X2=1)), 
                         transform=logit, invTransform=expit, 
                         mesh=getSPDEmeshSimStudy(), prior=getSPDEprior(mesh), 
                         addKDE=FALSE, esthFromSeismic=TRUE, kde.args=NULL, pProcMethod=c("kde", "inlabru"), 
                         significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                         nPostSamples=1000, useSeismic=TRUE, verbose=FALSE, seed=NULL, 
                         family="normal", doModAssess=FALSE, previousFit=NULL, 
                         fixedParameters=NULL, experimentalMode=FALSE, addLogProbs=FALSE, 
                         logProbsNoRep=NULL, bwRepel=NULL) {
  
  # set defaults
  # family = match.arg(family)
  pProcMethod = match.arg(pProcMethod)
  
  if(pProcMethod != "kde") {
    stop("currently only kde method is supported")
  }
  
  if(addKDE) {
    # add log kernel density estimator as a covariate
    
    if(esthFromSeismic) {
      # estimate bandwidth via a variogram estimator on the seismic data
      require(gstat)
      require(sf)
      # transformedGridEsts = transform(seismicDat[,3])
      varioMod = vgm(model="Gau")
      varioDat = data.frame(gridEsts=seismicDat[,3], x=seismicDat$east, y=seismicDat$north)
      # require(sp)
      # # coordinates(varioDat) = ~x+y
      # empVario = variogram(transformedGridEsts~x+y, data=varioDat)
      # varioFit = fit.variogram(empVario, varioMod)
      # 
      varioDatSF = sf::st_as_sf(varioDat, coords=c("x", "y"))
      empVario = variogram(gridEsts~1, data=varioDatSF)
      varioFit = fit.variogram(empVario, varioMod)
      hEst = varioFit$range^2
      
      if(is.null(kde.args)) {
        kde.args = list(H=diag(rep(hEst, 2)))
      } else {
        kde.args$H = diag(rep(hEst, 2))
      }
    }
    
    # fit the density estimator
    fitPProc = getDensityDueToWells(seismicDat=seismicDat, wellDat=wellDat, 
                                    predPts=as.matrix(wellDat[,1:2]), method=pProcMethod, 
                                    centerScale=TRUE, kde.args=kde.args)
    wGrid = fitPProc$wellEffectGrid
    wObs = fitPProc$wellEffectPredPts
  } else if(addLogProbs) {
    wGrid = logProbsNoRep
    wObs = logProbsNoRep[wellDat$gridI]
    
    control.fixed$prec$X3 = NULL
  } else {
    wGrid = NULL
    wObs = NULL
    
    control.fixed$prec$X3 = NULL
  }
  
  # construct prediction points and covariates
  predPts = matrix(unlist(seismicDat[,1:2]), ncol=2)
  
  if(!useSeismic) {
    xPred = cbind(1, wGrid)
    
    # construct well data covariates
    xObs = cbind(1, wObs)
  } else {
    xPred = cbind(1, transform(seismicDat$seismicEst), wGrid)
    
    # interpolate seismic data to the well points
    wellSeismicEsts = bilinearInterp(as.matrix(wellDat[,1:2]), seismicDat, 
                                     transform=transform, invTransform=invTransform)
    
    # construct well data covariates
    xObs = cbind(1, transform(wellSeismicEsts), wObs)
  }
  
  # set observations
  obsValues = wellDat$volFrac
  obsCoords = cbind(wellDat$east, wellDat$north)
  
  fitSPDE(obsCoords=obsCoords, obsValues=obsValues, xObs=xObs, 
          predCoords=predPts, xPred=xPred, control.fixed=control.fixed, 
          transform=transform, invTransform=invTransform, 
          mesh=mesh, prior=prior, addNugToPredCoords=FALSE, 
          significanceCI=significanceCI, int.strategy=int.strategy, strategy=strategy, 
          nPostSamples=nPostSamples, verbose=verbose, link=link, seed=seed, 
          family=family, doModAssess=doModAssess, previousFit=previousFit, 
          fixedParameters=fixedParameters, experimentalMode=experimentalMode)
}


# function for fitting the SPDE model to data
# 
# Inputs:
# obsCoords: matrix with easting and northing columns for observations
# xObs: matrix with intercept and covariate information for observations
# predCoords: data.frame with easting and northing columns for prediction grid
# xPred: matrix with intercept and covariate information for prediction grid
# transform: how to transform obsValues prior to modeling
# invTransform: inverse of transform. Used to backtransform predictions prior to 
#               aggregation over domain
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
# customFixedI: if not NULL, a vector of length ncol(xObs) determining a custom 
#               linear combination of fixed effects to add to SPDE effect in 
#               order to make a custom set of predictions
# previousFit: a previous INLA model fit used to initialize optimization
# fixedParameters: A list of parameters to fix in the model rather than infer. 
#                  Contains some of all of the elements: spde$effRange, 
#                  spde$margVar, familyPrec, clusterPrec, beta (NOT TESTED)
# experimentalMode: Whether to use INLA variational inference tools (NOT TESTED)
# 
# Outputs:
# INLA model, predictions, summary statistics, input data, posterior draws, etc.
fitSPDE = function(obsCoords, obsValues, xObs=matrix(rep(1, length(obsValues)), ncol=1), 
                   predCoords, xPred = matrix(rep(1, nrow(predCoords)), ncol=1), 
                   control.fixed = list(prec=list(default=0), mean=list(default=0)), 
                   transform=I, invTransform=I, 
                   mesh=getSPDEmesh(obsCoords), prior=getSPDEprior(mesh), 
                   significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                   nPostSamples=1000, verbose=TRUE, link=1, seed=NULL, 
                   family=c("normal", "binomial", "betabinomial"), 
                   doModAssess=FALSE, customFixedI=NULL, 
                   previousFit=NULL, addNugToPredCoords=TRUE, 
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
  # control.family = list(hyper = list(prec = list(prior="loggamma", param=c(1000,10))))
  control.family = list(hyper = list(prec = list(prior="loggamma", param=c(1.116960395, 0.009827238))))
  
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
    
    modeControl$result = previousFit
    modeControl$theta = previousFit$mode$theta
    modeControl$x = previousFit$mode$x
    modeControl$restart = TRUE
    modeControl$fixed = FALSE
  }
  
  # Set distributional quantiles we're interested in: median and based on significanceCI
  allQuantiles = c(0.5, (1-significanceCI) / 2, 1 - (1-significanceCI) / 2)
  
  # fixed effect priors: are they improper or not?
  controlFixed = c(list(quantiles=allQuantiles), control.fixed)
  
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
    hyperparNames = names(postSamples[[1]]$hyperpar)
    nuggetVars = sapply(postSamples, function(x) {1 / x$hyperpar[which(hyperparNames == "Precision for the Gaussian observations")]})
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
  
  # now make custom predictions
  if(!is.null(customFixedI)) {
    customFixedMat = diag(customFixedI)
    customFixedPredMat = xPred  %*% customFixedMat %*% latentMat[fixedIndices,]
    
    customPredMat = customFixedPredMat + spatialPredMat
    
    if(!is.null(fixedParameters$beta)) {
      stop("setting customFixedI and fixedParameters simultaneously not supported")
    }
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
  
  # now make custom predictions
  if(!is.null(customFixedI)) {
    customFixedObsMat = xObs  %*% customFixedMat %*% latentMat[fixedIndices,]
    
    customObsMat = customFixedObsMat + spatialObsMat
    
    if(!is.null(fixedParameters$beta)) {
      stop("setting customFixedI and fixedParameters simultaneously not supported")
    }
  }
  
  # add in nugget if necessary
  if(family == "normal") {
    # get cluster effect variance
    clusterVars = nuggetVars
    
    if(addNugToPredCoords) {
      predMatNugget = predMat + matrix(rnorm(length(predMat), sd=rep(sqrt(clusterVars), each=nrow(predMat))), nrow=nrow(predMat))
    } else {
      predMatNugget = NULL
    }
    obsMatNugget = obsMat + matrix(rnorm(length(obsMat), sd=rep(sqrt(clusterVars), each=nrow(obsMat))), nrow=nrow(obsMat))
  } else {
    stop("family not supported")
  }
  
  # back transform predictions
  obsMat = invTransform(obsMat)
  obsMatNugget = invTransform(obsMatNugget)
  predMat = invTransform(predMat)
  if(addNugToPredCoords) {
    predMatNugget = invTransform(predMatNugget)
  } else {
    predMatNugget = NULL
  }
  
  if(!is.null(customFixedI)) {
    customObsMat = invTransform(customObsMat)
    customPredMat = invTransform(customPredMat)
  }
  
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
  
  if(addNugToPredCoords) {
    predNuggetEst = predEst
    predNuggetSDs = apply(predMatNugget, 1, sd)
    predNuggetLower = apply(predMatNugget, 1, quantile, probs=(1-significanceCI)/2)
    predNuggetMedian = apply(predMatNugget, 1, median)
    predNuggetUpper = apply(predMatNugget, 1, quantile, probs=1-(1-significanceCI)/2)
  } else {
    predNuggetEst = NULL
    predNuggetSDs = NULL
    predNuggetLower = NULL
    predNuggetMedian = NULL
    predNuggetUpper = NULL
  }
  
  if(!is.null(customFixedI)) {
    customObsEst = rowMeans(customObsMat)
    customObsSDs = apply(customObsMat, 1, sd)
    customObsLower = apply(customObsMat, 1, quantile, probs=(1-significanceCI)/2)
    customObsMedian = apply(customObsMat, 1, median)
    customObsUpper = apply(customObsMat, 1, quantile, probs=1-(1-significanceCI)/2)
    
    customPredEst = rowMeans(customPredMat)
    customPredSDs = apply(customPredMat, 1, sd)
    customPredLower = apply(customPredMat, 1, quantile, probs=(1-significanceCI)/2)
    customPredMedian = apply(customPredMat, 1, median)
    customPredUpper = apply(customPredMat, 1, quantile, probs=1-(1-significanceCI)/2)
  } else {
    customObsEst = NULL
    customObsSDs = NULL
    customObsLower = NULL
    customObsMedian = NULL
    customObsUpper = NULL
    
    customPredEst = NULL
    customPredSDs = NULL
    customPredLower = NULL
    customPredMedian = NULL
    customPredUpper = NULL
  }
  
  # get aggregate summary statistics over prediction domain
  predAggMat = colMeans(predMat)
  predAggEst = mean(predAggMat)
  predAggSDs = sd(predAggMat)
  predAggLower = quantile(predAggMat, probs=(1-significanceCI)/2)
  predAggMedian = quantile(predAggMat, probs=.5)
  predAggUpper = quantile(predAggMat, probs=1-(1-significanceCI)/2)
  
  if(length(xPred) != 0) {
    interceptSummary=mod$summary.fixed[1,1:5]
    fixedEffectSummary=mod$summary.fixed[,1:5]
  } 
  else {
    interceptSummary = matrix(rep(0, 5), nrow=1)
    fixedEffectSummary = mod$summary.fixed
  }
  hyperNames = row.names(mod$summary.hyperpar)
  
  rangeSummary=mod$summary.hyperpar[which(hyperNames == "Range for field"),1:5]
  spatialSDSummary = mod$summary.hyperpar[which(hyperNames == "Stdev for field"),1:5]
  
  # get posterior hyperparameter samples and transform them as necessary
  hyperMat = sapply(postSamples, function(x) {x$hyperpar})
  hyperNames = row.names(hyperMat)
  if(family == "normal") {
    clusterVarI = which(hyperNames == "Precision for the Gaussian observations")
    spatialRangeI = which(hyperNames == "Range for field" )
    spatialSDI = which(hyperNames == "Stdev for field")
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
      summaryHyperNames = row.names(parameterSummaryTable)
      sdSummary=parameterSummaryTable[summaryHyperNames == "errorSD",]
      varSummary=parameterSummaryTable[summaryHyperNames == "errorVar",]
      rangeSummary=parameterSummaryTable[summaryHyperNames == "spatialRange",]
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
       predEst=predEst, predSDs=predSDs, predLower=predLower, predMedian=predMedian, predUpper=predUpper, predAggMat=predAggMat, 
       predAggEst=predAggEst, predAggSDs=predAggSDs, predAggLower=predAggLower, predAggMedian=predAggMedian, predAggUpper=predAggUpper, 
       mesh=mesh, prior=prior, stack=stack.full, 
       interceptSummary=interceptSummary, fixedEffectSummary=fixedEffectSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, 
       parameterSummaryTable=parameterSummaryTable, 
       fixedEffectDraws=latentMat[fixedIndices,], 
       spatialPredMat=spatialPredMat, fixedPredMat=fixedPredMat, 
       spatialObsMat=spatialObsMat, fixedObsMat=fixedObsMat, 
       obsMat=obsMat, obsMatNugget=obsMatNugget, predMat=predMat, predMatNugget=predMatNugget, 
       customObsEst=customObsEst, customObsSDs=customObsSDs, customObsLower=customObsLower, 
       customObsMedian=customObsMedian, customObsUpper=customObsUpper, 
       customPredEst=customPredEst, customPredSDs=customPredSDs, 
       customPredLower=customPredLower, customPredMedian=customPredMedian, 
       customPredUpper=customPredUpper, 
       hyperMat=hyperMat, timings=timings, sigmaEpsilonDraws=sqrt(clusterVars))
}







