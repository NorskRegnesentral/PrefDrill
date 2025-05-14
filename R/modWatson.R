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
# improperCovariatePrior: if TRUE, N(0, infty) prior on covariates (aside from 
#                         intercept, which already has this prior)
# fixedParameters: A list of parameters to fix in the model rather than infer. 
#                  Contains some of all of the elements: spde$effRange, 
#                  spde$margVar, familyPrec, clusterPrec, beta (NOT TESTED)
# experimentalMode: Whether to use INLA variational inference tools (NOT TESTED)
# 
# Outputs:
# INLA model, predictions, summary statistics, input data, posterior draws, etc.
fitWatsonSimDat = function(wellDat, seismicDat, 
                           predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
                           pseudoCoords=predGrid, 
                           control.fixed = list(prec=list(default=0, X.pp21=1/.5^2, X.y2=1/.5^2), mean=list(default=0, X.pp21=1, X.y2=1)), 
                           transform=logit, invTransform=expit, 
                           mesh=getSPDEmeshSimStudy(), prior=getSPDEprior(mesh), 
                           repelDist=500, sharedInt=FALSE, 
                           significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                           nPostSamples=1000, verbose=FALSE, seed=123, 
                           doModAssess=FALSE, previousFit=NULL, customFixedI=NULL, 
                           quadratureMethod=c("pseudoSites", "mesh"), 
                           fixedParameters=NULL, fixedRepelAmt=500, 
                           experimentalMode=FALSE) {
  
  # set defaults
  # family = match.arg(family)
  quadratureMethod = match.arg(quadratureMethod)
  
  pseudoOnSeismicGrid = all.equal(c(pseudoCoords[,1], pseudoCoords[,2]), c(seismicDat[,1], seismicDat[,2]))
  
  # construct prediction points and covariates
  predPts = matrix(unlist(seismicDat[,1:2]), ncol=2)
  xPred = cbind(1, transform(seismicDat$seismicEst))
  
  # interpolate seismic data to the well points
  wellSeismicEsts = bilinearInterp(wellDat[,1:2], seismicDat, 
                                   transform=transform, invTransform=invTransform)
  
  # construct well data covariates
  xObs = cbind(1, transform(wellSeismicEsts))
  
  # set observations
  obsValues = wellDat$volFrac
  obsCoords = cbind(wellDat$east, wellDat$north)
  
  # interpolate seismic data to the pseudo sites
  if(!pseudoOnSeismicGrid) {
    pseudoSeismicEsts = bilinearInterp(pseudoCoords, seismicDat, 
                                       transform=transform, invTransform=invTransform)
  } else {
    pseudoSeismicEsts = seismicDat[,3]
  }
  
  # construct covariates at pseudo site locations
  xPseudo = cbind(1, transform(pseudoSeismicEsts))
  
  fitWatson(obsCoords=obsCoords, obsValues=obsValues, xObs=xObs, 
            pseudoCoords=pseudoCoords, xPseudo=xPseudo, 
            predCoords=predPts, xPred=xPred, fixedRepelAmt=fixedRepelAmt, 
            control.fixed=control.fixed, 
            transform=transform, invTransform=invTransform, repelDist=repelDist, 
            sharedInt=sharedInt, mesh=mesh, prior=prior, significanceCI=significanceCI, 
            int.strategy=int.strategy, strategy=strategy, nPostSamples=nPostSamples, 
            verbose=verbose, link=link, seed=seed, doModAssess=doModAssess, 
            customFixedI=customFixedI, quadratureMethod=quadratureMethod, 
            previousFit=previousFit, 
            fixedParameters=fixedParameters, experimentalMode=experimentalMode)
}

# function for fitting the Watson et al. model to data
# 
# Inputs:
# obsCoords: matrix with easting and northing columns for observations
# xObs: matrix with intercept and covariate information for observations (NOT 
#       INCLUDING REPULSION COVARIATES!)
# predCoords: data.frame with easting and northing columns for prediction grid
# xPred: matrix with intercept and covariate information for prediction grid 
#        (NOT INCLUDING REPULSION COVARIATES!)
# pseudoCoords: matrix with easting and northing columns for pseudo site coords
# xPseudo: matrix with intercept and covariate information for pseudo sites (NOT 
#          INCLUDING REPULSION COVARIATES!)
# covPriors: list of lists. Contains priors for the non-intercept covariates. If 
#            length 1, repeats priors for all non-intercept covariates.
# transform: how to transform obsValues prior to modeling
# invTransform: inverse of transform. Used to backtransform predictions prior to 
#               aggregation over domain
# repelDist: radial bandwidth for cylindrical radial basis function controling 
#            repel distance
# sharedInt: if TRUE, includes just 1 shared intercept per iteration in the 
#            point process model. Otherwise, includes separate intercepts, 1 for 
#            each iteration, in the point process model
# mesh: SPDE mesh
# prior: SPDE prior
# significanceCI: the credible level of the CIs (e.g., .8 means 80% CI)
# int.strategy: inla integration strategy
# strategy: inla strategy
# nPostSamples: number posterior draws
# verbose: verbose argument to inla
# link: link=1 is the canonical link in inla
# seed: random seed. Not set if NULL
# doModAssess: whether or not to calculate CPO, DIC, and WAIC
# customFixedI: if not NULL, a vector of length ncol(xObs) determining a custom 
#               linear combination of fixed effects to add to SPDE effect in 
#               order to make a custom set of predictions
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
fitWatson = function(obsCoords, obsValues, xObs=matrix(rep(1, length(obsValues)), ncol=1), 
                     pseudoCoords=makePseudositesRect(), xPseudo=matrix(rep(1, nrow(pseudoCoords)), nrow=1), 
                     predCoords, xPred = matrix(rep(1, nrow(predCoords)), ncol=1), 
                     control.fixed = list(prec=list(default=0), mean=list(default=0)), 
                     transform=I, invTransform=I, repelDist=10, sharedInt=FALSE, 
                     mesh=getSPDEmeshSimStudy(), prior=getSPDEprior(mesh), 
                     significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                     nPostSamples=1000, verbose=TRUE, link=1, seed=NULL, 
                     doModAssess=FALSE, customFixedI=NULL, quadratureMethod=c("pseudoSites", "mesh"), 
                     previousFit=NULL, fixedRepelAmt=NULL, 
                     fixedParameters=NULL, experimentalMode=FALSE) {
  
  startTime = proc.time()[3]
  if(!is.null(seed))
    set.seed(seed)
  
  quadratureMethod = match.arg(quadratureMethod)
  
  if(experimentalMode) {
    if(strategy != "gaussian") {
      stop("only gaussian integration is possible if experimentalMode is TRUE")
    }
  }
  
  if(!is.null(fixedParameters)) {
    stop("fixedParameters not currently supported for Watson et al. model")
  }
  
  # check if an intercept is included as a parameter
  if(!is.null(xObs)) {
    if(all(xObs[,1] == 1)) {
      hasInt = TRUE
    } else {
      hasInt = FALSE
    }
  }
  
  if(!hasInt) {
    stop("currently, model must include an intercept")
  }
  
  if(!is.null(customFixedI)) {
    stop("customFixedI not currently supported")
  }
  
  if(quadratureMethod == "mesh") {
    stop("mesh quadrature currently not supported")
  }
  
  if(!all.equal(c(unlist(pseudoCoords)), c(unlist(predCoords)))) {
    stop("must debug this case")
    pseudoArePred = FALSE
  } else {
    pseudoArePred = TRUE
  }
  
  # begin defining the model
  startTimeDefineModel = proc.time()[3]
  
  nObs = nrow(obsCoords)
  nPseudoPerIter = nrow(pseudoCoords)
  
  # fix SPDE model parameters if necessary
  if(!is.null(fixedParameters)) {
    prior = getSPDEModelFixedPar(mesh, effRange=fixedParameters$spde$effRange, 
                                 margVar=fixedParameters$spde$margVar)
  }
  
  # set family prior
  control.family.gaussian = list(hyper = list(prec = list(prior="loggamma", param=c(1000,10))))
  
  if(!is.null(fixedParameters$familyPrec)) {
    # fix the family precision parameter on INLA's latent scale
    control.family.gaussian = list(initial=log(fixedParameters$familyPrec), fixed=TRUE)
  }
  
  # set preferentiality parameter prior
  prefPrior = list(prior="gaussian", param=c(0, 4))
  
  # if there is an intercept, we need to remove it and instead include a 
  # separate intercept for each iteration in the PProc model. Otherwise, no need 
  # to change xObs.pp. Will later expand xObs.pp to account for pseudo-sites
  if(hasInt && !sharedInt) {
    xObs.pp = cbind(Diagonal(n=nObs, x = rep(1, nrow=nrow(xObs))), xObs[,-1])
    xPseudoModInt = matrix(xPseudo[,-1], nrow=nrow(xPseudo))
  } else {
    xObs.pp = xObs
    xPseudoModInt = xPseudo
  }
  
  # make sure matrices are sparse
  xObs.pp = Matrix(xObs.pp, nrow=nrow(xObs.pp), sparse = TRUE)
  xPseudoModInt = Matrix(xPseudoModInt, nrow=nrow(xPseudo), sparse = TRUE)
  
  # construct A matrix for observations (response)
  AObs = inla.spde.make.A(mesh, loc = obsCoords)
  
  # construct A matrix for observations (Bernoulli regression)
  for(i in 1:nObs) {
    thisA.pp = inla.spde.make.A(mesh, loc = as.matrix(rbind(matrix(obsCoords[i,], nrow=1), 
                                                            pseudoCoords)))
    
    # expand xPseudoModInt if need be to include different intercepts per iteration. 
    # Also, make sure everything is stored as a sparse matrix
    if(hasInt && sharedInt) {
      thisXpseudo.pp = xPseudoModInt
    } else if(hasInt) {
      thisXpseudo.pp = cbind(matrix(0, nrow=nPseudoPerIter, ncol=nObs), xPseudoModInt)
      thisXpseudo.pp[,i] = 1
      thisXpseudo.pp = Matrix(thisXpseudo.pp, nrow=nrow(thisXpseudo.pp), sparse = TRUE)
    }
    
    # combine observation and pseudo-site covariates in covariates for this iteration
    thisX.pp = rbind(matrix(xObs.pp[i,], nrow=1), 
                     thisXpseudo.pp)
    
    # get repulsion covariate as a separate design matrix to it's easier to fix
    if(i == 1) {
      thisRepelMat = Matrix(matrix(rep(0, 1+nPseudoPerIter), ncol=1), ncol=1, sparse=TRUE)
    } else {
      thisRepelMat = getRepulsionCov(rbind(obsCoords[i,], pseudoCoords), obsCoords[1:(i-1),], 
                                     repelDist=repelDist, returnSparse=TRUE)
    }
    # thisX.pp = cbind(thisX.pp, thisRepelMat)
    
    if(i == 1) {
      A.pp = thisA.pp
      X.pp = thisX.pp
      X.pp.rep = thisRepelMat
    } else {
      A.pp = rbind(A.pp, thisA.pp)
      X.pp = rbind(X.pp, thisX.pp)
      X.pp.rep = rbind(X.pp.rep, thisRepelMat)
    }
  }
  
  # calculate fixed PProc offset for repulsion if specified
  if(!is.null(fixedRepelAmt)) {
    repelOffset = as.numeric(fixedRepelAmt * X.pp.rep)
  } else {
    repelOffset = rep(0, nrow(X.pp.rep))
  }
  
  # construct A matrix for predictions
  APred = inla.spde.make.A(mesh, loc = predCoords)
  
  # make inla stack
  ys = transform(obsValues)
  rs = rep(c(1, rep(0, nrow(pseudoCoords))), nObs)
  m = ncol(AObs) # number of basis elements
  nR = length(rs)
  nPreds = nrow(predCoords)
  latticeInds = 1:m
  cluster = 1:nObs
  
  # NOTE: INLA doesn't support sparse design matrices for some reason???
  X.ppDense = as.matrix(X.pp)
  
  # browser() # separate seismic covariate from others for a separate prior
  
  # construct the observation stack 
  stack.y = inla.stack(data = list(y=cbind(ys, NA), offset=rep(0, nObs)),
                       A = list(AObs, 1),
                       effects = list(field.y=latticeInds, X.y=xObs),
                       tag = "y",
                       remove.unused=FALSE)
  
  stack.pp = inla.stack(data = list(y=cbind(NA, rs), offset=repelOffset), 
                        A = list(A.pp, 1),
                        effects = list(field.pp=latticeInds, X.pp=X.ppDense),
                        tag = "pp", 
                        remove.unused=FALSE)
  
  endTimeDefineModel = proc.time()[3]
  totalTimeDefineModel = endTimeDefineModel - startTimeDefineModel
  
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
  controlFixed = c(list(quantiles=allQuantiles), control.fixed)
  
  # construct the stack
  stack.full = inla.stack(stack.y, stack.pp)
  stackDat = inla.stack.data(stack.full, spde=prior)
  
  # see: inla.doc("loggamma")
  # shape=.1, scale=10 for unit mean, variance 100 prior
  
  # setup the formula
  thisFormula = paste0("y ~ -1 + f(field.y, model=prior)", 
                       "+ f(field.pp, copy='field.y', fixed=FALSE, hyper=list(beta=prefPrior))")
  if(!is.null(xObs)) {
    
    thisFormula = paste0(thisFormula, " + X.y + X.pp")
  }
  
  thisFormula = as.formula(thisFormula)
  
  startModelFitTime = proc.time()[3]
  
  mod = inla(
    thisFormula, 
    data = stackDat, offset=offset, 
    control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles), 
    family=c("gaussian", "poisson"), verbose=verbose, control.inla=control.inla, 
    control.compute=list(config=TRUE, cpo=doModAssess, dic=doModAssess, waic=doModAssess), 
    control.mode=modeControl, 
    control.fixed=controlFixed, 
    control.family=list(control.family.gaussian, list())
  )
  
  endModelFitTime = proc.time()[3]
  totalModelFitTime = endModelFitTime - startModelFitTime
  
  # get prediction covariates
  predRepelMat = getRepulsionCov(predCoords, obsCoords, 
                                 repelDist=repelDist, returnSparse=TRUE)
  if(!is.null(fixedRepelAmt)) {
    predOffset.pp = predRepelMat * fixedRepelAmt
  } else {
    predOffset.pp = rep(0, nrow(predCoords))
  }
  
  # generate samples from posterior
  startTimePosteriorSampling = proc.time()[3]
  postSamples = inla.posterior.sample(nPostSamples, mod)
  endTimePosteriorSampling = proc.time()[3]
  totalTimePosteriorSampling = endTimePosteriorSampling - startTimePosteriorSampling
  
  latentMat = sapply(postSamples, function(x) {x$latent})
  
  hyperparNames = names(postSamples[[1]]$hyperpar)
  nuggetVars = sapply(postSamples, function(x) {1 / x$hyperpar[which(hyperparNames == "Precision for the Gaussian observations")]})
  
  latentVarNames = rownames(postSamples[[1]]$latent)
  field.yIndices = which(grepl("field.y", latentVarNames))
  fixed.yIndices = which(grepl("X.y", latentVarNames))
  field.ppIndices = which(grepl("field.pp", latentVarNames))
  fixed.ppIndices = which(grepl("X.pp", latentVarNames))
  
  if(!sharedInt) {
    # remove intercepts (1 from each iteration)
    fixed.ppIndices = fixed.ppIndices[-(1:nObs)]
  } else {
    browser() # make sure we're not supposed to do anything in this case
  }
  
  # integrate intensity to calculate point process model intercept. 
  # Approximate integral via sum over pseudosites
  if(quadratureMethod == "mesh") {
    stop("mesh quadrature not supported")
  } else {
    # fixed part
    if(!sharedInt && hasInt) {
      fixedPseudoMat = matrix(X.ppDense[2:(1+nPseudoPerIter),-(1:nObs)], nrow=nPseudoPerIter) %*% matrix(latentMat[fixed.ppIndices,], ncol=nPostSamples)
    } else {
      fixedPseudoMat = X.ppDense[2:(1+nPseudoPerIter),] %*% latentMat[fixed.ppIndices,]
    }
    
    # spatial part
    APseudo = A.pp[2:(1+nPseudoPerIter),]
    spatialPseudoMat = APseudo %*% latentMat[field.ppIndices,]
    
    # offset part
    pseudoOffset = X.pp.rep[2:(1+nPseudoPerIter)]
    
    quadMat = sweep(fixedPseudoMat + spatialPseudoMat, 1, pseudoOffset, "+")
  }
  totals = colSums(exp(quadMat))
  intercept.pp = -log(totals)
  
  # generate predictions: first without nugget then add the nugget in
  if(hasInt) {
    xPredNoInt = matrix(xPred[,-1], nrow=nrow(xPred))
  } else {
    xPredNoInt = xPred
  }
  
  if(length(xPred) != 0) {
    fixedPredMat.y = xPred  %*% latentMat[fixed.yIndices,]
    if(!sharedInt && hasInt) {
      fixedPredMat.pp = xPredNoInt %*% latentMat[fixed.ppIndices,]
    } else {
      fixedPredMat.pp = xPred %*% latentMat[fixed.ppIndices,]
    }
    fixedPredMat.pp = fixedPredMat.pp
  } else {
    fixedPredMat.y = 0
    fixedPredMat.pp = 0
  }
  fixedPredMat.pp = sweep(sweep(fixedPredMat.pp, 2, intercept.pp, "+"), 1, predOffset.pp, "+")
  
  spatialPredMat.y = APred %*% latentMat[field.yIndices,]
  spatialPredMat.pp = APred %*% latentMat[field.ppIndices,]
  
  predMat.y = fixedPredMat.y + spatialPredMat.y
  predMat.pp = fixedPredMat.pp + spatialPredMat.pp
  
  # fixed parameters (aside from repulsion) currently not supported
  # if(!is.null(offsetPred)) {
  #   predMat.y = sweep(predMat.y, 1, offsetPred.y, "+")
  # }
  
  # now make custom predictions (for responses)
  if(!is.null(customFixedI)) {
    browser()
    customFixedMat = diag(customFixedI)
    customFixedPredMat = xPred  %*% customFixedMat %*% latentMat[fixed.yIndices,]
    
    customPredMat = customFixedPredMat + spatialPredMat.y
    
    if(!is.null(fixedParameters$beta)) {
      stop("setting customFixedI and fixedParameters simultaneously not supported")
    }
  }
  
  # do the same for the observations
  if(length(xObs) != 0) {
    fixedObsMat.y = xObs  %*% matrix(latentMat[fixed.yIndices,], ncol=nPostSamples)
    if(!sharedInt && hasInt) {
      fixedObsMat.pp = xObs[,-1] %*% matrix(latentMat[fixed.ppIndices,], ncol=nPostSamples)
    } else {
      fixedObsMat.pp = xObs %*% matrix(latentMat[fixed.ppIndices,], ncol=nPostSamples)
    }
  } else {
    fixedObsMat.y = 0
    fixedObsMat.pp = 0
  }
  fixedObsMat.pp = sweep(fixedObsMat.pp, 2, intercept.pp, "+")
  
  # get repulsion design matrix at observations (at the time they were sampled). 
  # Add effect to observation predictions (not necessary if repulsion isn't 
  # fixed, since already included in fixedObsMat.pp)
  obsRepelMat = getRepulsionCovAtObs(obsCoords=obsCoords, repelDist=repelDist, returnSparse=TRUE)
  if(!is.null(fixedRepelAmt)) {
    offsetObs = obsRepelMat %*% fixedRepelAmt
    fixedObsMat.pp = sweep(fixedObsMat.pp, 1, offsetObs, "+")
  }
  
  spatialObsMat.y = AObs %*% latentMat[field.yIndices,]
  spatialObsMat.pp = AObs %*% latentMat[field.ppIndices,]
  
  obsMat.y = fixedObsMat.y + spatialObsMat.y
  obsMat.pp = fixedObsMat.pp + spatialObsMat.pp
  
  # currently, fixed parameters aside from fixed repulsion not supported
  # if(!is.null(offsetEst)) {
  #   obsMat = sweep(obsMat, 1, offsetEst, "+")
  # }
  
  # now make custom predictions
  if(!is.null(customFixedI)) {
    browser()
    if(!sharedInt && hasInt) {
      customFixedObsMat = xObs[,-1]  %*% customFixedMat %*% latentMat[fixedIndices,]
    } else {
      customFixedObsMat = xObs  %*% customFixedMat %*% latentMat[fixedIndices,]
    }
    
    customObsMat = customFixedObsMat + spatialObsMat
    
    if(!is.null(fixedParameters$beta)) {
      stop("setting customFixedI and fixedParameters simultaneously not supported")
    }
  }
  
  # add in nugget if necessary
  
  # get cluster effect variance
  clusterVars = nuggetVars
  
  predMat.yNugget = predMat.y + matrix(rnorm(length(predMat.y), sd=rep(sqrt(clusterVars), each=nrow(predMat.y))), nrow=nrow(predMat.y))
  obsMat.yNugget = obsMat.y + matrix(rnorm(length(obsMat.y), sd=rep(sqrt(clusterVars), each=nrow(obsMat.y))), nrow=nrow(obsMat.y))
  
  # back transform predictions
  obsMat.y = invTransform(obsMat.y)
  obsMat.yNugget = invTransform(obsMat.yNugget)
  predMat.y = invTransform(predMat.y)
  predMat.yNugget = invTransform(predMat.yNugget)
  
  if(!is.null(customFixedI)) {
    browser()
    customObsMat = invTransform(customObsMat)
    customPredMat = invTransform(customPredMat)
  }
  
  # get summary statistics
  obs.yEst = rowMeans(obsMat.y)
  obs.ySDs = apply(obsMat.y, 1, sd)
  obs.yLower = apply(obsMat.y, 1, quantile, probs=(1-significanceCI)/2)
  obs.yMedian = apply(obsMat.y, 1, median)
  obs.yUpper = apply(obsMat.y, 1, quantile, probs=1-(1-significanceCI)/2)
  
  obs.ppEst = rowMeans(obsMat.pp)
  obs.ppSDs = apply(obsMat.pp, 1, sd)
  obs.ppLower = apply(obsMat.pp, 1, quantile, probs=(1-significanceCI)/2)
  obs.ppMedian = apply(obsMat.pp, 1, median)
  obs.ppUpper = apply(obsMat.pp, 1, quantile, probs=1-(1-significanceCI)/2)
  
  obs.yNuggetEst = obs.yEst
  obs.yNuggetSDs = apply(obsMat.yNugget, 1, sd)
  obs.yNuggetLower = apply(obsMat.yNugget, 1, quantile, probs=(1-significanceCI)/2)
  obs.yNuggetMedian = apply(obsMat.yNugget, 1, median)
  obs.yNuggetUpper = apply(obsMat.yNugget, 1, quantile, probs=1-(1-significanceCI)/2)
  
  pred.yEst = rowMeans(predMat.y)
  pred.ySDs = apply(predMat.y, 1, sd)
  pred.yLower = apply(predMat.y, 1, quantile, probs=(1-significanceCI)/2)
  pred.yMedian = apply(predMat.y, 1, median)
  pred.yUpper = apply(predMat.y, 1, quantile, probs=1-(1-significanceCI)/2)
  
  pred.ppEst = rowMeans(predMat.pp)
  pred.ppSDs = apply(predMat.pp, 1, sd)
  pred.ppLower = apply(predMat.pp, 1, quantile, probs=(1-significanceCI)/2)
  pred.ppMedian = apply(predMat.pp, 1, median)
  pred.ppUpper = apply(predMat.pp, 1, quantile, probs=1-(1-significanceCI)/2)
  
  pred.yNuggetEst = pred.yEst
  pred.yNuggetSDs = apply(predMat.yNugget, 1, sd)
  pred.yNuggetLower = apply(predMat.yNugget, 1, quantile, probs=(1-significanceCI)/2)
  pred.yNuggetMedian = apply(predMat.yNugget, 1, median)
  pred.yNuggetUpper = apply(predMat.yNugget, 1, quantile, probs=1-(1-significanceCI)/2)
  
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
  
  
  # get aggregate summary statistics over prediction domain (no need for the 
  # aggregate summary statistics for the PProc model)
  pred.yAggMat = colMeans(predMat.y)
  pred.yAggEst = mean(pred.yAggMat)
  pred.yAggSDs = sd(pred.yAggMat)
  pred.yAggLower = quantile(pred.yAggMat, probs=(1-significanceCI)/2)
  pred.yAggMedian = quantile(pred.yAggMat, probs=.5)
  pred.yAggUpper = quantile(pred.yAggMat, probs=1-(1-significanceCI)/2)
  
  if(length(xPred) != 0) {
    fixedEffectSummary=mod$summary.fixed[,1:5]
    if(hasInt) {
      interceptSummary = fixedEffectSummary[1,]
    } else {
      interceptSummary = matrix(rep(0, 5), nrow=1)
    }
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
  
  clusterVarI = which(hyperNames == "Precision for the Gaussian observations")
  spatialRangeI = which(hyperNames == "Range for field.y" )
  spatialSDI = which(hyperNames == "Stdev for field.y")
  prefParI = which(hyperNames == "Beta for field.pp")
  if(!is.matrix(hyperMat)) {
    mat = NULL
  } else {
    mat = apply(hyperMat, 2, function(x) {c(totalVar=x[spatialSDI]^2+1/x[clusterVarI], spatialVar=x[spatialSDI]^2, errorVar=1/x[clusterVarI], 
                                            totalSD=sqrt(x[spatialSDI]^2+1/x[clusterVarI]), spatialSD=x[spatialSDI], errorSD=sqrt(1/x[clusterVarI]), 
                                            spatialRange=x[spatialRangeI], prefPar=x[prefParI])})
  }
  
  hyperNames = c("totalVar", "spatialVar", "errorVar", "totalSD", "spatialSD", "errorSD", "spatialRange", "prefPar")
  
  getSummaryStatistics = function(draws) {
    c(Est=mean(draws, na.rm=TRUE), SD=sd(draws, na.rm=TRUE), 
      Qlower=quantile(probs=(1 - significanceCI) / 2, draws, na.rm=TRUE), 
      Q50=quantile(probs=0.5, draws, na.rm=TRUE), 
      Qupper=quantile(probs=1 - (1 - significanceCI) / 2, draws, na.rm=TRUE))
  }
  
  if(is.matrix(hyperMat)) {
    rownames(mat) = hyperNames
    
    summaryNames = c("Est", "SD", "Qlower", "Q50", "Qupper")
    parameterSummaryTable = t(apply(mat, 1, getSummaryStatistics))
    colnames(parameterSummaryTable) = summaryNames
    
    # separate out default parameter summaries
    summaryHyperNames = row.names(parameterSummaryTable)
    sdSummary=parameterSummaryTable[summaryHyperNames == "errorSD",]
    varSummary=parameterSummaryTable[summaryHyperNames == "errorVar",]
    rangeSummary=parameterSummaryTable[summaryHyperNames == "spatialRange",]
    prefParSummary = parameterSummaryTable[summaryHyperNames == "prefPar",]
  } else {
    parameterSummaryTable = NULL
    sdSummary = NULL
    varSummary = NULL
    rangeSummary = NULL
    overdispersionSummary = NULL
    prefParSummary = NULL
  }
  
  interceptSummary.pp = getSummaryStatistics(intercept.pp)
  
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
       obs.yEst=obs.yEst, obs.ySDs=obs.ySDs, obs.yLower=obs.yLower, obs.yMedian=obs.yMedian, obs.yUpper=obs.yUpper, 
       pred.yEst=pred.yEst, pred.ySDs=pred.ySDs, pred.yLower=pred.yLower, pred.yMedian=pred.yMedian, pred.yUpper=pred.yUpper, pred.yAggMat=pred.yAggMat, 
       pred.yAggEst=pred.yAggEst, pred.yAggSDs=pred.yAggSDs, pred.yAggLower=pred.yAggLower, pred.yAggMedian=pred.yAggMedian, pred.yAggUpper=pred.yAggUpper, 
       obs.ppEst=obs.ppEst, obs.ppSDs=obs.ppSDs, obs.ppLower=obs.ppLower, obs.ppMedian=obs.ppMedian, obs.ppUpper=obs.ppUpper, 
       pred.ppEst=pred.ppEst, pred.ppSDs=pred.ppSDs, pred.ppLower=pred.ppLower, pred.ppMedian=pred.ppMedian, pred.ppUpper=pred.ppUpper, 
       mesh=mesh, prior=prior, stack=stack.full, 
       interceptSummary=interceptSummary, interceptSummary.pp=interceptSummary.pp, 
       fixedEffectSummary=fixedEffectSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, prefParSummary=prefParSummary, 
       parameterSummaryTable=parameterSummaryTable, 
       fixedEffect.yDraws=latentMat[fixed.yIndices,], 
       spatialPredMat.y=spatialPredMat.y, fixedPredMat.y=fixedPredMat.y, 
       spatialObsMat.y=spatialObsMat.y, fixedObsMat.y=fixedObsMat.y, 
       obsMat.y=obsMat.y, obsMat.yNugget=obsMat.yNugget, 
       predMat.y=predMat.y, predMat.yNugget=predMat.yNugget, 
       fixedEffect.ppDraws=latentMat[fixed.ppIndices,], 
       spatialPredMat.pp=spatialPredMat.pp, fixedPredMat.pp=fixedPredMat.pp, 
       spatialObsMat.pp=spatialObsMat.pp, fixedObsMat.pp=fixedObsMat.pp, 
       obsMat.pp=obsMat.pp, predMat.pp=predMat.pp, 
       customObsEst=customObsEst, customObsSDs=customObsSDs, customObsLower=customObsLower, 
       customObsMedian=customObsMedian, customObsUpper=customObsUpper, 
       customPredEst=customPredEst, customPredSDs=customPredSDs, 
       customPredLower=customPredLower, customPredMedian=customPredMedian, 
       customPredUpper=customPredUpper, 
       hyperMat=hyperMat, timings=timings, sigmaEpsilonDraws=sqrt(clusterVars))
}

# get a reasonable default set of pseudosites for a rectangular spatial domain
# 
# inputs: 
# lowerLeft: a 2-vector of form c(x, y) of lower left corner of the rectangular domain
# width: width of domain
# height: height of the domain
# n: target number of pseudosites
# doPlot: whether to plot the resulting mesh for testing purposes
# ...: other arguments to inla.mesh.2d
# 
# Outputs:
# 2 column matrix of eastings and northings of grid points
makePseudositesRect = function(lowerLeft=c(0,0), width=1, height=1, n=60^2, 
                               doPlot=TRUE, ...) {
  
  # a h = w
  a = width/height
  
  # a n_h = n_w
  # n_h n_w = n
  # a n_h^2 = n
  # n_h = sqrt(n/a)
  nh = sqrt(n/a)
  nw = a * nh
  
  # make the 1D grid sequences in each direction
  deltaW = width/nw
  deltaH = height/nh
  ws = seq(deltaW/2, width-deltaW/2, by=deltaW)
  hs = seq(deltaH/2, width-deltaH/2, by=deltaH)
  
  # get the full 2D grid
  pts = make.surface.grid(list(east=ws, north=hs))
  
  # return results
  pts
}

# construct B_pred matrix where B_ij is whether predPt i is within repelDist of 
# any wellPts numbered 1 to (j-1)

# construct B_resp matrix where B_ij = I_ij is whether predPt i is within 
# repelDist of any other wellPts numbered 1 to (j-1)

# get Poisson process likelihood quadrature weights from mesh, domain
# 
# Based on:
# https://becarioprecario.bitbucket.io/spde-gitbook/ch-lcox.html#the-mesh-and-the-weights
# but modified to use sf instead of sp
getQuadWeights = function(mesh=getSPDEmeshSimStudy(), domain=getDomainSimStudy()) {
  require(INLAspacetime)
  
  # calculate the dual (Voronoi) mesh
  dualMesh = INLAspacetime::mesh.dual(mesh, returnclass="list", mc.cores=1)
  
  # calculate the centroids of the Voronoi cells
  centroidXs = sapply(dualMesh, function(x) {mean(x[,1])})
  centroidYs = sapply(dualMesh, function(x) {mean(x[,2])})
  centroids = cbind(centroidXs, centroidYs)
  
  # convert to sf
  dualMesh = lapply(dualMesh, function(x) {
    x = rbind(x, x[1,])
    sf::st_polygon(list(x))
  })
  
  # calclulate mesh weights based on area of polygon (after intersecting with domain)
  weights = sapply(dualMesh, function(x) {
    inter = sf::st_intersection(x, domain)
    sf::st_area(inter)
  })
  
  if(FALSE) {
    # plot dual mesh and weights for testing
    xr = range(mesh$loc[,1])
    yr = range(mesh$loc[,2])
    
    splot(centroids[,1], centroids[,2], weights, pch=19, cex=.1, xlim=xr, ylim=yr)
    for(i in 1:length(dualMesh)) {
      thisPoly = dualMesh[[i]]
      
      plot(thisPoly, add=TRUE, border=rgb(.5, .5, .5, .2), lwd=.3)
    }
    
  }
  
  list(dualMesh=dualMesh, centroids=centroids, weights=weights)
}

# Constructs design matrix for the 1-0 repulsion within some distance. Takes 
# value -1 near repulsion distance of obsCoords, and 0 otherwise
# 
# predCoords: new point locations
# obsCoords: observation locations new points are repelled from
getRepulsionCov = function(predCoords, obsCoords, repelDist=10, returnSparse=FALSE) {
  dists = rdist(predCoords, obsCoords)
  out = matrix(apply(dists, 1, function(x) {-as.numeric(any(x < repelDist))}), ncol=1)
  
  if(returnSparse) {
    out = Matrix(out, sparse = TRUE)
  } else {
    out
  }
}

# same as getRepulsionCov, but gives the values of the repuslion covariate for 
# each observation when it was sampled
getRepulsionCovAtObs = function(obsCoords, repelDist=10, returnSparse=FALSE) {
  dists = rdist(obsCoords)
  diag(dists) = Inf
  repelInd = apply(dists, 1, function(x) {min(which(x < repelDist))})
  out = matrix(-as.numeric(repelInd < 1:nrow(obsCoords)), ncol=1)
  
  if(returnSparse) {
    out = Matrix(out, sparse = TRUE)
  } else {
    out
  }
}












