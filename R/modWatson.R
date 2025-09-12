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
fitWatsonSimDat = function(wellDat, seismicDat, nPseudo=2500, 
                           predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
                           pseudoCoords=getPseudoCoordsSimStudy(maxPts=nPseudo, predGrid=predGrid), 
                           control.fixed = list(prec=list(default=0, X.pp21=1/.5^2, X.y2=1/.5^2), mean=list(default=0, X.pp21=1, X.y2=1)), 
                           transform=logit, invTransform=expit, 
                           mesh=getSPDEmeshSimStudy(), prior=getSPDEprior(mesh), 
                           repelDist=500, sharedInt=FALSE, prefMean=0, 
                           significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                           nPostSamples=1000, verbose=FALSE, seed=NULL, 
                           doModAssess=FALSE, previousFit=NULL, customFixedI=NULL, 
                           quadratureMethod=c("pseudoSites", "mesh"), 
                           fixedParameters=NULL, fixedRepelAmt=500, 
                           addNugToPredCoords=FALSE, getPPres=FALSE, anisFac=1, 
                           experimentalMode=FALSE, bernApprox=FALSE, controlvb=control.vb(), 
                           debugMode=FALSE) {
  
  # set defaults
  # family = match.arg(family)
  quadratureMethod = match.arg(quadratureMethod)
  
  # check if pseudoCoords are the same as the seismicDat coordinates
  if(all.equal(dim(pseudoCoords), dim(seismicDat[,1:2])) == TRUE) {
    pseudoOnSeismicGrid = all.equal(c(pseudoCoords[,1], pseudoCoords[,2]), c(seismicDat[,1], seismicDat[,2]))
  } else {
    pseudoOnSeismicGrid = FALSE
  }
  
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
            control.fixed=control.fixed, obsGridI=wellDat$gridI, 
            transform=transform, invTransform=invTransform, repelDist=repelDist, 
            sharedInt=sharedInt, mesh=mesh, prior=prior, significanceCI=significanceCI, 
            int.strategy=int.strategy, strategy=strategy, nPostSamples=nPostSamples, 
            verbose=verbose, link=link, seed=seed, doModAssess=doModAssess, 
            customFixedI=customFixedI, quadratureMethod=quadratureMethod, prefMean=prefMean, 
            previousFit=previousFit, addNugToPredCoords=addNugToPredCoords, getPPres=getPPres, 
            fixedParameters=fixedParameters, experimentalMode=experimentalMode, bernApprox=bernApprox, 
            controlvb=controlvb, anisFac=anisFac, debugMode=debugMode)
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
# addNugToPredCoords: whether or not to compute results regarding adding nugget 
#                     at predictive coords
# getPPres: Whether or not to compute information on point process estimation
# 
# Outputs:
# INLA model, predictions, summary statistics, input data, posterior draws, etc.
fitWatson = function(obsCoords, obsValues, xObs=matrix(rep(1, length(obsValues)), ncol=1), 
                     pseudoCoords=makePseudositesRect(), obsGridI, 
                     xPseudo=matrix(rep(1, nrow(pseudoCoords)), nrow=1), 
                     predCoords, xPred = matrix(rep(1, nrow(predCoords)), ncol=1), 
                     control.fixed = list(prec=list(default=0), mean=list(default=0)), 
                     transform=I, invTransform=I, repelDist=10, sharedInt=FALSE, 
                     mesh=getSPDEmeshSimStudy(), prior=getSPDEprior(mesh), 
                     significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                     nPostSamples=1000, verbose=TRUE, link=1, seed=NULL, 
                     doModAssess=FALSE, customFixedI=NULL, quadratureMethod=c("pseudoSites", "mesh"), 
                     previousFit=NULL, fixedRepelAmt=NULL, addNugToPredCoords=TRUE, 
                     getPPres=TRUE, prefMean=0, anisFac=1, 
                     fixedParameters=NULL, experimentalMode=FALSE, bernApprox=FALSE, 
                     controlvb=control.vb(), debugMode=FALSE) {
  
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
  
  if(ncol(pseudoCoords) == 3) {
    pseudoPredI = pseudoCoords[,3] # indices of predGrid closest to pseudoCoords
    pseudoCoords = pseudoCoords[,1:2]
  }
  
  # check if pseudoCoords are the same as the prediction coordinates
  if(all.equal(dim(pseudoCoords), dim(predCoords)) == TRUE) {
    # stop("must debug this case")
    pseudoArePred = all.equal(c(pseudoCoords[,1], pseudoCoords[,2]), c(predCoords[,1], predCoords[,2]))
  } else {
    pseudoArePred = FALSE
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
  prefPrior = list(prior="gaussian", param=c(prefMean, 1/4))
  
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
  xObs.pp = Matrix(xObs.pp, sparse = TRUE)
  xPseudoModInt = Matrix(xPseudoModInt, sparse = TRUE)
  
  # construct A matrix for observations (response)
  AObs = inla.spde.make.A(mesh, loc = obsCoords)
  
  # get gridCoords associated with well data and pseudoCoords
  obsGridCoords = predCoords[obsGridI,]
  if(pseudoArePred) {
    pseudoGridCoords = pseudoCoords
  } else {
    pseudoGridCoords = predCoords[pseudoPredI,]
  }
  
  # construct A and design matrices for observations (Bernoulli regression)
  for(i in 1:nObs) {
    thisA.pp = inla.spde.make.A(mesh, loc = as.matrix(rbind(matrix(obsCoords[i,], nrow=1), 
                                                            pseudoCoords)))
    
    # get repulsion covariate as a separate design matrix so it's easier to fix effect size
    if(i == 1) {
      thisRepelMat = Matrix(matrix(rep(0, 1+nPseudoPerIter), ncol=1), sparse=TRUE)
    } else {
      # thisRepelMat = getRepulsionCov(rbind(obsCoords[i,], pseudoCoords), obsCoords[1:(i-1),], 
      #                                repelDist=repelDist, returnSparse=TRUE)
      
      thisRepelMat = getRepulsionCov(rbind(obsGridCoords[i,], pseudoGridCoords), matrix(obsCoords[1:(i-1),], ncol=2), 
                                     repelDist=repelDist, returnSparse=TRUE, anisFac=anisFac)
    }
    
    if(i == 1) {
      A.pp = thisA.pp
      X.pp.rep = thisRepelMat
      
      if(!is.null(xObs)) {
        X.pp = rbind(matrix(xObs[i,], nrow=1), 
                     xPseudo)
      }
      
    } else {
      A.pp = rbind(A.pp, thisA.pp)
      X.pp.rep = rbind(X.pp.rep, thisRepelMat)
      
      if(!is.null(xObs)) {
        X.pp = rbind(X.pp, 
                     matrix(xObs[i,], nrow=1), 
                     xPseudo)
      }
    }
    
    if(debugMode) {
      
      if(i != 1) {
        if(FALSE) {
          # checked this already and it looks good
          browser()
          testGridCoords = rbind(obsGridCoords[i,], pseudoGridCoords)
          testCoords = rbind(obsCoords[i,], pseudoCoords)
          
          pseudoGridEast = sort(unique(pseudoCoords[,1]))
          pseudoGridNorth = sort(unique(pseudoCoords[,2]))
          predGridEast = sort(unique(predCoords[,1]))
          predGridNorth = sort(unique(predCoords[,2]))
          
          squilt(testGridCoords[,1], testGridCoords[,2], as.numeric(thisRepelMat), grid=list(x=predGridEast, y=predGridNorth))
          points(obsCoords[1:(i-1),1], obsCoords[1:(i-1),2], col="red")
          # points(obsGridCoords[1:(i-1),1], obsGridCoords[1:(i-1),2], col="purple")
          points(pseudoCoords[,1], pseudoCoords[,2], col="lightblue", cex=.3)
          points(obsGridCoords[i,1], obsGridCoords[i,2])
          points(obsCoords[i,1], obsCoords[i,2], col="green")
          
          # squilt(obsGridCoords[1:(i-1),1], obsGridCoords[1:(i-1),2], as.numeric(thisRepelMat)[1:(i-1)], grid=list(x=pseudoGridEast, y=pseudoGridNorth))
          # points(obsCoords[1:(i-1),1], obsCoords[1:(i-1),2], col="red")
          # points(obsGridCoords[1:(i-1),1], obsGridCoords[1:(i-1),2], col="purple")
          # points(pseudoCoords[,1], pseudoCoords[,2], col="lightblue", cex=.3)
          # 
          # squilt(testCoords[,1], testCoords[,2], as.numeric(thisRepelMat), grid=list(x=pseudoGridEast, y=pseudoGridNorth))
          # points(obsCoords[1:(i-1),1], obsCoords[1:(i-1),2], col="red")
          # points(obsGridCoords[1:(i-1),1], obsGridCoords[1:(i-1),2], col="purple")
          # points(pseudoCoords[,1], pseudoCoords[,2], col="lightblue", cex=.3)
          # points(predCoords[,1], predCoords[,2], col="blue", pch=".")
        }
      }
    }
  }
  
  # make indices for intercepts
  idxIter = rep(1:nObs, each=1+nPseudoPerIter)
  
  # calculate fixed PProc offset for repulsion if specified
  if(!is.null(fixedRepelAmt)) {
    repelOffset = as.numeric(fixedRepelAmt * X.pp.rep)
  } else {
    repelOffset = rep(0, nrow(X.pp.rep))
  }
  
  if(debugMode) {
    browser()
    
    testGridCoords = rbind(obsGridCoords[i,], pseudoGridCoords)
    testCoords = rbind(obsCoords[i,], pseudoCoords)
    
    pseudoGridEast = sort(unique(pseudoCoords[,1]))
    pseudoGridNorth = sort(unique(pseudoCoords[,2]))
    predGridEast = sort(unique(predCoords[,1]))
    predGridNorth = sort(unique(predCoords[,2]))
    
    endI = length(repelOffset)
    startI = endI - nPseudoPerIter + 1
    
    squilt(pseudoGridCoords[,1], pseudoGridCoords[,2], as.numeric(repelOffset)[startI:endI], 
           grid=list(x=predGridEast, y=predGridNorth))
    points(obsCoords[1:(i-1),1], obsCoords[1:(i-1),2], col="red")
    # points(obsGridCoords[1:(i-1),1], obsGridCoords[1:(i-1),2], col="purple")
    points(pseudoCoords[,1], pseudoCoords[,2], col="lightblue", cex=.3)
    
    range(as.numeric(repelOffset)[obsInds]) # should just be 0's, which it is
    # obsInds = seq(from=1, to=length(repelOffset), by=nPseudoPerIter+1)
    # squilt(obsGridCoords[,1], obsGridCoords[,2], as.numeric(repelOffset)[obsInds], 
    #        grid=list(x=predGridEast, y=predGridNorth))
    # points(obsCoords[,1], obsCoords[,2], col="red")
    # points(obsGridCoords[,1], obsGridCoords[,2], col="purple")
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
  
  # # NOTE: INLA doesn't support sparse design matrices for some reason???
  # X.pp.inla = inla.as.sparse(X.pp)
  # X.ppDense = as.matrix(X.pp)
  
  # browser() # separate seismic covariate from others for a separate prior
  
  # construct the observation stack 
  if(!is.null(xObs)) {
    Alist.y = list(AObs, 1)
    effects.y = list(field.y=latticeInds, X.y=xObs)
    Alist.pp = list(A.pp, 1, 1)
    effects.pp = list(field.pp=latticeInds, idx=idxIter, X.pp=X.pp)
  } else {
    Alist.y = list(AObs)
    effects.y = list(field.y=latticeInds)
    Alist.pp = list(A.pp, 1)
    effects.pp = list(field.pp=latticeInds, idx=idxIter)
  }
  stack.y = inla.stack(data = list(y=cbind(ys, NA), offset=rep(0, nObs)),
                       A = Alist.y,
                       effects = effects.y,
                       tag = "y",
                       remove.unused=FALSE)
  
  # stack.pp = inla.stack(data = list(y=cbind(NA, rs), offset=repelOffset), 
  #                       A = list(A.pp, 1),
  #                       effects = list(field.pp=latticeInds, X.pp=X.ppDense),
  #                       tag = "pp", 
  #                       remove.unused=FALSE)
  stack.pp = inla.stack(data = list(y=cbind(NA, rs), offset=repelOffset), 
                        A = Alist.pp,
                        effects = effects.pp,
                        tag = "pp", 
                        remove.unused=FALSE)
  
  endTimeDefineModel = proc.time()[3]
  totalTimeDefineModel = endTimeDefineModel - startTimeDefineModel
  
  # fit model
  control.inla = list(strategy=strategy, int.strategy=int.strategy, control.vb=controlvb)
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
  
  
  # make prior to convert iid effects to fixed effects. Fix the precision at 0
  prec.prior <- list(prec = list(prior = "gaussian", param = c(0, 1)))
  
  # see: inla.doc("loggamma")
  # shape=.1, scale=10 for unit mean, variance 100 prior
  
  # setup the formula
  # thisFormula = paste0("y ~ -1 + f(field.y, model=prior)", 
  #                      "+ f(field.pp, copy='field.y', fixed=FALSE, hyper=list(beta=prefPrior))")
  # if(!is.null(xObs)) {
  #   
  #   thisFormula = paste0(thisFormula, " + X.y + X.pp")
  # }
  
  thisFormula = paste0("y ~ -1 + f(field.y, model=prior)", 
                       "+ f(field.pp, copy='field.y', fixed=FALSE, hyper=list(beta=prefPrior))", 
                       "+ f(idx, model='iid', hyper=prec.prior, initial = 0, fixed = TRUE)")
  if(!is.null(xObs)) {
    thisFormula = paste0(thisFormula, " + X.y + X.pp")
  }
  
  thisFormula = as.formula(thisFormula)
  
  fam.pp = ifelse(bernApprox, "binomial", "poisson")
  
  startModelFitTime = proc.time()[3]
  
  mod = inla(
    thisFormula, 
    data = stackDat, offset=offset, 
    control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles), 
    family=c("gaussian", fam.pp), verbose=verbose, control.inla=control.inla, 
    control.compute=list(config=TRUE, cpo=doModAssess, dic=doModAssess, waic=doModAssess), 
    control.mode=modeControl, 
    control.fixed=controlFixed, 
    control.family=list(control.family.gaussian, list())
  )
  
  endModelFitTime = proc.time()[3]
  totalModelFitTime = endModelFitTime - startModelFitTime
  
  # get prediction covariates
  predRepelMat = getRepulsionCov(predCoords, obsCoords, 
                                 repelDist=repelDist, returnSparse=TRUE, anisFac=anisFac)
  if(!is.null(fixedRepelAmt)) {
    predOffset.pp = predRepelMat * fixedRepelAmt
  } else {
    predOffset.pp = rep(0, nrow(predCoords))
  }
  
  if(debugMode) {
    browser()
    
    testGridCoords = rbind(obsGridCoords[i,], pseudoGridCoords)
    testCoords = rbind(obsCoords[i,], pseudoCoords)
    
    predGridEast = sort(unique(predCoords[,1]))
    predGridNorth = sort(unique(predCoords[,2]))
    
    squilt(predCoords[,1], predCoords[,2], as.numeric(predOffset.pp), 
           grid=list(x=predGridEast, y=predGridNorth))
    points(obsCoords[,1], obsCoords[,2], col="red")
  }
  
  # generate samples from posterior
  startTimePosteriorSampling = proc.time()[3]
  postSamples = inla.posterior.sample(nPostSamples, mod)
  endTimePosteriorSampling = proc.time()[3]
  totalTimePosteriorSampling = endTimePosteriorSampling - startTimePosteriorSampling
  
  latentMat = sapply(postSamples, function(x) {x$latent})
  latentMatMatrix = Matrix(latentMat, sparse = FALSE)
  
  hyperparNames = names(postSamples[[1]]$hyperpar)
  nuggetVars = sapply(postSamples, function(x) {1 / x$hyperpar[which(hyperparNames == "Precision for the Gaussian observations")]})
  
  latentVarNames = rownames(postSamples[[1]]$latent)
  field.yIndices = which(grepl("field.y", latentVarNames))
  fixed.yIndices = which(grepl("X.y", latentVarNames))
  field.ppIndices = which(grepl("field.pp", latentVarNames))
  fixed.ppIndices = which(grepl("X.pp", latentVarNames))
  
  # integrate intensity to calculate point process model intercept. 
  # Approximate integral via sum over pseudosites
  if(quadratureMethod == "mesh") {
    stop("mesh quadrature not supported")
  } else {
    # fixed part
    if(!is.null(xObs)) {
      # fixedPseudoMat = matrix(X.ppDense[2:(1+nPseudoPerIter),-(1:nObs)], nrow=nPseudoPerIter) %*% matrix(latentMat[fixed.ppIndices,], ncol=nPostSamples)
      fixedPseudoMat = matrix(X.pp[2:(1+nPseudoPerIter),], nrow=nPseudoPerIter) %*% matrix(latentMat[fixed.ppIndices,], ncol=nPostSamples)
    } else {
      # fixedPseudoMat = X.ppDense[2:(1+nPseudoPerIter),] %*% latentMat[fixed.ppIndices,]
      fixedPseudoMat = 0
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
    if(getPPres) {
      fixedPredMat.pp = xPred %*% latentMat[fixed.ppIndices,]
    } else {
      fixedPredMat.pp = 0
    }
  } else {
    fixedPredMat.y = 0
    fixedPredMat.pp = 0
  }
  
  spatialPredMat.y = APred %*% latentMatMatrix[field.yIndices,]
  spatialPredMat.pp = APred %*% latentMatMatrix[field.ppIndices,]
  
  predMat.y = fixedPredMat.y + spatialPredMat.y
  if(getPPres) {
    predMat.pp = fixedPredMat.pp + spatialPredMat.pp
  } else {
    predMat.pp = 0
  }
  
  # fixed parameters (aside from repulsion) currently not supported
  # if(!is.null(offsetPred)) {
  #   predMat.y = sweep(predMat.y, 1, offsetPred.y, "+")
  # }
  
  if(!is.null(predOffset.pp) && getPPres) {
    predMat.pp = sweep(predMat.pp, 1, predOffset.pp, "+")
  }
  
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
    if(getPPres) {
      fixedObsMat.pp = xObs %*% matrix(latentMat[fixed.ppIndices,], ncol=nPostSamples)
    } else {
      fixedObsMat.pp = NULL
    }
  } else {
    fixedObsMat.y = 0
    if(getPPres) {
      fixedObsMat.pp = 0
    } else {
      fixedObsMat.pp = NULL
    }
  }
  if(getPPres) {
    # do nothing since intercept is included in fixed effects now
    # fixedObsMat.pp = sweep(fixedObsMat.pp, 2, intercept.pp, "+")
  } else {
    fixedObsMat.pp = NULL
  }
  
  # get repulsion design matrix at observations (at the time they were sampled). 
  # Add effect to observation predictions (not necessary if repulsion isn't 
  # fixed, since already included in fixedObsMat.pp)
  obsRepelMat = getRepulsionCovAtObs(obsCoords=obsCoords, obsGridCoords=obsGridCoords, repelDist=repelDist, 
                                     returnSparse=TRUE, anisFac=anisFac)
  if(!is.null(fixedRepelAmt) && getPPres) {
    offsetObs = obsRepelMat %*% fixedRepelAmt
    fixedObsMat.pp = sweep(fixedObsMat.pp, 1, offsetObs, "+")
  }
  
  spatialObsMat.y = AObs %*% latentMat[field.yIndices,]
  if(getPPres) {
    spatialObsMat.pp = AObs %*% latentMat[field.ppIndices,]
  } else {
    spatialObsMat.pp = NULL
  }
  
  obsMat.y = fixedObsMat.y + spatialObsMat.y
  if(getPPres) {
    obsMat.pp = fixedObsMat.pp + spatialObsMat.pp
  } else {
    obsMat.pp = NULL
  }
  
  # currently, fixed parameters aside from fixed repulsion not supported
  # if(!is.null(offsetEst)) {
  #   obsMat = sweep(obsMat, 1, offsetEst, "+")
  # }
  
  # now make custom predictions
  if(!is.null(customFixedI)) {
    browser()
    customFixedObsMat = xObs  %*% customFixedMat %*% latentMat[fixedIndices,]
    
    customObsMat = customFixedObsMat + spatialObsMat
    
    if(!is.null(fixedParameters$beta)) {
      stop("setting customFixedI and fixedParameters simultaneously not supported")
    }
  }
  
  # add in nugget if necessary
  
  # get cluster effect variance
  clusterVars = nuggetVars
  
  if(addNugToPredCoords) {
    predMat.yNugget = predMat.y + matrix(rnorm(length(predMat.y), sd=rep(sqrt(clusterVars), each=nrow(predMat.y))), nrow=nrow(predMat.y))
  } else {
    predMat.yNugget = NULL
  }
  obsMat.yNugget = obsMat.y + matrix(rnorm(length(obsMat.y), sd=rep(sqrt(clusterVars), each=nrow(obsMat.y))), nrow=nrow(obsMat.y))
  
  # back transform predictions
  obsMat.y = invTransform(obsMat.y)
  obsMat.yNugget = invTransform(obsMat.yNugget)
  predMat.y = invTransform(predMat.y)
  if(addNugToPredCoords) {
    predMat.yNugget = invTransform(predMat.yNugget)
  } else {
    predMat.yNugget = NULL
  }
  
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
  
  if(getPPres) {
    obs.ppEst = rowMeans(obsMat.pp)
    obs.ppSDs = apply(obsMat.pp, 1, sd)
    obs.ppLower = apply(obsMat.pp, 1, quantile, probs=(1-significanceCI)/2)
    obs.ppMedian = apply(obsMat.pp, 1, median)
    obs.ppUpper = apply(obsMat.pp, 1, quantile, probs=1-(1-significanceCI)/2)
  } else {
    obs.ppEst = NULL
    obs.ppSDs = NULL
    obs.ppLower = NULL
    obs.ppMedian = NULL
    obs.ppUpper = NULL
  }
  
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
  
  if(getPPres) {
    pred.ppEst = rowMeans(predMat.pp)
    pred.ppSDs = apply(predMat.pp, 1, sd)
    pred.ppLower = apply(predMat.pp, 1, quantile, probs=(1-significanceCI)/2)
    pred.ppMedian = apply(predMat.pp, 1, median)
    pred.ppUpper = apply(predMat.pp, 1, quantile, probs=1-(1-significanceCI)/2)
  } else {
    pred.ppEst = NULL
    pred.ppSDs = NULL
    pred.ppLower = NULL
    pred.ppMedian = NULL
    pred.ppUpper = NULL
  }
  
  if(addNugToPredCoords) {
    pred.yNuggetEst = pred.yEst
    pred.yNuggetSDs = apply(predMat.yNugget, 1, sd)
    pred.yNuggetLower = apply(predMat.yNugget, 1, quantile, probs=(1-significanceCI)/2)
    pred.yNuggetMedian = apply(predMat.yNugget, 1, median)
    pred.yNuggetUpper = apply(predMat.yNugget, 1, quantile, probs=1-(1-significanceCI)/2)
  } else {
    pred.yNuggetEst = NULL
    pred.yNuggetSDs = NULL
    pred.yNuggetLower = NULL
    pred.yNuggetMedian = NULL
    pred.yNuggetUpper = NULL
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

# function for fitting the Watson et al. model to data (with sparse covariate effects)
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
# addNugToPredCoords: whether or not to compute results regarding adding nugget 
#                     at predictive coords
# getPPres: Whether or not to compute information on point process estimation
# 
# Outputs:
# INLA model, predictions, summary statistics, input data, posterior draws, etc.
fitWatsonOld = function(obsCoords, obsValues, xObs=matrix(rep(1, length(obsValues)), ncol=1), 
                        pseudoCoords=makePseudositesRect(), xPseudo=matrix(rep(1, nrow(pseudoCoords)), nrow=1), 
                        predCoords, xPred = matrix(rep(1, nrow(predCoords)), ncol=1), 
                        control.fixed = list(prec=list(default=0), mean=list(default=0)), 
                        transform=I, invTransform=I, repelDist=10, sharedInt=FALSE, 
                        mesh=getSPDEmeshSimStudy(), prior=getSPDEprior(mesh), 
                        significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                        nPostSamples=1000, verbose=TRUE, link=1, seed=NULL, 
                        doModAssess=FALSE, customFixedI=NULL, quadratureMethod=c("pseudoSites", "mesh"), 
                        previousFit=NULL, fixedRepelAmt=NULL, addNugToPredCoords=TRUE, 
                        getPPres=TRUE, 
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
  
  # check if pseudoCoords are the same as the prediction coordinates
  if(all.equal(dim(pseudoCoords), dim(predCoords)) == TRUE) {
    # stop("must debug this case")
    pseudoArePred = all.equal(c(pseudoCoords[,1], pseudoCoords[,2]), c(predCoords[,1], predCoords[,2]))
  } else {
    pseudoArePred = FALSE
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
  xObs.pp = Matrix(xObs.pp, sparse = TRUE)
  xPseudoModInt = Matrix(xPseudoModInt, sparse = TRUE)
  
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
      thisXpseudo.pp = Matrix(thisXpseudo.pp, sparse = TRUE)
    }
    
    # combine observation and pseudo-site covariates in covariates for this iteration
    thisX.pp = rbind(matrix(xObs.pp[i,], nrow=1), 
                     thisXpseudo.pp)
    
    # get repulsion covariate as a separate design matrix to it's easier to fix
    if(i == 1) {
      thisRepelMat = Matrix(matrix(rep(0, 1+nPseudoPerIter), ncol=1), sparse=TRUE)
    } else {
      thisRepelMat = getRepulsionCov(rbind(obsCoords[i,], pseudoCoords), obsCoords[1:(i-1),], 
                                     repelDist=repelDist, returnSparse=TRUE, anisFac=anisFac)
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
  X.pp.inla = inla.as.sparse(X.pp)
  X.ppDense = as.matrix(X.pp)
  
  # browser() # separate seismic covariate from others for a separate prior
  
  # construct the observation stack 
  stack.y = inla.stack(data = list(y=cbind(ys, NA), offset=rep(0, nObs)),
                       A = list(AObs, 1),
                       effects = list(field.y=latticeInds, X.y=xObs),
                       tag = "y",
                       remove.unused=FALSE)
  
  # stack.pp = inla.stack(data = list(y=cbind(NA, rs), offset=repelOffset), 
  #                       A = list(A.pp, 1),
  #                       effects = list(field.pp=latticeInds, X.pp=X.ppDense),
  #                       tag = "pp", 
  #                       remove.unused=FALSE)
  stack.pp = inla.stack(data = list(y=cbind(NA, rs), offset=repelOffset), 
                        A = list(A.pp, 1),
                        effects = list(field.pp=latticeInds, X.pp=X.pp.inla),
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
    if(getPPres) {
      if(!sharedInt && hasInt) {
        fixedObsMat.pp = xObs[,-1] %*% matrix(latentMat[fixed.ppIndices,], ncol=nPostSamples)
      } else {
        fixedObsMat.pp = xObs %*% matrix(latentMat[fixed.ppIndices,], ncol=nPostSamples)
      }
    } else {
      fixedObsMat.pp = NULL
    }
  } else {
    fixedObsMat.y = 0
    if(getPPres) {
      fixedObsMat.pp = 0
    } else {
      fixedObsMat.pp = NULL
    }
  }
  if(getPPres) {
    fixedObsMat.pp = sweep(fixedObsMat.pp, 2, intercept.pp, "+")
  } else {
    fixedObsMat.pp = NULL
  }
  
  # get repulsion design matrix at observations (at the time they were sampled). 
  # Add effect to observation predictions (not necessary if repulsion isn't 
  # fixed, since already included in fixedObsMat.pp)
  obsRepelMat = getRepulsionCovAtObs(obsCoords=obsCoords, repelDist=repelDist, returnSparse=TRUE)
  if(!is.null(fixedRepelAmt) && getPPres) {
    offsetObs = obsRepelMat %*% fixedRepelAmt
    fixedObsMat.pp = sweep(fixedObsMat.pp, 1, offsetObs, "+")
  }
  
  spatialObsMat.y = AObs %*% latentMat[field.yIndices,]
  if(getPPres) {
    spatialObsMat.pp = AObs %*% latentMat[field.ppIndices,]
  } else {
    spatialObsMat.pp = NULL
  }
  
  obsMat.y = fixedObsMat.y + spatialObsMat.y
  if(getPPres) {
    obsMat.pp = fixedObsMat.pp + spatialObsMat.pp
  } else {
    obsMat.pp = NULL
  }
  
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
  
  if(addNugToPredCoords) {
    predMat.yNugget = predMat.y + matrix(rnorm(length(predMat.y), sd=rep(sqrt(clusterVars), each=nrow(predMat.y))), nrow=nrow(predMat.y))
  } else {
    predMat.yNugget = NULL
  }
  obsMat.yNugget = obsMat.y + matrix(rnorm(length(obsMat.y), sd=rep(sqrt(clusterVars), each=nrow(obsMat.y))), nrow=nrow(obsMat.y))
  
  # back transform predictions
  obsMat.y = invTransform(obsMat.y)
  obsMat.yNugget = invTransform(obsMat.yNugget)
  predMat.y = invTransform(predMat.y)
  if(addNugToPredCoords) {
    predMat.yNugget = invTransform(predMat.yNugget)
  } else {
    predMat.yNugget = NULL
  }
  
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
  
  if(getPPres) {
    obs.ppEst = rowMeans(obsMat.pp)
    obs.ppSDs = apply(obsMat.pp, 1, sd)
    obs.ppLower = apply(obsMat.pp, 1, quantile, probs=(1-significanceCI)/2)
    obs.ppMedian = apply(obsMat.pp, 1, median)
    obs.ppUpper = apply(obsMat.pp, 1, quantile, probs=1-(1-significanceCI)/2)
  } else {
    obs.ppEst = NULL
    obs.ppSDs = NULL
    obs.ppLower = NULL
    obs.ppMedian = NULL
    obs.ppUpper = NULL
  }
  
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
  
  if(getPPres) {
    pred.ppEst = rowMeans(predMat.pp)
    pred.ppSDs = apply(predMat.pp, 1, sd)
    pred.ppLower = apply(predMat.pp, 1, quantile, probs=(1-significanceCI)/2)
    pred.ppMedian = apply(predMat.pp, 1, median)
    pred.ppUpper = apply(predMat.pp, 1, quantile, probs=1-(1-significanceCI)/2)
  } else {
    pred.ppEst = NULL
    pred.ppSDs = NULL
    pred.ppLower = NULL
    pred.ppMedian = NULL
    pred.ppUpper = NULL
  }
  
  if(addNugToPredCoords) {
    pred.yNuggetEst = pred.yEst
    pred.yNuggetSDs = apply(predMat.yNugget, 1, sd)
    pred.yNuggetLower = apply(predMat.yNugget, 1, quantile, probs=(1-significanceCI)/2)
    pred.yNuggetMedian = apply(predMat.yNugget, 1, median)
    pred.yNuggetUpper = apply(predMat.yNugget, 1, quantile, probs=1-(1-significanceCI)/2)
  } else {
    pred.yNuggetEst = NULL
    pred.yNuggetSDs = NULL
    pred.yNuggetLower = NULL
    pred.yNuggetMedian = NULL
    pred.yNuggetUpper = NULL
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
getRepulsionCov = function(predCoords, obsCoords, repelDist=10, returnSparse=FALSE, anisFac=1) {
  if(anisFac != 1) {
    # invert the geometric transformation
    obsCoords[,1] = obsCoords[,1]*anisFac
    predCoords[,1] = predCoords[,1]*anisFac
  }
  dists = rdist(predCoords, obsCoords)
  out = matrix(apply(dists, 1, function(x) {-as.numeric(any(x < repelDist))}), ncol=1)
  
  if(returnSparse) {
    out = Matrix(out, sparse = TRUE)
  } else {
    out
  }
}

# same as getRepulsionCov, but gives the values of the repulsion covariate for 
# each observation when it was sampled
getRepulsionCovAtObs = function(obsCoords, obsGridCoords, repelDist=10, returnSparse=FALSE, anisFac=1) {
  if(anisFac != 1) {
    # invert the geometric transformation
    obsCoords[,1] = obsCoords[,1] * anisFac
    obsGridCoords[,1] = obsGridCoords[,1] * anisFac
  }
  dists = rdist(obsGridCoords, obsCoords)
  diag(dists) = Inf
  repelInd = apply(dists, 1, function(x) {
    inRangeI = which(x < repelDist)
    ifelse(identical(integer(0), inRangeI), Inf, min(inRangeI))
  })
  out = matrix(-as.numeric(repelInd < 1:nrow(obsCoords)), ncol=1)
  
  if(returnSparse) {
    out = Matrix(out, sparse = TRUE)
  } else {
    out
  }
}

# function for subsetting prediction grid to get a different regular grid of 
# pseudo sites for the Watson model for the simulation study
getPseudoCoordsSimStudy = function(maxPts=2500, predGrid=NULL, getPredGridI=TRUE) {
  
  if(is.null(predGrid)) {
    out = readSurfaceRMS("data/seisTruthReplicates/RegularizedPred_1.txt")
    predGrid = out$surfFrame
  }
  
  xs = sort(unique(predGrid[,1]))
  ys = sort(unique(predGrid[,2]))
  # xlim = c(12.5, 5008.344)
  # ylim = c(0, 15012.5)
  xlim = range(xs)
  ylim = range(ys)
  xwid = diff(xlim)
  ywid = diff(ylim)
  
  # xwid / delt + 1 = newnx
  # ywid / delt + 1 = newny
  # maxPts = newnx * newny
  # (xwid / delt + 1) * (ywid / delt + 1) = maxPts
  
  solve_delt <- function(xwid, ywid) {
    f <- function(delt) {
      (xwid / delt + 1) * (ywid / delt + 1) - maxPts
    }
    uniroot(f, lower = 1e-6, upper = max(c(xwid, ywid)))$root
  }
  delt = solve_delt(xwid=xwid, ywid=ywid)
  
  newnx = floor(xwid/delt)
  newny = floor(ywid/delt)
  
  newxs = seq(xlim[1], xlim[2], l=newnx)
  newys = seq(ylim[1], ylim[2], l=newny)
  
  outGrid = make.surface.grid(list(east=newxs, north=newys))
  
  if(getPredGridI) {
    # calculate indices of points in original predGrid closest to output pseudosite grid points
    
    if(is.null(predGrid)) {
      stop("predGrid must be specified if getPredGridI == TRUE")
    }
    
    # first get distances to east and north grid points separately
    dists = rdist(xs, newxs)
    nearEast = xs[apply(dists, 2, which.min)]
    dists = rdist(ys, newys)
    nearNorth = ys[apply(dists, 2, which.min)]
    
    if(!is.list(predGrid)) {
      colnames(predGrid) = c("east", "north")
      predGrid = as.data.frame(predGrid)
    }
    
    nearOutGrid = make.surface.grid(list(east=nearEast, north=nearNorth))
    
    # now use closest east and north grid coords to match indices in predGrid
    require(prodlim)
    predGridIs = row.match(as.data.frame(nearOutGrid), predGrid)
    
    # append info to outGrid
    outGrid = cbind(outGrid, predGridIs=predGridIs)
  }
  
  outGrid
}

# head(modelFitCombs[(modelFitCombs$fitModFunI == 4) & (modelFitCombs$n == 250) & (modelFitCombs$repelAreaProp == 0.001) & (modelFitCombs$propVarCase == "diggle"),])
testPseudoConvergence = function(i=31626, nPseudos=c(500, 750, 1000, 1500, 2000, 2500, 3500, 5000), 
                                 significance=c(.8, .95), 
                                 adaptScen=c("batch", "adaptPref", "adaptVar"), 
                                 regenData=FALSE, verbose=FALSE, doPlot=TRUE, anisFac=3) {
  startT = proc.time()[3]
  adaptScen = match.arg(adaptScen)
  
  pseudoFile=paste0("savedOutput/testPseudo/testPseudo_", i, ".RData")
  
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
  
  if(fitModFunI != 4) {
    stop("fitModFun isn't Watson model")
  }
  
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
  
  
  # subset well data based on n for this run
  if("wellDat" %in% names(wellDat)) {
    wellDat = wellDat$wellDat[1:n,]
  } else {
    wellDat = wellDat[1:n,]
  }
  wellDat[,1] = wellDat[,1]/anisFac
  
  # interpolate truth to well points
  truthWells = bilinearInterp(wellDat[,1:2], truth, transform=logit, invTransform=expit)
  
  # make sure prior on preferentiality is centered at truth
  
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
  
  predGrid = cbind(east=seismicDat$east, north=seismicDat$north)
  mesh = getSPDEmeshSimStudy(anisFac=anisFac)
  
  # loop through pseudosite resolutions
  
  if(!file.exists(pseudoFile) || regenData) {
    fullPredMat = c()
    fixedList = list()
    parList = list()
    allTs = numeric(length(nPseudos))
    for(j in 1:length(nPseudos)) {
      nPseudo = nPseudos[j]
      
      print(paste0("fitting Watson model for res ", nPseudo, " and j=", j, "/", length(nPseudos)))
      
      # Fit model and calculate scores if need be
      
      if(fitModFunI == 4) {
        repDist = repAreaToDist(repelAreaProp)
        
        inputList = list(wellDat, seismicDat, repelDist=repDist, nPseudo=nPseudo, prefMean=prefPar, mesh=mesh, anisFac=anisFac, getPPres=TRUE)
      }
      
      
      out = do.call("fitModFun", inputList)
      
      predMat = out$predMat.y # doesn't include nugget
      predAggMat = out$pred.yAggMat # doesn't include nugget
      obsMat = out$obsMat.y # doesn't include nugget
      
      fixedEffectSummary = out$fixedEffectSummary
      parameterSummaryTable = out$parameterSummaryTable
      
      # concatenate results
      fullPredMat = cbind(fullPredMat, rowMeans(predMat))
      fixedList = c(fixedList, list(fixedEffectSummary))
      parList = c(parList, list(parameterSummaryTable))
      allTs[j] = out$timings$totalTime
      
      if(doPlot) {
        
        if("logitProbsNoRep" %in% names(wellDat)) {
          logitProbsNoRepWells = wellDat$logitProbsNoRep
        }
        
        
        
        preds = rowMeans(predMat)
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
        summary(lm(logitProbsNoRepWells ~ I(logit(wellDat$volFrac)) + I(logit(pSeismic))))
        summary(lm(logitProbsNoRep ~ I(logit(gTruth)) + I(logit(gSeismic))))
        
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
        
        pdf(file=paste0("figures/testPseudo/testPreds_simStudy_", adaptScen, "_par", 
                        sampleParI, "_rep", repI, "_nPseudo", nPseudo, "_nPseudo", nPseudo, ".pdf"), width=8, height=5)
        par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
        squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
               xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
               zlim=range(c(preds, gTruth)), asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
        points(pEast, pNorth, cex=.5)
        
        squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
               xlab="Easting", ylab="Northing", main="Seismic Estimate", 
               asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
        points(pEast, pNorth, cex=.5)
        
        splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
              resetGraphics=FALSE, 
              zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
              asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
        
        squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
               colScale=seqCols, zlim=range(c(preds, gTruth)), 
               xlab="Easting", ylab="Northing", main="Estimate", 
               asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
        points(pEast, pNorth, cex=.5)
        dev.off()
        
        pdf(file=paste0("figures/testPseudo/testQuants_simStudy_", adaptScen, "_par", 
                        sampleParI, "_rep", repI, "_nPseudo", nPseudo, "_nPseudo", nPseudo, ".pdf"), width=8, height=5)
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
        
        pdf(file=paste0("figures/testPseudo/testLogitSDs_simStudy_", adaptScen, "_par", 
                        sampleParI, "_rep", repI, "_nPseudo", nPseudo, "_nPseudo", nPseudo, ".pdf"), width=8, height=5)
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
          probsNoRep = expit(logitProbsNoRep)
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
          
          pdf(file=paste0("figures/testPseudo/testLambdas_simStudy_", adaptScen, "_par", 
                          sampleParI, "_rep", repI, "_nPseudo", nPseudo, ".pdf"), width=16, height=7.5)
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
          
          
          
          
          pdf(file=paste0("figures/testPseudo/testLambdas2_simStudy_", adaptScen, "_par", 
                          sampleParI, "_rep", repI, "_nPseudo", nPseudo, ".pdf"), width=16, height=7.5)
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
        
      }
    }
    
    save(fullPredMat, fixedList, parList, allTs, nPseudos, file=pseudoFile)
  } else {
    out = load(pseudoFile)
  }
  
  # plot results
  browser()
  
  
  
  invisible(NULL)
}








