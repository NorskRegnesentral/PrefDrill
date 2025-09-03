# This script fits Diggle model to data and generates predictions

# function for fitting the Diggle model to data
# 
# Inputs:
# wellDat: data.frame with columns: east, north
# seismicDat: data.frame with columns: east, north, [additional columns with covariate info]
# predGrid: grid over which to make predictions
# transform: how to transform obsValues prior to modeling
# invTransform: inverse of transform. Used to backtransform predictions prior to aggregation
# mesh: SPDE mesh
# prior: SPDE prior
# significanceCI: the credible level of the CIs (e.g., .8 means 80% CI) [Not useful so far]
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
# INLA model, predictions, summary statistics, input data, (posterior draws), etc.

fitDigglesimDat = function(wellDat, seismicDat,
                           predGrid=cbind(east=seismicDat$east, north=seismicDat$north),
                           control.fixed = list(prec=list(default=0, X_y=1/.5^2), mean=list(default=0, X_y=1)), ##NBNB: Review
                           transform=logit, invTransform=expit, 
                           mesh=getSPDEmeshSimStudy(), prior=getSPDEprior(mesh), 
                           significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                           nPostSamples=1000, verbose=FALSE, seed=NULL, 
                           family="normal", doModAssess=FALSE, previousFit=NULL, 
                           addNugToPredCoords=FALSE, prefMean=0, 
                           fixedParameters=NULL, experimentalMode=FALSE) {
  
  # construct prediction points
  predPts = st_as_sf(seismicDat[,1:2], coords = c("east", "north"))
  
  # interpolate seismic data to the well points
  wellSeismicEsts = bilinearInterp(wellDat[,1:2], seismicDat,
                                   transform = transform, invTransform = invTransform)

  # construct well data covariates
  xObs = data.frame(X = transform(wellSeismicEsts))
  
  # set observations
  obsValues = wellDat$volFrac
  obsCoords = cbind(wellDat$east, wellDat$north)
  
  seismicDat$seismicEst = transform(seismicDat$seismicEst)
  
  fitDiggle(obsCoords=obsCoords, obsValues=obsValues, xObs=xObs,
            covs = list(X= seismicDat),
            predCoords=predPts, 
            control.fixed = control.fixed,
            transform=transform, invTransform=invTransform, prefMean=prefMean, 
            mesh=mesh, prior=prior, significanceCI=significanceCI, prefMean=prefMean, 
            int.strategy=int.strategy, strategy=strategy, nPostSamples=nPostSamples, 
            verbose=verbose, link=link, seed=seed, addNugToPredCoords=addNugToPredCoords, 
            family=family, doModAssess=doModAssess, previousFit=previousFit, 
            fixedParameters=fixedParameters, experimentalMode=experimentalMode)
  
}


#Document!

fitDiggle = function(obsCoords, obsValues, xObs=matrix(rep(1, length(obsValues)), ncol=1), 
                     covs,predCoords,
                     control.fixed = list(prec=list(default=0), mean=list(default=0)), 
                     transform=I, invTransform=I, 
                     mesh=getSPDEmesh(obsCoords), prior=getSPDEprior(mesh), 
                     significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                     nPostSamples=1000, verbose=TRUE, link=1, seed=NULL,
                     family=c("normal", "binomial", "betabinomial"), 
                     doModAssess=FALSE, customFixedI = NULL, prefMean=0, 
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
  control.family = list(hyper = list(prec = list(prior="loggamma", param=c(1000,10))))

  if(!is.null(fixedParameters$familyPrec)) {
    # fix the family precision parameter on INLA's latent scale
    control.family = list(initial=log(fixedParameters$familyPrec), fixed=TRUE)
  }
  
  ## inlabru code
  prefPrior = list(prior="gaussian", param=c(prefMean, 1/4))
  ys = as.numeric(transform(obsValues))
  # Making X terra
  X_terra = terra::rast(covs$X, type="xyz")
  
  wellDat = data.frame(X=as.numeric(xObs[,1]),
                       y = ys,
                       volFrac = obsValues,
                       east=obsCoords[,1],
                       north=obsCoords[,2]
                       )

  well_data_sf =st_as_sf(wellDat, coords = c("east", "north"))  
  
  cmp <- ~ field_pp(geometry, copy = "field_y", fixed =F, hyper = list(beta=prefPrior)) + Intercept_y(1) + Intercept_pp(1) + X_pp(X_terra, model="linear") +
    field_y(geometry, model=prior) + X_y(X, model="linear")
  
  ## Construct study area
  ext_rast <- ext(X_terra)
  studyArea <- sf::st_as_sf(as.polygons(ext_rast))
  
  ## Check for covariate extent
  lik1_vect <- vect(mesh$loc[,1:2])
  ext_pts <- ext(lik1_vect)
  combined_ext <- ext(
    min(ext_rast$xmin, ext_pts$xmin),
    max(ext_rast$xmax, ext_pts$xmax),
    min(ext_rast$ymin, ext_pts$ymin),
    max(ext_rast$ymax, ext_pts$ymax)
  )
  expanded_rast <- extend(X_terra, combined_ext)
  while(any(is.na(values(expanded_rast)))){
    w <- matrix(1, 3, 3)
    expanded_rast <- focal(expanded_rast, w = w, fun = mean, na.policy = "only", na.rm = TRUE)
    X_terra <- expanded_rast
  }
  
  lik1 <- bru_obs("cp",
                  formula = geometry ~  Intercept_pp + X_pp + field_pp,
                  data = well_data_sf[,c("volFrac","geometry")],
                  domain = list(geometry = mesh),
                  samplers = studyArea
  )
  
  lik2 <- bru_obs(family,
                  formula = y ~  Intercept_y + X_y + field_y,
                  data = well_data_sf,
                  control.family = control.family
  )
  
  allQuantiles = c(0.5, (1-significanceCI) / 2, 1 - (1-significanceCI) / 2)

  modeControl = inla.set.control.mode.default()

  endTimeDefineModel = proc.time()[3]
  totalTimeDefineModel = endTimeDefineModel - startTimeDefineModel
  
  startModelFitTime = proc.time()[3]
  
  mod <- bru(cmp, lik1,lik2,
             options = list(
               bru_verbose=0,bru_max_iter=1,
               control.inla = list(strategy = strategy, int.strategy = int.strategy)
               )
             )
  endModelFitTime = proc.time()[3]
  totalModelFitTime = endModelFitTime - startModelFitTime
  
  
  startTimePosteriorSampling = proc.time()[3]
  postSamples = inla.posterior.sample(nPostSamples, mod)
  endTimePosteriorSampling = proc.time()[3]
  totalTimePosteriorSampling = endTimePosteriorSampling - startTimePosteriorSampling
  
  latentMat = sapply(postSamples, function(x) {x$latent})
  
  if(family == "normal") {
    hyperparNames = names(postSamples[[1]]$hyperpar)
    nuggetVars = sapply(postSamples, function(x) {1 / x$hyperpar[which(hyperparNames == "Precision for the Gaussian observations[2]")]})
  }
  
  latentVarNames = rownames(postSamples[[1]]$latent)
  field_y_Indices = which(grepl("field_y", latentVarNames))
  field_pp_Indices = which(grepl("field_pp", latentVarNames))
  fixedIndices = which(grepl("X_y", latentVarNames))
  fixedInt_y_Indices = which(grepl("Intercept_y", latentVarNames))
  fixedInt_pp_Indices = which(grepl("Intercept_pp", latentVarNames))
  
  ## Spatial GRF predictions
  APred = inla.spde.make.A(mesh, loc = predCoords)
  AObs = inla.spde.make.A(mesh,loc=well_data_sf)
  spatial_y_PredMat = APred %*% latentMat[field_y_Indices,]
  spatial_pp_PredMat = APred %*% latentMat[field_pp_Indices,]
  spatial_y_ObsMat = AObs %*% latentMat[field_y_Indices,]

  fixed_y_pred = latentMat[fixedInt_y_Indices,]
  fixed_pp_pred = latentMat[fixedInt_pp_Indices,]
  
  X_pred_covs = extract(X_terra,predCoords)[,2]
  X_pred_cov_mat = matrix(rep(X_pred_covs,ncol(latentMat)), 
                          ncol = ncol(latentMat))
  X_obs_mat = matrix(rep(well_data_sf$X,ncol(latentMat)), 
                     ncol = ncol(latentMat))
  
  X_samps = matrix(latentMat[fixedIndices,], nrow = 1, ncol = ncol(latentMat))
  fixedPredMat = sweep(X_pred_cov_mat, MARGIN = 2, STATS = X_samps, FUN = "*")
  fixedObsMat = sweep(X_obs_mat, MARGIN = 2, STATS = X_samps, FUN = "*")
  
  Predbase_int_mat = matrix(1,nrow=nrow(predCoords),ncol=nPostSamples)
  IntPred_y_samps = sweep(Predbase_int_mat, MARGIN = 2, STATS = latentMat[fixedInt_y_Indices,], FUN = "*")
  IntPred_pp_samps = sweep(Predbase_int_mat, MARGIN = 2, STATS = latentMat[fixedInt_pp_Indices,], FUN = "*")
  
  Obsbase_int_mat = matrix(1,nrow=nrow(well_data_sf),ncol=nPostSamples)
  IntObs_y_samps = sweep(Obsbase_int_mat, MARGIN = 2, STATS = latentMat[fixedInt_y_Indices,], FUN = "*")
 
  predMat = IntPred_y_samps + fixedPredMat + spatial_y_PredMat
  obsMat = IntObs_y_samps + fixedObsMat + spatial_y_ObsMat
  
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
  
  obsMat = invTransform(obsMat)
  obsMatNugget = invTransform(obsMatNugget)
  predMat = invTransform(predMat)
  if(addNugToPredCoords) {
    predMatNugget = invTransform(predMatNugget)
  } else {
    predMatNugget = NULL
  }
  
  # summary statistic
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
  
  # aggregate summary statistics
  predAggMat = colMeans(predMat)
  predAggEst = mean(predAggMat)
  predAggSDs = sd(predAggMat)
  predAggLower = quantile(predAggMat, probs=(1-significanceCI)/2)
  predAggMedian = quantile(predAggMat, probs=.5)
  predAggUpper = quantile(predAggMat, probs=1-(1-significanceCI)/2)
  
  interceptSummary=mod$summary.fixed[3, 1:5]
  fixedEffectSummary=mod$summary.fixed[grepl("_y", row.names(mod$summary.fixed)), 1:5]
  
  rangeSummary=mod$summary.hyperpar[which(hyperparNames == "Range for field_y"), 1:5]
  spatialSDSummary = mod$summary.hyperpar[which(hyperparNames == "Stdev for field_y"), 1:5]
  prefParSummary = mod$summary.hyperpar[which(hyperparNames == "Beta for field_pp"), 1:5]
  
  # posterior hyperparameter samples
  hyperMat = sapply(postSamples, function(x) {x$hyperpar})
  if(family == "normal") {
    clusterVarI = which(hyperparNames == "Precision for the Gaussian observations[2]")
    spatialRangeI = which(hyperparNames == "Range for field_y" )
    spatialSDI = which(hyperparNames == "Stdev for field_y")
    prefParI = which(hyperparNames == "Beta for field_pp")
    
    if(!is.matrix(hyperMat)) {
      mat = NULL
    } else {
      mat = apply(hyperMat, 2, function(x) {c(totalVar=x[spatialSDI]^2+1/x[clusterVarI], spatialVar=x[spatialSDI]^2, errorVar=1/x[clusterVarI], 
                                              totalSD=sqrt(x[spatialSDI]^2+1/x[clusterVarI]), spatialSD=x[spatialSDI], errorSD=sqrt(1/x[clusterVarI]), 
                                              spatialRange=x[spatialRangeI], preferentialParameter=x[prefParI])})
    }
  } else {
    stop("family not supported")
  }

  if(family == "normal")
    hyperNames = c("totalVar", "spatialVar", "errorVar", "totalSD", "spatialSD", "errorSD", "spatialRange", "prefPar")
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
      prefParSummary = parameterSummaryTable[summaryHyperNames == "prefPar",]
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
       obsCoords=obsCoords, xObs=xObs, obsValues=obsValues, predCoords=predCoords,
       obsEst=obsEst, obsSDs=obsSDs, obsLower=obsLower, obsMedian=obsMedian, obsUpper=obsUpper, 
       predEst=predEst, predSDs=predSDs, predLower=predLower, predMedian=predMedian, predUpper=predUpper, predAggMat=predAggMat, 
       predAggEst=predAggEst, predAggSDs=predAggSDs, predAggLower=predAggLower, predAggMedian=predAggMedian, predAggUpper=predAggUpper, 
       mesh=mesh, prior=prior,
       interceptSummary=interceptSummary, fixedEffectSummary=fixedEffectSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, 
       parameterSummaryTable=parameterSummaryTable, 
       fixedEffectDraws=latentMat[fixedIndices,], 
       spatialPredMat=spatial_y_PredMat, fixedPredMat=fixedPredMat, 
       spatialObsMat=spatial_y_ObsMat, fixedObsMat=fixedObsMat, 
       obsMat=obsMat, obsMatNugget=obsMatNugget, predMat=predMat, predMatNugget=predMatNugget, 
       hyperMat=hyperMat, timings=timings, sigmaEpsilonDraws=sqrt(clusterVars))

}