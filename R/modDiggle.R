# This script fits Diggle model to data and generates predictions

# function for fitting the Diggle model to data
# 
# Inputs:
# wellDat: data.frame with columns: east, north
# seismicDat: data.frame with columns: east, north, [additional columns with covariate info]
# studyArea: Polygons that delimits the area where observations are made
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


fitDigglesimDat = function(wellDat, seismicDat, studyArea,
                           predGrid=cbind(east=seismicDat$east, north=seismicDat$north),
                           control.fixed = list(prec=list(default=0, X2=1/.5^2, X3=1), mean=list(default=0, X2=1)), ##NBNB: Review
                           transform=logit, invTransform=expit, 
                           mesh=getSPDEmeshSimStudy(), prior=getSPDEprior(mesh), 
                           #addKDE=FALSE, esthFromSeismic=TRUE, kde.args=NULL, pProcMethod=c("kde", "inlabru"), 
                           significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                           nPostSamples=1000, verbose=TRUE, seed=123, 
                           family="normal", doModAssess=FALSE, previousFit=NULL, 
                           fixedParameters=NULL, experimentalMode=FALSE) {
  
  #mesh = inla.mesh.2d(boundary=studyArea,loc.domain = testData$wellDat[,c("east","north")],max.edge=c(1000,10000),cutoff=100,offset=c(10,-.1))
  #Not sure about kde stuff
  
  predPts = st_as_sf(seismicDat[,1:2], coords = c("east", "north"))
  #xPred = data.frame(X =  transform(seismicDat$seismicEst))
  
  # interpolate seismic data to the well points
  wellSeismicEsts = bilinearInterp(wellDat[,1:2], seismicDat, transform = transform, invTransform = invTransform)
  
  # construct well data covariates
  xObs = data.frame(X = transform(wellSeismicEsts))
  
  # set observations
  obsValues = wellDat$volFrac
  obsCoords = cbind(wellDat$east, wellDat$north)
  
  fitDiggle(obsCoords=obsCoords, obsValues=obsValues, xObs=xObs, studyArea=studyArea,
            covs = list(X= seismicDat),
            predCoords=predPts, control.fixed = control.fixed,#xPred=xPred, 
            transform=transform, invTransform=invTransform, 
            mesh=mesh, prior=prior, 
            significanceCI=significanceCI, int.strategy=int.strategy, strategy=strategy, 
            nPostSamples=nPostSamples, verbose=verbose, link=link, seed=seed, 
            family=family, doModAssess=doModAssess, previousFit=previousFit, 
            fixedParameters=fixedParameters, experimentalMode=experimentalMode)
  
}


#Document

fitDiggle = function(obsCoords, obsValues, xObs=matrix(rep(1, length(obsValues)), ncol=1), 
                     studyArea=studyArea, ## Do it inside, not user input
                     covs= NULL,
                     predCoords, #xPred = matrix(rep(1, nrow(predCoords)), ncol=1), 
                     control.fixed = list(prec=list(default=0), mean=list(default=0)), 
                     transform=I, invTransform=I, 
                     mesh=getSPDEmesh(obsCoords), prior=getSPDEprior(mesh), 
                     significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                     nPostSamples=1000, verbose=TRUE, link=1, seed=NULL,
                     family=c("normal", "binomial", "betabinomial"), 
                     doModAssess=FALSE, customFixedI = NULL,
                     previousFit=NULL,
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
  wellDat$X = as.numeric(xObs[,1]) ##NBNB
  ys = as.numeric(transform(obsValues))
  # Making X terra
  X_terra = terra::rast(covs$X, type="xyz")
  
  
  #prior = inla.spde2.pcmatern(mesh, prior.range=c(1000, 0.5), prior.sigma = c(0.5, 0.5))
  #wellDat$volFrac = as.numeric(wellDat$volFrac) #Review this part please!
  wellDat$y = ys
  well_data_sf =st_as_sf(wellDat, coords = c("east", "north"), crs = crs(studyArea))  
  cmp <- ~ field_pp(geometry, copy = "field_y", fixed =F) + Intercept_y(1) + Intercept_pp(1) + X(X_terra, model="linear") +
    field_y(geometry, model=prior)
  cmp_alt <- ~ field_pp(geometry, model = prior) + Intercept_y(1) +
    Intercept_pp(1) + X(X_terra, model="linear") +
    field_y(geometry, model = prior)
  
  mesh$crs <- st_crs(studyArea)$wkt

  ## Check for covariate extent
  lik1_vect <- vect(mesh$loc[,1:2])
  ext_rast <- ext(X_terra)
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
                  formula = geometry ~  Intercept_pp + X + field_pp,
                  data = well_data_sf,
                  domain = list(geometry = mesh),
                  samplers = studyArea
  )
  
  lik2 <- bru_obs(family,
                  formula = y ~  Intercept_y + X + field_y,
                  data = well_data_sf,
                  domain = list(geometry = mesh),
                  control.family = control.family
  )
  
  allQuantiles = c(0.5, (1-significanceCI) / 2, 1 - (1-significanceCI) / 2)
  controlFixed=list(quantiles=allQuantiles)
  
  modeControl = inla.set.control.mode.default()

    mod <- bru(cmp, lik1, lik2,
             options = list(
               bru_verbose=3,
               control.inla = list(strategy = strategy, int.strategy = int.strategy, verbose=TRUE),
               control.mode = modeControl,
               control.fixed = controlFixed
             )
  )
  
  interceptSummary = mod$summary.fixed[2,1:5]
  fixedEffectSummary = mod$summary.fixed[3,1:5]
  rangeSummary = mod$summary.hyperpar[2,1:5]  
  
  pred <- predict(
    mod,
    predCoords,
    ~ data.frame(
      #lambda = exp(eco_intercept + eco_cov + w1),
      log_lambda = Intercept_pp + field_pp,
      volFrac = as.numeric(invTransform(Intercept_y  + X + field_y)),
      AggObs = mean(Intercept_y + X + field_y)
    ))
  
  predEst = pred$volFrac$mean
  predSDs = pred$volFrac$sd
  predLower = pred$volFrac$q0.025
  predMedian = pred$volFrac$median
  predUpper = pred$volFrac$q0.975
  
  predAggEst = pred$AggObs$mean[1]
  predAggSDs = pred$AggObs$sd[1]
  predAggLower = pred$AggObs$q0.025[1]
  predAggMedian = pred$AggObs$median[1]
  predAggUpper = pred$AggObs$q0.975[1]
  
  
  obsCoords_sf = st_as_sf(as.data.frame(obsCoords), coords = c(1, 2), crs = "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs")
  
  
  obs_pred <- predict(
    mod,
    obsCoords_sf,
    ~ data.frame(
      #lambda = exp(eco_intercept + eco_cov + w1),
      fixedObs = Intercept_y + X,
      spatialObs = field_y,
      Obs = Intercept_y + X + field_y,
      AggObs = mean(Intercept_y + X + field_y)
      #volFrac = Intercept_y  + X + field_y
    ))
  
  
  obsEst = obs_pred$Obs$mean
  obsSDs = obs_pred$Obs$sd
  obsLower = obs_pred$Obs$q0.025
  obsMedian = obs_pred$Obs$median
  obsUpper = obs_pred$Obs$q0.975
  
  postSamples = inla.posterior.sample(nPostSamples, mod) #Use this then
  
  if(family == "normal") {
    hyperparNames = names(postSamples[[1]]$hyperpar)
    nuggetVars = sapply(postSamples, function(x) {1 / x$hyperpar[which(hyperparNames == "Precision for the Gaussian observations[2]")]})
  }
  
  hyperMat = sapply(postSamples, function(x) {x$hyperpar})
  hyperNames = row.names(hyperMat)
  if(family == "normal") {
    clusterVarI = which(hyperNames == "Precision for the Gaussian observations[2]")
    spatialRangeI = which(hyperNames == "Range for field_y" )
    spatialSDI = which(hyperNames == "Stdev for field_y")
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
  timings = data.frame(totalTime = totalTime)
  list(mod=mod, 
       obsCoords=obsCoords, xObs=xObs, obsValues=obsValues, #predCoords=predCoords,
       #xPred=xPred, 
       obsEst=obsEst, obsSDs=obsSDs, obsLower=obsLower, obsMedian=obsMedian, obsUpper=obsUpper, 
       predEst=predEst, predSDs=predSDs, predLower=predLower, predMedian=predMedian, predUpper=predUpper, #predAggMat=predAggMat, 
       predAggEst=predAggEst, predAggSDs=predAggSDs, predAggLower=predAggLower, predAggMedian=predAggMedian, predAggUpper=predAggUpper, 
       mesh=mesh, prior=prior, #stack=stack.full, 
       interceptSummary=interceptSummary, fixedEffectSummary=fixedEffectSummary, rangeSummary=rangeSummary, 
       sdSummary=sdSummary, varSummary=varSummary, 
       parameterSummaryTable=parameterSummaryTable, 
       #fixedEffectDraws=latentMat[fixedIndices,], 
       #spatialPredMat=spatialPredMat, fixedPredMat=fixedPredMat, 
       #spatialObsMat=spatialObsMat, fixedObsMat=fixedObsMat, 
       #obsMat=obsMat, obsMatNugget=obsMatNugget, predMat=predMat, predMatNugget=predMatNugget, 
       hyperMat=hyperMat, timings=timings#, sigmaEpsilonDraws=sqrt(clusterVars)
  )
  
  
}
