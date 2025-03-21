# This script fits the successive (one-at-a-time) SPDE model to data and generates predictions

# function for fitting the SPDE model to data
# 
# Inputs:
# ptDat: data.frame with columns: east, north
# gridDat: data.frame with columns: east, north, [additional columns with covariate info]
# transform: how to transform obsValues prior to modeling
# invTransform: inverse of transform. Used to backtransform predictions prior to aggregation
# mesh: SPDE mesh
# prior: SPDE prior
# significanceCI: the credible level of the CIs (e.g., .8 means 80% CI)
# int.strategy: inla integration strategy
# strategy: inla strategy
# pProcMethod: Method for fitting the one of "kde" (default) or "inlabru"
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
# minN: minimum number of points required to fit the point process density
# 
# Outputs:
# INLA model, predictions, summary statistics, input data, posterior draws, etc.
fitSuccessiveSPDE = function(ptDat, gridDat, 
                             transform=I, invTransform=I, 
                             mesh=getSPDEmesh(obsCoords), prior=getSPDEprior(mesh), 
                             significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                             pProcMethod=c("kde", "inlabru"), 
                             nPostSamples=1000, verbose=TRUE, link=1, seed=NULL, 
                             family=c("normal", "binomial", "betabinomial"), 
                             doModAssess=FALSE, 
                             previousFit=NULL, improperCovariatePrior=TRUE, 
                             fixedParameters=NULL, experimentalMode=FALSE, 
                             minN=2) {
  
  # select point process density fitting method
  pProcMethod = match.arg(pProcMethod)
  
  # set up data at point data locations
  obsPts = ptDat[,1:2]
  obsVals = ptDat[,3]
  xObs = matrix(nrow=nrow(obsPts), ncol=ncol(gridDat)-2+2) # (-2 for east, north, +2 for intercept, successive ests)
  xObs[,1] = 1 # add in intercept
  xObs[,2] = NA # add in w (successive spatial estimates that start at 0)
  xObs[1:minN, 2] = 0
  
  # interpolate gridded covarates onto the point locations via bilinear interpolation
  gridColNames = names(gridDat)
  for(i in 3:ncol(gridDat)) {
    interpVals = bilinearInterp(obsPts, gridDat[,i])
    xObs[,i] = interpVals
  }
  col.names(xObs) = c("int", "w", gridColNames[-(1:2)])
  
  # set up info needed for final predictions (fill in values of w later)
  predPts = gridDat[,1:2]
  xPred = cbind(1, NA, gridDat[,-(1:2)])
  
  # Successively fit the pointp process density model, each time adding one more 
  # datapoint.
  
  print(paste0("starting to fit successive SPDE model..."))
  for(i in startN:nrow(obsCoords)) {
    
    # use only the first i well points for model fitting
    thisObsCoords = matrix(obsCoords[1:i,], ncol=2)
    thisObsValues = obsValues[1:i]
    thisXobs = matrix(xObs[1:i,], ncol=ncol(xObs))
    
    
    if(i != nrow(obsCoords)) {
      # use only the next well point for prediction unless it's the last iteration
      thisPredCoords = matrix(obsCoords[i+1,], nrow=1)
      thisXpred = matrix(xObs[i+1,], nrow=1)
    } else {
      # for the last fit, generate predictions over whole grid
      thisPredCoords = predPts
      thisXpred = xPred
      
      # add in w based on 
    }
    
    fit = fitSPDE(obsCoords=thisObsCoords, obsValues=thisObsValues, xObs=thisXobs, 
                  predCoords=thisPredCoords, xPred=thisXpred, 
                  transform=transform, invTransform=invTransform, 
                  mesh=mesh, prior=prior, 
                  significanceCI=significanceCI, int.strategy=int.strategy, strategy=strategy, 
                  nPostSamples=nPostSamples, verbose=verbose, link=link, seed=seed, 
                  family=family, doModAssess=doModAssess, 
                  previousFit=previousFit, improperCovariatePrior=improperCovariatePrior, 
                  fixedParameters=fixedParameters, experimentalMode=experimentalMode)
    
    if(i != nrow(obsCoords)) {
      # If this isn't the final fit:
      
      # update previousFit for faster optimization/model fitting in next iteration
      previousFit = fit$mod
      
      # update w (covariate representing information from previous point data)
      browser()
      nextW = mean(fit$spatialPredMat)
      xObs[i+1,2] = nextW
    } else {
      # If this is the 
    }
    
    
    
    
  }
  
}









