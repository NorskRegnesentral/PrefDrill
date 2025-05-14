# This script fits the successive (one-at-a-time) SPDE model to data and generates predictions

# function for fitting the SPDE model to data
# 
# Inputs:
# ptDat: data.frame with columns: east, north
# gridDat: data.frame with columns: east, north, [additional columns with covariate info]
# transform: how to transform obsVals prior to modeling
# invTransform: inverse of transform. Used to backtransform predictions prior to aggregation
# mesh: SPDE mesh
# prior: SPDE prior
# significanceCI: the credible level of the CIs (e.g., .8 means 80% CI)
# int.strategy: inla integration strategy
# strategy: inla strategy
# pProcMethod: Method for fitting the one of "kde" (default) or "inlabru"
# priorW: prior mean for w, on the transformed scale
# useLastPtPProc: if successiveMethod=="pproc", uses i-th location to estimate 
#                 w_i in addition to locations 1, ..., i-1.
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
# kde.args: arguments to getDensityDueToWells
# esthFromSeismic: whether to estimate h, kde's bandwidth parameter, from 
#                  seismic/grid data
# ...: additional arguments passed to fitSPDE
# 
# Outputs:
# allWs matrix of estimates of w ovver grid at each iterations, obsWs, INLA 
# model, predictions, summary statistics, input data, posterior draws, etc.
fitSuccessiveSPDEsimDat = function(ptDat, gridDat, 
                                   transform=logit, invTransform=expit, 
                                   mesh=getSPDEmeshSimStudy(), prior=getSPDEprior(mesh), 
                                   significanceCI=.8, int.strategy="ccd", strategy="simplified.laplace", 
                                   successiveMethod=c("pproc", "spde"), pProcMethod=c("kde", "inlabru"), includeInt=FALSE, 
                                   priorW=ifelse(includeInt, 0, 1), useLastPtPProc=TRUE, 
                                   nPostSamples=1000, verbose=FALSE, link=1, seed=NULL, 
                                   family=c("normal", "binomial", "betabinomial"), 
                                   doModAssess=FALSE, storeSuccessiveWs=TRUE, 
                                   previousFit=NULL, improperCovariatePrior=TRUE, 
                                   fixedParameters=NULL, experimentalMode=FALSE, 
                                   minN=4, kde.args=NULL, esthFromSeismic=TRUE, ...) {
  
  # select default methods
  successiveMethod = match.arg(successiveMethod)
  pProcMethod = match.arg(pProcMethod)
  family = match.arg(family)
  
  # check if seismic data is included
  seismicIncluded = ncol(gridDat) >= 3
  if(!seismicIncluded) {
    stop("Not yet tested without seismic data")
  }
  
  # set up data at point data locations
  obsPts = as.matrix(ptDat[,1:2])
  obsVals = ptDat[,3]
  nCovs = ifelse(includeInt, ncol(gridDat)-2+2, ncol(gridDat)-2+1)
  xObs = matrix(nrow=nrow(obsPts), ncol=nCovs) # (-2 for east, north, +2 for intercept, successive ests)
  gridColNames = names(gridDat)
  if(includeInt) {
    wCol = 2
    seismicCols = 3:ncol(gridDat)
    xObs[,1] = 1 # add in intercept
    xObs[,wCol] = priorW # add in w (successive spatial estimates that start at 0)
    colnames(xObs) = c("int", "w", gridColNames[-(1:2)])
    customFixedI = c(rep(1, 2), rep(0, ncol(gridDat)-2))
  } else {
    wCol = 1
    seismicCols = (3:ncol(gridDat)) - 1
    xObs[,wCol] = priorW # add in w (successive spatial estimates that start at 0)
    colnames(xObs) = c("w", gridColNames[-(1:2)])
    customFixedI = c(rep(1, 1), rep(0, ncol(gridDat)-2))
  }
  
  # interpolate gridded covariates onto the point locations via bilinear interpolation
  colDiff = ifelse(includeInt, 0, 1)
  for(i in 3:ncol(gridDat)) {
    interpVals = bilinearInterp(obsPts, cbind(gridDat[,1:2], gridDat[,i]))
    xObs[,i-colDiff] = interpVals
  }
  
  # set up info needed for final predictions (fill in values of w later)
  predPts = gridDat[,1:2]
  xPred = cbind(priorW, gridDat[,-(1:2)])
  if(includeInt) {
    xPred = cbind(1, xPred)
  }
  
  # if need be, estimate bandwidth based on seismic data
  if(esthFromSeismic && successiveMethod == "pproc") {
    # estimate bandwidth via a variogram estimator on the seismic data
    require(gstat)
    require(sf)
    # transformedGridEsts = transform(gridDat[,3])
    varioMod = vgm(model="Gau")
    varioDat = data.frame(gridEsts=gridDat[,3], x=gridDat$east, y=gridDat$north)
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
  
  # Successively fit the point process density model, each time adding one more 
  # datapoint.
  
  
  if(storeSuccessiveWs) {
    # Initialize matrix containing w_i(s) values over grid, for each successive estimate
    extraWs = ifelse(useLastPtPProc && (successiveMethod == "pproc"), 0, 1)
    allWs = matrix(nrow=nrow(gridDat), ncol=(nrow(xObs)-minN+1 + extraWs))
  } else {
    allWs = NULL
  }
  
  print(paste0("starting to fit successive model..."))
  for(i in minN:nrow(obsPts)) {
    print(paste0("iteration ", i, "/", nrow(xObs), "..."))
    
    # use only the first i well points for model fitting/estimation
    thisObsCoords = matrix(obsPts[1:i,], nrow=i)
    thisObsValues = obsVals[1:i]
    thisXobs = matrix(xObs[1:i,], nrow=i)
    
    if(successiveMethod == "spde") {
      # fit a full SPDE model each iteration
      
      fit = fitSPDE(obsCoords=as.matrix(thisObsCoords), obsValues=thisObsValues, xObs=thisXobs, 
                    predCoords=as.matrix(predPts), xPred=xPred, 
                    transform=transform, invTransform=invTransform, 
                    mesh=mesh, prior=prior, customFixedI=customFixedI, 
                    significanceCI=significanceCI, int.strategy=int.strategy, strategy=strategy, 
                    nPostSamples=nPostSamples, verbose=verbose, link=link, seed=seed, 
                    family=family, doModAssess=doModAssess, 
                    previousFit=previousFit, improperCovariatePrior=improperCovariatePrior, 
                    fixedParameters=fixedParameters, experimentalMode=experimentalMode, ...)
      
      # update previousFit for faster optimization/model fitting in next iteration
      previousFit = fit$mod
      
      nextWs = transform(fit$customPredEst) # predictions minus effect due to seismic data
    } else {
      # successiveMethod == "pproc". I.e., we fit a density estimator
      
      if(i != nrow(xObs)) {
        # If this isn't the final fit:
        
        # fit the density estimator
        fitPProc = getDensityDueToWells(seismicDat=gridDat, wellDat=thisObsCoords, 
                                        method=pProcMethod, 
                                        centerScale=TRUE, kde.args=kde.args)
        # list(wellEffectGrid=wellEffect, wellEffectPts=wellEffectWells)
        nextWs = fitPProc$wellEffectGrid
      } else {
        # If this is the final fit:
        
        # first fit the density estimator one last time
        fitPProc = getDensityDueToWells(seismicDat=gridDat, wellDat=thisObsCoords, 
                                        method=pProcMethod, 
                                        centerScale=TRUE, kde.args=kde.args)
        nextWs = fitPProc$wellEffectGrid
        
        # final update of w at prediction points before final model
        xPred[,wCol] = fitPProc$wellEffectGrid
        
        if(useLastPtPProc) {
          # In this case, we still need to set w for the last observation
          xObs[i,wCol] = bilinearInterp(matrix(obsPts[i,], nrow=1), cbind(predPts, nextWs))
        }
        
        # get final model fit
        fit = fitSPDE(obsCoords=as.matrix(obsPts), obsValues=obsVals, xObs=xObs, 
                      predCoords=as.matrix(predPts), xPred=xPred, 
                      transform=transform, invTransform=invTransform, 
                      mesh=mesh, prior=prior, customFixedI=customFixedI, 
                      significanceCI=significanceCI, int.strategy=int.strategy, strategy=strategy, 
                      nPostSamples=nPostSamples, verbose=verbose, link=link, seed=seed, 
                      family=family, doModAssess=doModAssess, 
                      previousFit=previousFit, improperCovariatePrior=improperCovariatePrior, 
                      fixedParameters=fixedParameters, experimentalMode=experimentalMode, ...)
      }
    }
    
    # update running definition of w over grid
    xPred[,wCol] = nextWs
    
    # add w_i estimates over prediction grid to allWs
    if(storeSuccessiveWs) {
      allWs[,i-minN+1+extraWs] = xPred[,wCol]
    }
    
    # update xObs w column. First interpolate prediction ws from grid to pts
    if(successiveMethod == "pproc" && useLastPtPProc) {
      if(i == minN) {
        pts = obsPts
      } else {
        pts = matrix(obsPts[i:nrow(obsPts),], ncol=2)
      }
    } else {
      if(i != nrow(obsPts)) {
        pts = matrix(obsPts[(i+1):nrow(obsPts),], ncol=2)
      } else {
        pts = NULL
      }
    }
    
    if(!is.null(pts)) {
      # interpolate grid ws onto the rest of the dataset sampled later
      interpW = bilinearInterp(pts, cbind(predPts, nextWs))
      
      # update xObs
      inds = (1:length(interpW)) + nrow(xObs) - length(interpW)
      xObs[inds, wCol] = interpW
    }
    
  }
  
  c(list(allWs=allWs, obsWs=xObs[,wCol]), fit)
}









