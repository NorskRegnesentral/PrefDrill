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
                                   nPostSamples=1000, verbose=TRUE, link=1, seed=NULL, 
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
    xObs[,2] = NA # add in w (successive spatial estimates that start at 0)
    xObs[1:(minN+1), 2] = priorW
    colnames(xObs) = c("int", "w", gridColNames[-(1:2)])
    customFixedI = c(rep(1, 2), rep(0, ncol(gridDat)-2))
  } else {
    wCol = 1
    seismicCols = (3:ncol(gridDat)) - 1
    xObs[,1] = NA # add in w (successive spatial estimates that start at 0)
    xObs[1:(minN+2), 1] = priorW
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
    
    # set prediction points
    if(i != nrow(obsPts)) {
      if(useLastPtPProc && (successiveMethod == "pproc")) {
        # If we're predicting at a well point location given we know we've 
        # selected that location. I.e., predict at point i based on spatial info 
        # from points 1-i
        
        # predict at the latest observed well point location(s) and at full pred Grid
        if(i == minN) {
          thisPredCoords = matrix(obsPts[1:minN,], nrow=minN)
          colnames(thisPredCoords) = colnames(predPts)
        } else {
          thisPredCoords = matrix(obsPts[i,], nrow=1)
          colnames(thisPredCoords) = colnames(predPts)
        }
      } else {
        # predict at the next (unobserved) well point and at full pred Grid
        
        if(successiveMethod == "spde") {
          thisPredCoords = matrix(obsPts[i+1,], nrow=1)
          thisXpred = matrix(xObs[i+1,], nrow=1)
          if(i != nrow(obsPts)-1) {
            # in addition to predicting at the next point, also need to predict 
            # at the point after to get
            # w_i(s_{i+1}) = eta_{i-1}(s_{i+1}) - beta_{i-1} u(s_{i+1})
            thisPredCoords = rbind(thisPredCoords, obsPts[i+2,])
            thisXpred = rbind(thisXpred, xObs[i+2,])
          } else {
            # second to last iteration
            # in addition to predicting at the last point, also need to predict 
            # at the full prediction grid to get
            # w_n(s_{n+1}) = eta_{n-1}(s_{n+1}) - beta_{n-1} u(s_{n+1})
            if(! storeSuccessiveWs) {
              thisPredCoords = rbind(thisPredCoords, predPts)
              thisXpred = rbind(thisXpred, xPred)
            }
          }
          colnames(thisPredCoords) = colnames(predPts)
          
          if(storeSuccessiveWs && i != nrow(obsPts)-1) {
            # if storing w for each iteration, and it's not second to last iter 
            # (where we've already added grid points to pred pts), add grid pts 
            # to pred pts
            thisPredCoords = rbind(thisPredCoords, predPts)
            thisXpred = rbind(thisXpred, xPred)
          }
        } else {
          # If using a kde
          
          # no need to include grid points with density estimator, since already 
          # being estimated there
          thisPredCoords = matrix(obsPts[i+1,], nrow=1)
          thisXpred = matrix(xObs[i+1,], nrow=1)
          
          if(i != nrow(obsPts)-1) {
            # in addition to predicting at the next point, also need to predict 
            # at the point after to get
            # w_i(s_{i+1}) = eta_{i-1}(s_{i+1}) - beta_{i-1} u(s_{i+1})
            thisPredCoords = rbind(thisPredCoords, obsPts[i+2,])
            thisXpred = rbind(thisXpred, matrix(xObs[i+2,], nrow=1))
          }
        }
        
        colnames(thisPredCoords) = colnames(predPts)
        colnames(thisXpred) = colnames(xPred)
        
        # 
        # thisXpred = rbind(matrix(xObs[i+1,], nrow=1), xPred)
      }
      
      thisCustomFixedI = customFixedI
    } else {
      # This is the last iteration (i == n)
      
      if(useLastPtPProc && (successiveMethod == "pproc")) {
        # for the last fit, predict at the latest observed well point location(s) 
        # and at full pred Grid. Don't need to include the full pred grid here 
        # for the kde fit, since we add them in for the final SPDE model fit
        thisPredCoords = matrix(obsPts[i,], nrow=1)
        colnames(thisPredCoords) = colnames(predPts)
        thisXpred = matrix(xObs[i,], nrow=1)
        colnames(thisXpred) = colnames(xPred)
        thisCustomFixedI = NULL
      } else if(successiveMethod == "pproc") {
        # for the last fit, generate predictions over whole grid
        thisPredCoords = NULL # no need to include grid points explicitly
        thisXpred = xPred
        thisCustomFixedI = NULL
      } else {
        # for the last fit, generate predictions over whole grid
        thisPredCoords = predPts
        thisXpred = xPred
        thisCustomFixedI = NULL
      }
      
    }
    
    if(successiveMethod == "spde") {
      # fit a full SPDE model each iteration
      
      fit = fitSPDE(obsCoords=as.matrix(thisObsCoords), obsValues=thisObsValues, xObs=thisXobs, 
                    predCoords=as.matrix(thisPredCoords), xPred=thisXpred, 
                    transform=transform, invTransform=invTransform, 
                    mesh=mesh, prior=prior, customFixedI=thisCustomFixedI, 
                    significanceCI=significanceCI, int.strategy=int.strategy, strategy=strategy, 
                    nPostSamples=nPostSamples, verbose=verbose, link=link, seed=seed, 
                    family=family, doModAssess=doModAssess, 
                    previousFit=previousFit, improperCovariatePrior=improperCovariatePrior, 
                    fixedParameters=fixedParameters, experimentalMode=experimentalMode, ...)
      
      # update previousFit for faster optimization/model fitting in next iteration
      previousFit = fit$mod
      
      nextWs = transform(fit$customPredEst) # predictions minus effect due to seismic data
      
      if(i != nrow(obsPts)) {
        # If this isn't the final fit:
        
        # update w (covariate representing information from previous point data)
        browser()
        xObs[i+1,wCol] = nextWs[1]
        nextWs = nextWs[-1]
        if(i != (nrow(obsPts)-1)) {
          # in addition to predicting at the next point, also need to predict 
          # at the point after to get
          # w_i(s_{i+1}) = eta_{i-1}(s_{i+1}) - beta_{i-1} u(s_{i+1})
          xObs[i+2,wCol] = nextWs[1]
          nextWs = nextWs[-1]
        } else {
          # in the second to last iteration, must save all predictions over the grid. 
          # no need to set xPred here, since this is done at the end of the for loop. 
          # nextWs already consists of predictions over the grid.
          # DO NOTHING
        }
        
      } else {
        # If this is the final fit:
        
        # nextWs already consists of only w at prediction coords.
        # DO NOTHING
      }
    } else {
      # successiveMethod == "pproc". I.e., we fit a density estimator
      
      if(i != nrow(xObs)) {
        # If this isn't the final fit:
        
        # fit the density estimator
        fitPProc = getDensityDueToWells(seismicDat=gridDat, wellDat=thisObsCoords, 
                                        predPts=thisPredCoords, method=pProcMethod, 
                                        centerScale=TRUE, kde.args=kde.args)
        # list(wellEffectGrid=wellEffect, wellEffectPts=wellEffectWells)
        nextWs = fitPProc$wellEffectGrid
        
        # update w (covariate representing information from previous point data)
        # browser()
        if(!useLastPtPProc) {
          xObs[i+1,wCol] = fitPProc$wellEffectPredPts
        } else if(i == minN) {
          xObs[1:minN,wCol] = fitPProc$wellEffectPredPts
        } else {
          xObs[i,wCol] = fitPProc$wellEffectPredPts[1]
        }
      } else {
        # If this is the final fit:
        
        # first fit the density estimator one last time
        fitPProc = getDensityDueToWells(seismicDat=gridDat, wellDat=thisObsCoords, 
                                        predPts=thisPredCoords, method=pProcMethod, 
                                        centerScale=TRUE, kde.args=kde.args)
        nextWs = fitPProc$wellEffectGrid
        
        # final update of w at prediction points before final model
        xPred[,wCol] = fitPProc$wellEffectGrid
        
        if(useLastPtPProc) {
          # In this case, we still need to set w for the last observation
          xObs[i,wCol] = fitPProc$wellEffectPredPts
        }
        
        # get final model fit
        fit = fitSPDE(obsCoords=as.matrix(obsPts), obsValues=obsVals, xObs=xObs, 
                      predCoords=as.matrix(predPts), xPred=xPred, 
                      transform=transform, invTransform=invTransform, 
                      mesh=mesh, prior=prior, customFixedI=thisCustomFixedI, 
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
  }
  
  c(list(allWs=allWs, obsWs=xObs[,wCol]), fit)
}









