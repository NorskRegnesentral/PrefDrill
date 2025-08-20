

# # for sampling wells successively, iterating between sampling and predicting, 
# # and where samples are taken preferentially wrt the prediction mean.
# wellSampler = function(truthDat, seismicDat, modelFitter, nWells=20, minN=4, 
#                        predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
#                        transform=logit, invTransform=expit, prefPar=0, 
#                        samplingModel=c("ipp"), sigmaSqErr=.1^2, 
#                        repelType=c("none", "rbf"), 
#                        bwRepel=NULL, rbf=c("uniform", "gaussian", "exp"), 
#                        repelAmount=NULL, seed=NULL, ...) {
#   
#   if(!is.null(seed)) {
#     set.seed(seed)
#   }
#   if(is.null(bwRepel)) {
#     bwRepel = max(diff(range(seismicDat$east)), diff(range(seismicDat$north)))/20
#   }
#   
#   samplingModel = match.arg(samplingModel)
#   repelType = match.arg(repelType)
#   rbf = match.arg(rbf)
#   
#   # get first datapoints based on only seismic data
#   print("Initial well sampling...")
#   initWellDat = basicWellSampler(nWells=minN, wellDat=NULL, seismicDat=seismicDat, 
#                                  truthDat=truthDat, modelFitter=modelFitter, 
#                                  predGrid=predGrid, prefPar=prefPar, 
#                                  transform=transform, invTransform=invTransform, 
#                                  samplingModel=samplingModel, sigmaSqErr=sigmaSqErr, 
#                                  repelType=repelType, bwRepel=bwRepel, rbf=rbf, 
#                                  repelAmount=repelAmount, seed=NULL, ...)$wellDat
#   wellDat = data.frame(east=initWellDat[,1], north=initWellDat[,2], volFrac=initWellDat[,3])
#   
#   # # now loop through, sampling and fitting models. Capture previous model fits 
#   # # to hopefully make things go faster
#   # print("Successive sampling...")
#   # startT = proc.time()[3]
#   # prevFit = NULL
#   # for(i in (minN+1):nWells) {
#   #   print(paste0("sampling well ", i, "/", nWells))
#   #   sampT = system.time(out <- basicWellSampler(nWells=1, wellDat=wellDat, seismicDat=seismicDat, 
#   #                                 truthDat=truthDat, modelFitter=modelFitter, 
#   #                                 predGrid=predGrid, prefPar=prefPar, 
#   #                                 transform=transform, invTransform=invTransform, 
#   #                                 samplingModel=samplingModel, sigmaSqErr=sigmaSqErr, 
#   #                                 repelType=repelType, bwRepel=bwRepel, rbf=rbf, 
#   #                                 repelAmount=repelAmount, seed=NULL, previousFit=prevFit))[3]
#   #   newWellDat = out$wellDat
#   #   prevFit = out$fit$mod
#   #   print(paste0("took ", sampT, " seconds"))
#   #   
#   #   wellDat = rbind(wellDat, newWellDat)
#   # }
#   # endT = proc.time()[3]
#   # print(paste0("total time: ", (endT - startT)/60, " minutes"))
#   
#   sampT = system.time(out <- sequentialWellSampler(nWells=nWells-minN, wellDat=wellDat, seismicDat=seismicDat, 
#                                 truthDat=truthDat, modelFitter=modelFitter, 
#                                 predGrid=predGrid, prefPar=prefPar, 
#                                 transform=transform, invTransform=invTransform, 
#                                 samplingModel=samplingModel, sigmaSqErr=sigmaSqErr, 
#                                 repelType=repelType, bwRepel=bwRepel, rbf=rbf, 
#                                 repelAmount=repelAmount, seed=NULL, previousFit=prevFit))[3]
#   
#   wellDat
# }

# for sampling wells from a fixed inhomogeneous Poisson process
# 
# NOTE: must sample points one at a time if there is any repulsion
# 
# Inputs:
# modelFitter: If not null, uses predictions based on wellDat to get sampling 
#              intensity. Else use seismicDat
basicWellSampler = function(nWells=1, wellDat=NULL, seismicDat, truthDat=NULL, modelFitter=NULL, 
                            predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
                            preds=NULL, transform=logit, invTransform=expit, prefPar=0, 
                            samplingModel=c("ipp"), sigmaSqErr=.1^2, 
                            repelType=c("none", "rbf"), bwRepel=NULL, 
                            rbf=c("uniform", "gaussian", "exp"), repelAmount=NULL, 
                            seed=NULL, isWatson=FALSE, ...) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  if(is.null(bwRepel)) {
    bwRepel = max(diff(range(seismicDat$east)), diff(range(seismicDat$north)))/20
  }
  
  samplingModel = match.arg(samplingModel)
  repelType = match.arg(repelType)
  rbf = match.arg(rbf)
  
  # in the case where we have repulsion and we need to sample multiple wells, 
  # sample them one at a time
  if((nWells > 1) && (repelType != "none")) {
    tmpWellDat = wellDat
    for(i in 1:nWells) {
      out = basicWellSampler(nWells=1, wellDat=wellDat, seismicDat=seismicDat, 
                             truthDat=truthDat, modelFitter=modelFitter, 
                             predGrid=predGrid, preds=preds, 
                             transform=transform, invTransform=invTransform, prefPar=prefPar, 
                             samplingModel=samplingModel, sigmaSqErr=sigmaSqErr, 
                             repelType=repelType, bwRepel=bwRepel, 
                             rbf=rbf, repelAmount=repelAmount, 
                             seed=seed, isWatson=isWatson, ...)
      # list(wellDat = as.data.frame(res), preds=preds, predAggMat=predAggMat, fit=fit)
      
      thisWell = out$wellDat
      
      if(i == 1) {
        allWells = thisWell
      } else {
        allWells = rbind(allWells, thisWell)
      }
      wellDat = rbind(wellDat, thisWell)
    }
    
    out = list(wellDat = wellDat, logitProbsNoRep = out$logitProbsNoRep)
    
    return(out)
  }
  
  if(is.null(repelAmount)) {
    repelAmount = ifelse(rbf == "uniform", 10, 1)
  }
  
  repelKern = repelKernel(repelType=repelType, bw=bwRepel, rbf=rbf)
  
  # first fit the model and get predictions over the grid (or just take seismic estimates)
  if(is.null(preds)) {
    if(!is.null(wellDat) && !is.null(modelFitter)) {
      if(!isWatson) {
        fit = modelFitter(wellDat=wellDat, seismicDat=seismicDat, predGrid=predGrid, 
                          transform=transform, invTransform=invTransform, ...)
      } else {
        fit = modelFitter(wellDat=wellDat, seismicDat=seismicDat, predGrid=predGrid, 
                          transform=transform, invTransform=invTransform, repelDist=bwRepel, ...)
      }
      
      preds = fit$predEst
      predAggMat = fit$predAggMat
    } else if(!is.null(seismicDat)) {
      fit = NULL
      preds = seismicDat[,3]
      predAggMat = mean(preds)
    } else if(prefPar == 0) {
      fit = NULL
      preds = rep(0.5, nrow(predGrid))
      predAggMat = mean(preds)
    } else {
      stop("preds, seismicDat, and modelFitter are all NULL, and prefPar != 0")
    }
  } else {
    fit = NULL
    predAggMat = mean(preds)
  }
  
  if(any(preds <= 0)) {
    stop("some sampling probabilities are negative. Check data and predictions on correct scale.")
  }
  
  # calculate selection probabilities over the grid based on predictions
  if(samplingModel == "ipp") {
    # in the inhomogeneous Poisson process, sampling probability is proportional 
    # to the predictions
    logitProbsNoRep = logit(preds)*prefPar
  } else {
    stop(paste0("sampling model ", samplingModel, " not yet supported"))
  }
  
  # update selection probabilities to account for repulsion effects
  if(!is.null(wellDat) && repelType != "none") {
    # calculate cross-distance matrix from grid to well points
    distMat = rdist(predGrid[,1:2], as.matrix(wellDat[,1:2]))
    
    # calculate repulsions from each well point to the grid
    repelEffectsMat = matrix(repelKern(distMat), ncol=nrow(wellDat))
    
    # for each grid point, take sum of all repulsion effects together to get 
    # total repulsion (define 0 * Inf = 0)
    tempSum = apply(repelEffectsMat, 1, sum)
    repelEffectsGrid = tempSum * repelAmount
    forceZero = (tempSum == 0) & (is.infinite(repelAmount))
    repelEffectsGrid[forceZero] = 0
    
    # update sampling probabilities
    logitProbs = logitProbsNoRep - repelEffectsGrid
  }
  
  # transform probabilities back to probability scale, ensure sum to 1 over grid
  probs = expit(logitProbs)
  probs = probs * (1 / sum(probs))
  
  # sample the grid cell
  nextIs = sample(1:length(probs), nWells, prob = probs, replace=TRUE)
  
  # sample the point uniformly within the grid cell
  # NOTE: with strict repulsion may need to modify this step
  deltaX = diff(sort(unique(predGrid[,1])))[1]
  deltaY = diff(sort(unique(predGrid[,2])))[1]
  delta = c(deltaX, deltaY)
  centerPts = matrix(predGrid[nextIs,], nrow=nWells)
  ptErrs = matrix((runif(2*nWells) - 0.5) * rep(delta, nWells), byrow=TRUE, ncol=2)
  nextPts = centerPts + ptErrs
  
  # get the sampling probabilities without repulsion for the well points
  nextLogitProbsNoRep = logitProbsNoRep[nextIs]
  
  if(!is.null(truthDat)) {
    # get seismic estimate at that point
    nextTruth = bilinearInterp(nextPts, truthDat, 
                               transform=transform, invTransform=invTransform)
    nextVolFrac = invTransform(transform(nextTruth) + rnorm(nWells, sd=sqrt(sigmaSqErr)))
    res = cbind(nextPts, nextVolFrac, nextIs, nextLogitProbsNoRep)
    colnames(res) = c("east", "north", "volFrac", "gridI", "logitProbsNoRep")
  } else {
    res = cbind(nextPts, nextIs, nextLogitProbsNoRep)
    colnames(res) = c("east", "north", "gridI", "logitProbsNoRep")
  }
  
  list(wellDat = as.data.frame(res), preds=preds, predAggMat=predAggMat, 
       logitProbsNoRep=logitProbsNoRep, fit=fit)
}


# for either batch or sequential well sampling depending on setting
# 
# NOTE: must sample points one at a time if there is any repulsion
# 
# Inputs:
# modelFitter: If not null, uses predictions based on wellDat to get sampling 
#              intensity. Else use seismicDat
wellSampler = function(nWells=1, wellDat=NULL, seismicDat, truthDat=NULL, modelFitter=NULL, 
                       predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
                       preds=NULL, transform=logit, invTransform=expit, prefPar=0, 
                       samplingModel=c("ipp"), sigmaSqErr=.1^2, 
                       repelType=c("none", "rbf"), bwRepel=NULL, 
                       rbf=c("uniform", "gaussian", "exp"), repelAmount=NULL, 
                       seed=NULL, batchSize=1, minN=4, isWatson=FALSE, verbose=FALSE, 
                       saveAllPredsProbs=TRUE, getProbsNoRepOnly=TRUE, ...) {
  
  # set up defaults
  if(!is.null(seed)) {
    set.seed(seed)
  }
  if(is.null(bwRepel)) {
    bwRepel = max(diff(range(seismicDat$east)), diff(range(seismicDat$north)))/20
  }
  
  samplingModel = match.arg(samplingModel)
  repelType = match.arg(repelType)
  rbf = match.arg(rbf)
  
  if(is.null(repelAmount)) {
    repelAmount = ifelse(rbf == "uniform", 10, 1)
  }
  
  repelKern = repelKernel(repelType=repelType, bw=bwRepel, rbf=rbf)
  
  # if we don't need to do one at a time sampling, don't
  if((repelType == "none" || bwRepel == 0) && is.null(modelFitter)) {
    return(basicWellSampler(nWells=nWells, wellDat=wellDat, seismicDat=seismicDat, 
                            truthDat=truthDat, modelFitter=modelFitter, 
                            predGrid=predGrid, preds=preds, 
                            transform=transform, invTransform=invTransform, prefPar=prefPar, 
                            samplingModel=samplingModel, sigmaSqErr=sigmaSqErr, 
                            repelType=repelType, bwRepel=bwRepel, 
                            rbf=rbf, repelAmount=repelAmount, isWatson=isWatson, 
                            seed=seed, ...))
  }
  
  # calculate grid resolution
  deltaX = diff(sort(unique(predGrid[,1])))[1]
  deltaY = diff(sort(unique(predGrid[,2])))[1]
  delta = c(deltaX, deltaY)
  
  # precompute sampling probabilities based on preds if given as input
  if(!is.null(preds)) {
    
    # calculate selection probabilities over the grid based on predictions
    if(samplingModel == "ipp") {
      
      # in the inhomogeneous Poisson process, sampling probability is proportional 
      # to the predictions
      logitProbs = logit(preds) * prefPar
      
    } else {
      stop(paste0("sampling model ", samplingModel, " not yet supported"))
    }
    
    # base sampling intensity on input preds field
    fit = NULL
    predAggMat = mean(preds)
    
  } else {
    logitProbs = NULL
  }
  
  print("Successive sampling...")
  startT = proc.time()[3]
  
  # in the case where we have repulsion and we need to sample multiple wells, 
  # sample them one at a time, saving computations along the way
  tmpWellDat = wellDat
  distMat = NULL
  prevFit = NULL
  actualMinN = min(max(c(minN, batchSize)), nWells)
  nbatch = ceiling((nWells-actualMinN)/batchSize) + 1
  if(saveAllPredsProbs) {
    if(getProbsNoRepOnly) {
      allPreds = NULL
      allProbs = NULL
    } else {
      allPreds = matrix(nrow=nrow(predGrid), ncol=nbatch)
      allProbs = matrix(nrow=nrow(predGrid), ncol=nbatch)
    }
    
    allLogitProbsNoRep = matrix(nrow=nrow(predGrid), ncol=nbatch)
  }
  for(i in 1:nbatch) {
    if(verbose) {
      startBatchT = proc.time()[3]
      print(paste0("sampling batch ", i, "/", nbatch))
    }
    
    if(i == 1) {
      # the first batch is of size actualMinN
      thisBatchSize = min(c(actualMinN, nWells))
      
    } else if(i == nbatch) {
      # the last batch size is leftover wells. Must account for actualMinN batch and nbatch-2 regular batches
      thisBatchSize = nWells - actualMinN - (nbatch-2)*batchSize
      
    } else {
      
      thisBatchSize = batchSize
    }
    
    if(i == 2) {
      lastBatchSize = actualMinN
    } else {
      lastBatchSize = batchSize # this variable not used if i==1, so this is always true
    }
    
    batchInds = 1:thisBatchSize
    if(!is.null(wellDat)) {
      lastBatchInds = (nrow(wellDat) - lastBatchSize + 1):nrow(wellDat)
    } else  {
      lastBatchInds = NULL
    }
    
    # first get sampling intensity for this batch without repulsion
    
    if(!is.null(wellDat) && !is.null(modelFitter)) {
      wellDatDF = as.data.frame(wellDat)
      
      # base sampling intensity on model predictions
      if(!isWatson) {
        fit = modelFitter(wellDat=wellDatDF, seismicDat=seismicDat, predGrid=predGrid, 
                          transform=transform, invTransform=invTransform, previousFit=prevFit, 
                          ...)
        
        preds = fit$predEst
        predAggMat = fit$predAggMat
      } else {
        fit = modelFitter(wellDat=wellDatDF, seismicDat=seismicDat, predGrid=predGrid, repelDist=bwRepel, 
                          transform=transform, invTransform=invTransform, previousFit=prevFit, 
                          ...)
        
        preds = fit$pred.yEst
        predAggMat = fit$pred.yAggMat
      }
      
      
      
      # save this fit to use next time
      prevFit = fit$mod
      
    } else if(!is.null(seismicDat)) {
      
      # base sampling intensity on seismic data
      fit = NULL
      preds = seismicDat[,3]
      predAggMat = mean(preds)
      
    } else if(prefPar == 0) {
      
      # no preferentiality, so start sampling intensity as uniform
      fit = NULL
      preds = rep(0.5, nrow(predGrid))
      predAggMat = mean(preds)
      
    } else {
      stop("preds, seismicDat, and modelFitter are all NULL, and prefPar != 0")
    }
    
    # calculate selection probabilities over the grid based on predictions
    if(samplingModel == "ipp") {
      
      # in the inhomogeneous Poisson process, sampling probability is proportional 
      # to the predictions
      logitProbsNoRep = logit(preds) * prefPar
      
    } else {
      stop(paste0("sampling model ", samplingModel, " not yet supported"))
    }
    
    if(any(preds <= 0)) {
      stop("some sampling probabilities are negative. Check data and predictions on correct scale.")
    }
    
    if(saveAllPredsProbs && !getProbsNoRepOnly) {
      allPreds[,i] = preds
    }
    
    # while selection probabilities prior to repulsion are calculated in batches, 
    # repulsion must be accounted for one at a time
    for(j in 1:thisBatchSize) {
      
      if((repelType == "none" || bwRepel == 0)) {
        
        # no repel effects
        repelEffectsGrid = 0
        
      } else {
        
        # calculate repulsion effects
        if((i == 1) && (j == 1)) {
          
          if(!is.null(wellDat)) {
            
            # calculate cross-distance matrix from grid to well points
            distMat = rdist(predGrid[,1:2], as.matrix(wellDat[,1:2]))
            
            # calculate repulsions from each well point to the grid
            repelEffectsMat = matrix(repelKern(distMat), ncol=nrow(wellDat))
            
            # for each grid point, take sum of all repulsion effects together to get 
            # total repulsion (define 0 * Inf = 0)
            tempSum = apply(repelEffectsMat, 1, sum)
            repelEffectsGrid = tempSum * repelAmount
            forceZero = (tempSum == 0) & (is.infinite(repelAmount))
            repelEffectsGrid[forceZero] = 0
          } else {
            
            # no wells to have any repulsion effects from yet
            repelEffectsGrid = rep(0, nrow(predGrid))
          }
          
          
        } else {
          
          # only need to calculate distances from the last point to the grid
          newDistMat = rdist(predGrid[,1:2], matrix(as.matrix(wellDat[nrow(wellDat), 1:2]), nrow=1))
          
          # calculate repulsions from each well point to the grid
          newRepelEffectsMat = matrix(repelKern(newDistMat), ncol=1)
          
          # for each grid point, take sum of all repulsion effects together to get 
          # total repulsion (define 0 * Inf = 0)
          tempSum = newRepelEffectsMat
          newRepelEffectsGrid = tempSum * repelAmount
          forceZero = (tempSum == 0) & (is.infinite(repelAmount))
          newRepelEffectsGrid[forceZero] = 0
          
          # add on the new repulsion effects to the running total
          repelEffectsGrid = repelEffectsGrid + newRepelEffectsGrid
        }
      }
      
      
      # update selection probabilities to account for repulsion effects
      thisLogitProbs = logitProbsNoRep - repelEffectsGrid
      
      # transform probabilities back to probability scale, ensure sum to 1 over grid
      probs = expit(thisLogitProbs)
      probs = probs * (1 / sum(probs))
      
      if((j == 1) && saveAllPredsProbs) {
        if(!getProbsNoRepOnly) {
          allProbs[,i] = probs
        }
        allLogitProbsNoRep[,i] = logitProbsNoRep
      }
      
      # sample the grid cell
      nextIs = sample(1:length(probs), 1, prob = probs, replace=TRUE)
      
      # sample the point uniformly within the grid cell
      # NOTE: with strict repulsion may need to modify this step
      centerPts = matrix(predGrid[nextIs,], nrow=1)
      ptErrs = matrix((runif(2*1) - 0.5) * rep(delta, 1), byrow=TRUE, ncol=2)
      nextPts = centerPts + ptErrs
      
      # get the sampling probabilities without repulsion for the well points
      nextLogitProbsNoRep = logitProbsNoRep[nextIs]
      
      if(!is.null(truthDat)) {
        # get seismic estimate at that point
        nextTruth = bilinearInterp(nextPts, truthDat, 
                                   transform=transform, invTransform=invTransform)
        
        nextVolFrac = invTransform(transform(nextTruth) + rnorm(1, sd=sqrt(sigmaSqErr)))
        thisWell = cbind(nextPts, nextVolFrac, nextIs, nextLogitProbsNoRep)
        colnames(thisWell) = c("east", "north", "volFrac", "gridI", "logitProbsNoRep")
      } else {
        thisWell = cbind(nextPts, nextIs, nextLogitProbsNoRep)
        colnames(thisWell) = c("east", "north", "gridI", "logitProbsNoRep")
      }
      
      # update well data with this batch of samples
      wellDat = rbind(wellDat, thisWell)
    }
    
    if(verbose) {
      endBatchT = proc.time()[3]
      print(paste0("took ", endBatchT - startBatchT, " seconds"))
    }
  }
  
  endT = proc.time()[3]
  print(paste0("total time: ", (endT - startT)/60, " minutes"))
  
  names(preds) = NULL # remove names to save space
  if(saveAllPredsProbs) {
    list(wellDat = as.data.frame(wellDat), preds=preds, predAggMat=predAggMat, 
         allPreds=allPreds, allProbs=allProbs, allLogitProbsNoRep=allLogitProbsNoRep)
  } else {
    list(wellDat = as.data.frame(wellDat), preds=preds, predAggMat=predAggMat)
  }
}

# 0 if no repulsion, 1 if max repulsion
repelKernel = function(repelType=c("none", "rbf"), 
                       bw=1, rbf=c("uniform", "gaussian", "exp")) {
  
  repelType = match.arg(repelType)
  rbf = match.arg(rbf)
  
  if(repelType == "none") {
    kern = function(d) {
      rep(0, length(d))
    }
  } else if(repelType == "rbf") {
    if(rbf == "gaussian") {
      subKern = function(d) { exp(-(d/bw)^2) }
    } else if(rbf == "uniform") {
      kern = function(d) {
        out = rep(0, length(d))
        out[d < bw] = 1
        out
      }
      subKern = NULL
    } else if(rbf == "exponential") {
      subKern = function(d) { exp(-d/bw) }
    }
    if(!is.null(subKern)) {
      kern = function(d) {
        out = subKern(d)
      }
    }
  }
  
  kern
}








