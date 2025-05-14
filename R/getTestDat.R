# make a single simple set of datasets for testing purposes


# NOTE: this is not supposed to look like the simulation study data...
getTestDat = function(n=20, xlim=simStudyXlims, ylim=simStudyYlims, 
                      muSpat=0, sigmaSqSpat=1, corrSeismic=.6, sigmaSqErr=.01, 
                      aRange=diff(xlim)/4, prefPar=2.5, gridWidth=diff(xlim)/50, 
                      grid=make.surface.grid(list(east=seq(xlim[1]+gridWidth/2, xlim[2]-gridWidth/2, l=50), 
                                                  north=seq(ylim[1]+gridWidth/2, ylim[2]-gridWidth/2, l=50))), 
                      prefType=c("seismic", "truth"), 
                      repelType=c("none", "hardcore", "rbf"), 
                      bwRepel=max(diff(range(seismicDat$east)), diff(range(seismicDat$north)))/100, 
                      rbf=c("gaussian", "exp"), seed=1, saveDat=TRUE) {
  require(fields)
  set.seed(seed)
  prefType = match.arg(prefType)
  repelType = match.arg(repelType)
  rbf = match.arg(rbf)
  
  # simulate spatial fields for calculating truth and seismic data
  Sigma = stationary.cov(grid, aRange=aRange)*sigmaSqSpat
  L = t(chol(Sigma))
  truth = L %*% rnorm(nrow(L))
  indepPt = L %*% rnorm(nrow(L))
  seismicEst = corrSeismic * truth + sqrt(1-corrSeismic) * indepPt
  
  # add in mean
  truth = truth + muSpat
  seismicEst = seismicEst + muSpat
  
  # sample well point locations at grid centroids based on seismic estimates
  if(prefType == "seismic") {
    lambdas = exp(prefPar * seismicEst)
  } else if (prefType == "truth") {
    lambdas = exp(prefPar * truth)
  } else {
    stop("prefType must be 'seismic' or 'truth'")
  }
  sampleProbs = lambdas * (1/sum(lambdas))
  sampleI = sample(1:length(lambdas), size = n, replace=T, prob=sampleProbs)
  
  # now sample well point locations uniformly within the cells
  gridWidth = diff(sort(unique(grid[,1]))[1:2])
  pts = grid[sampleI,] + matrix(runif(2*n)*gridWidth - gridWidth/2, ncol=2)
  
  # bilinearly interpolate the truth at the well points and add nugget effect
  wellVolFrac = bilinearInterp(pts, cbind(grid, truth)) + rnorm(n, sd=sqrt(sigmaSqErr))
  
  # transform the truth and seismic data to correspond to fractions
  truthTrans = expit(truth)
  seismicEstTrans = expit(seismicEst)
  wellVolFracTrans = expit(wellVolFrac)
  
  # construct final datasets and truth
  truthTestDat = data.frame(east=grid[,1], north=grid[,2], truth=truthTrans)
  seismicTestDat = data.frame(east=grid[,1], north=grid[,2], seismicEst=seismicEstTrans)
  wellTestDat = data.frame(east=pts[,1], north=pts[,2], volFrac=wellVolFracTrans)
  
  if(saveDat) {
    if(prefType == "seismic") {
      save(truthTestDat, seismicTestDat, wellTestDat, file=paste0(globalDir, "/testDat.RData"))
    } else if(prefType == "truth") {
      save(truthTestDat, seismicTestDat, wellTestDat, file=paste0(globalDir, "/testDat_truthPref.RData"))
    } else {
      stop("chosen prefType not supported")
    }
  }
  
  invisible(list(wellDat=wellTestDat, seismicDat=seismicTestDat, truthDat=truthTestDat))
}

# NOTE: this is not supposed to look like the simulation study data...
getSuccessiveTestDat = function(n=20, truthDat=NULL, seismicDat=NULL, 
                                xlim=simStudyXlims, ylim=simStudyYlims, 
                                muSpat=0, sigmaSqSpat=1, corrSeismic=.6, sigmaSqErr=.01, 
                                modelFitter=fitSPDEsimDat, aRange=diff(xlim)/4, prefPar=4, 
                                gridWidth=diff(xlim)/50, transform=logit, invTransform=expit, 
                                grid=make.surface.grid(list(east=seq(xlim[1]+gridWidth/2, xlim[2]-gridWidth/2, l=50), 
                                                            north=seq(ylim[1]+gridWidth/2, ylim[2]-gridWidth/2, l=50))), 
                                samplingModel=c("ipp"), 
                                repelType=c("rbf", "none"), bwRepel=NULL, 
                                rbf=c("uniform", "gaussian", "exp"), repelAmount=NULL, 
                                seed=1, saveDat=TRUE, minN=4, ...) {
  
  require(fields)
  set.seed(seed)
  samplingModel = match.arg(samplingModel)
  repelType = match.arg(repelType)
  rbf = match.arg(rbf)
  
  if(is.null(bwRepel)) {
    # one twentieth of the domain diameter by default
    if(!is.null(seismicDat)) {
      bwRepel = max(diff(range(seismicDat[,1])), diff(range(seismicDat[,2])))/20
    } else {
      bwRepel = max(diff(xlim), diff(ylim))/20
    }
  }
  
  # simulate spatial fields for calculating truth and seismic data
  if(is.null(truthDat)) {
    Sigma = stationary.cov(grid, aRange=aRange)*sigmaSqSpat
    L = t(chol(Sigma))
    truth = L %*% rnorm(nrow(L))
    indepPt = L %*% rnorm(nrow(L))
    seismicEst = corrSeismic * truth + sqrt(1-corrSeismic) * indepPt
    
    # add in mean
    truth = truth + muSpat
    seismicEst = seismicEst + muSpat
    
    # transform the truth and seismic data to correspond to fractions
    truthTrans = invTransform(truth)
    seismicEstTrans = invTransform(seismicEst)
    
    # construct seismic and truth data frames
    truthTestDat = data.frame(east=grid[,1], north=grid[,2], truth=truthTrans)
    seismicTestDat = data.frame(east=grid[,1], north=grid[,2], seismicEst=seismicEstTrans)
  } else {
    # set seismic and truth data frames
    truthTestDat = truthDat
    seismicTestDat = seismicDat
  }
  
  # construct well data via successive sampling process
  wellTestDat = wellSampler(truthDat=truthTestDat, seismicDat=seismicTestDat, 
                            modelFitter=modelFitter, nWells=n, minN=minN, 
                            predGrid=grid, prefPar=prefPar, 
                            transform=transform, invTransform=invTransform, 
                            samplingModel=samplingModel, sigmaSqErr=sigmaSqErr, 
                            repelType=repelType, bwRepel=bwRepel, rbf=rbf, 
                            repelAmount=repelAmount, seed=NULL, ...)
  
  # save datasets if need be
  if(saveDat) {
    if(repelType == "none") {
      repelStr = ""
    } else {
      repelStr = paste0("_", rbf)
    }
    save(truthTestDat, seismicTestDat, wellTestDat, file=paste0(globalDir, "/testSuccessiveDat_", samplingModel, repelStr, ".RData"))
  }
  
  invisible(list(wellDat=wellTestDat, seismicDat=seismicTestDat, truthDat=truthTestDat))
}

wellSampler = function(truthDat, seismicDat, modelFitter, nWells=20, minN=4, 
                       predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
                       transform=logit, invTransform=expit, prefPar=0, 
                       samplingModel=c("ipp"), sigmaSqErr=.1^2, 
                       repelType=c("none", "rbf"), 
                       bwRepel=NULL, rbf=c("uniform", "gaussian", "exp"), 
                       repelAmount=NULL, seed=NULL, ...) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  if(is.null(bwRepel)) {
    bwRepel = max(diff(range(seismicDat$east)), diff(range(seismicDat$north)))/20
  }
  
  samplingModel = match.arg(samplingModel)
  repelType = match.arg(repelType)
  rbf = match.arg(rbf)
  
  # get first datapoints based on only seismic data
  print("Initial well sampling...")
  initWellDat = basicWellSampler(nWells=minN, wellDat=NULL, seismicDat=seismicDat, 
                                 truthDat=truthDat, modelFitter=modelFitter, 
                                 predGrid=predGrid, prefPar=prefPar, 
                                 transform=transform, invTransform=invTransform, 
                                 samplingModel=samplingModel, sigmaSqErr=sigmaSqErr, 
                                 repelType=repelType, bwRepel=bwRepel, rbf=rbf, 
                                 repelAmount=repelAmount, seed=NULL, ...)
  wellDat = data.frame(east=initWellDat[,1], north=initWellDat[,2], volFrac=initWellDat[,3])
  
  # now loop through, sampling and fitting models
  print("Successive sampling...")
  for(i in (minN+1):nWells) {
    print(paste0("sampling well ", i, "/", nWells))
    newWellDat = basicWellSampler(nWells=1, wellDat=wellDat, seismicDat=seismicDat, 
                                  truthDat=truthDat, modelFitter=modelFitter, 
                                  predGrid=predGrid, prefPar=prefPar, 
                                  transform=transform, invTransform=invTransform, 
                                  samplingModel=samplingModel, sigmaSqErr=sigmaSqErr, 
                                  repelType=repelType, bwRepel=bwRepel, rbf=rbf, 
                                  repelAmount=repelAmount, seed=NULL)
    
    wellDat = rbind(wellDat, newWellDat)
  }
  
  wellDat
}

basicWellSampler = function(nWells=1, wellDat=NULL, seismicDat, truthDat=NULL, modelFitter, 
                            predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
                            transform=logit, invTransform=expit, prefPar=0, 
                            samplingModel=c("ipp"), sigmaSqErr=.1^2, 
                            repelType=c("none", "rbf"), bwRepel=NULL, 
                            rbf=c("uniform", "gaussian", "exp"), repelAmount=NULL, 
                            seed=NULL, ...) {
  
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
      thisWell = basicWellSampler(nWells=1, wellDat=wellDat, seismicDat=seismicDat, 
                                  truthDat=truthDat, modelFitter=modelFitter, 
                                  predGrid=predGrid, 
                                  transform=transform, invTransform=invTransform, prefPar=prefPar, 
                                  samplingModel=samplingModel, sigmaSqErr=sigmaSqErr, 
                                  repelType=repelType, bwRepel=bwRepel, 
                                  rbf=rbf, repelAmount=repelAmount, 
                                  seed=seed, ...)
      
      if(i == 1) {
        allWells = thisWell
      } else {
        allWells = rbind(allWells, thisWell)
      }
      wellDat = rbind(wellDat, thisWell)
    }
    
    return(allWells)
  }
  
  if(is.null(repelAmount)) {
    repelAmount = ifelse(rbf == "uniform", 10, 1)
  }
  
  repelKern = repelKernel(repelType=repelType, bw=bwRepel, rbf=rbf)
  
  # first fit the model and get predictions over the grid (or just take seismic estimates)
  if(!is.null(wellDat)) {
    mod = modelFitter(wellDat=wellDat, seismicDat=seismicDat, predGrid=predGrid, 
                      transform=transform, invTransform=invTransform, ...)
    preds = mod$predEst
    predAggMat = mod$predAggMat
  } else {
    preds = seismicDat[,3]
    predAggMat = mean(preds)
  }
  
  if(any(preds <= 0)) {
    stop("some sampling probabilities are negative. Check data and predictions on correct scale.")
  }
  
  # calculate selection probabilities over the grid based on predictions
  if(samplingModel == "ipp") {
    # in the inhomogeneous Poisson process, sampling probability is proportional 
    # to the predictions
    logitProbs = logit(preds)*prefPar
  } else {
    stop(paste0("sampling model ", samplingModel, " not yet supported"))
  }
  
  # update selection probabilities to account for repulsion effects
  if(!is.null(wellDat)) {
    # calculate cross-distance matrix from grid to well points
    distMat = rdist(predGrid[,1:2], as.matrix(wellDat[,1:2]))
    
    # calculate repulsions from each well point to the grid
    repelEffectsMat = matrix(repelKern(distMat), ncol=nrow(wellDat))
    
    # for each grid point, take sum of all repulsion effects together to get 
    # total repulsion
    repelEffectsGrid = apply(repelEffectsMat, 1, sum) * repelAmount
    
    # update sampling probabilities
    logitProbs = logitProbs - repelEffectsGrid
  }
  
  # transform probabilities back to probability scale, ensure sum to 1 over grid
  probs = expit(logitProbs)
  probs = probs * (1 / sum(probs))
  
  # sample the grid cell
  nextIs = sample(1:length(probs), nWells, prob = probs, replace=TRUE)
  
  # sample the point uniformly within the grid cell
  # NOTE: with strict repulsion may need to modify this step
  delta = diff(sort(unique(predGrid[,1]))[1:2])
  centerPts = matrix(predGrid[nextIs,], nrow=nWells)
  nextPts = matrix(centerPts + (runif(2*nWells) - 0.5) * delta, nrow=nWells)
  
  if(!is.null(truthDat)) {
    # get seismic estimate at that point
    nextTruth = bilinearInterp(nextPts, truthDat, 
                               transform=transform, invTransform=invTransform)
    nextVolFrac = invTransform(transform(nextTruth) + rnorm(nWells, sd=sqrt(sigmaSqErr)))
    res = cbind(nextPts, nextVolFrac)
    colnames(res) = c("east", "north", "volFrac")
  } else {
    res = nextPts
    colnames(res) = c("east", "north")
  }
  
  list(wellDat = as.data.frame(res), preds=preds, predAggMat=predAggMat)
}

repelKernel = function(repelType=c("none", "rbf"), 
                       bw=1, rbf=c("uniform", "gaussian", "exp")) {
  
  repelType = match.arg(repelType)
  rbf = match.arg(rbf)
  
  if(repelType == "none") {
    kern = function(d) {
      rep(1, length(d))
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








