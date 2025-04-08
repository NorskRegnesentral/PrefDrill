# make a single simple set of datasets for testing purposes


# NOTE: this is not supposed to look like the simulation study data...
getTestDat = function(n=20, xlim=simStudyXlims, ylim=simStudyYlims, 
                      muSpat=0, sigmaSqSpat=1, corrSeismic=.6, sigmaSqErr=.01, 
                      aRange=diff(xlim)/4, prefPar=2.5, gridWidth=diff(xlim)/50, 
                      grid=make.surface.grid(list(east=seq(xlim[1]+gridWidth/2, xlim[2]-gridWidth/2, l=50), 
                                                  north=seq(ylim[1]+gridWidth/2, ylim[2]-gridWidth/2, l=50))), 
                      prefType=c("seismic", "truth"), seed=1, saveDat=TRUE) {
  require(fields)
  set.seed(seed)
  prefType = match.arg(prefType)
  
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
getSuccessiveTestDat = function(n=20, xlim=simStudyXlims, ylim=simStudyYlims, 
                                muSpat=0, sigmaSqSpat=1, corrSeismic=.6, sigmaSqErr=.01, 
                                modelFitter=fitSPDEsimDat, aRange=diff(xlim)/4, prefPar=2.5, 
                                gridWidth=diff(xlim)/50, transform=logit, invTransform=expit, 
                                grid=make.surface.grid(list(east=seq(xlim[1]+gridWidth/2, xlim[2]-gridWidth/2, l=50), 
                                                            north=seq(ylim[1]+gridWidth/2, ylim[2]-gridWidth/2, l=50))), 
                                samplingModel=c("ipp"), seed=1, saveDat=TRUE, minN=4, ...) {
  
  require(fields)
  set.seed(seed)
  samplingModel = match.arg(samplingModel)
  
  # simulate spatial fields for calculating truth and seismic data
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
  
  # construct well data via successive sampling process
  wellTestDat = wellSampler(truthDat=truthTestDat, seismicDat=seismicTestDat, 
                            modelFitter=modelFitter, nWells=n, minN=minN, 
                            predGrid=grid, prefPar=prefPar, 
                            transform=transform, invTransform=invTransform, 
                            samplingModel=samplingModel, sigmaSqErr=sigmaSqErr, seed=NULL, ...)
  
  # save datasets if need be
  if(saveDat) {
    save(truthTestDat, seismicTestDat, wellTestDat, file=paste0(globalDir, "/testSuccessiveDat_", samplingModel, ".RData"))
  }
  
  invisible(list(wellDat=wellTestDat, seismicDat=seismicTestDat, truthDat=truthTestDat))
}

wellSampler = function(truthDat, seismicDat, modelFitter, nWells=20, minN=4, 
                       predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
                       transform=logit, invTransform=expit, prefPar=0, 
                       samplingModel=c("ipp"), sigmaSqErr=.1^2, seed=NULL, ...) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  samplingModel = match.arg(samplingModel)
  
  # get first datapoints based on only seismic data
  print("Initial well sampling...")
  initWellDat = basicWellSampler(nWells=minN, wellDat=NULL, seismicDat=seismicDat, 
                                 truthDat=truthDat, modelFitter=modelFitter, 
                                 predGrid=predGrid, prefPar=prefPar, 
                                 transform=transform, invTransform=invTransform, 
                                 samplingModel=samplingModel, sigmaSqErr=sigmaSqErr, 
                                 seed=NULL, ...)
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
                                  seed=NULL)
    
    wellDat = rbind(wellDat, newWellDat)
  }
  
  wellDat
}

basicWellSampler = function(nWells=1, wellDat=NULL, seismicDat, truthDat=NULL, modelFitter, 
                            predGrid=cbind(east=seismicDat$east, north=seismicDat$north), 
                            transform=logit, invTransform=expit, prefPar=0, 
                            samplingModel=c("ipp"), sigmaSqErr=.1^2, seed=NULL, 
                            ...) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  samplingModel = match.arg(samplingModel)
  
  # first fit the model and get predictions over the grid (or just take seismic estimates)
  if(!is.null(wellDat)) {
    mod = modelFitter(wellDat=wellDat, seismicDat=seismicDat, predGrid=predGrid, 
                      transform=transform, invTransform=invTransform, ...)
    preds = mod$predEst
  } else {
    preds = seismicDat[,3]
  }
  
  if(any(preds <= 0)) {
    stop("some sampling probabilities are negative. Check data and predictions on correct scale.")
  }
  
  # pick the next grid point based on the sampling model
  if(samplingModel == "ipp") {
    # in the inhomogeneous Poisson process, sampling probability is proportional 
    # to the predictions
    probs = preds^prefPar
    probs = probs * (1 / sum(probs))
  } else {
    stop(paste0("sampling model ", samplingModel, " not yet supported"))
  }
  
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
  
  res
}






