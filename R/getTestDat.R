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
                                seed=1, saveDat=TRUE, minN=4, batchSize=1, ...) {
  
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
                            repelAmount=repelAmount, seed=NULL, batchSize=batchSize, ...)
  
  # save datasets if need be
  if(saveDat) {
    if(repelType == "none") {
      repelStr = ""
    } else {
      repelStr = paste0("_", rbf)
    }
    batchStr = ""
    if(batchSize != 1) {
      batchStr = paste0("batch", batchSize)
    }
    save(truthTestDat, seismicTestDat, wellTestDat, file=paste0(globalDir, "/testSuccessiveDat_", samplingModel, repelStr, batchStr, ".RData"))
  }
  
  invisible(list(wellDat=wellTestDat, seismicDat=seismicTestDat, truthDat=truthTestDat))
}








