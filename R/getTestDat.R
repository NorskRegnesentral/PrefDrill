# make a single simple set of datasets for testing purposes


# NOTE: this is not supposed to look like the simulation study data...
getTestDat = function(n=20, xlim=simStudyXlims, ylim=simStudyYlims, 
                      muSpat=0, sigmaSqSpat=1, corrSeismic=.6, sigmaSqErr=.1, 
                      aRange=diff(xlim)/4, seismicPrefPar=2.5, 
                      grid=make.surface.grid(list(east=seq(xlim[1], xlim[2], l=50), 
                                                  north=seq(ylim[1], ylim[2], l=50))), 
                      seed=123, 
                      saveDat=TRUE) {
  require(fields)
  set.seed(seed)
  
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
  lambdas = exp(seismicPrefPar * seismicEst)
  sampleProbs = lambdas * (1/sum(lambdas))
  sampleI = sample(1:length(lambdas), size = n, replace=T, prob=sampleProbs)
  
  # now sample well point locations uniformly within the cells
  gridWidth = diff(sort(grid[1,])[1:2])
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
  
  browser() # run checks to see if it looks reasonable
  if(saveDat) {
    save(truthTestDat, seismicTestDat, wellTestDat, file=paste0(globalDir, "/testDat.RData"))
  }
  
  list(wellDat=wellTestDat, seismicDat=seismicTestDat, truthDat=truthTestDat)
}








