# make a single simple set of datasets for testing purposes


# NOTE: this is not supposed to look like the simulation study data...
getTestDat = function(n=20, xlim=simStudyXlims, ylim=simStudyYlims, sigmaSqTruth=1, 
                      corrSeismic=.6, sigmaSqErr=.1, aRange=diff(xlim)/4, 
                      seismicPrefPar=.5, 
                      grid=make.surface.grid(list(east=seq(xlim[1], xlim[2], l=100), 
                                                  north=seq(ylim[1], ylim[2], l=100))), 
                      seed=123) {
  require(fields)
  set.seed(seed)
  
  # simulate spatial fields for calculating truth and seismic data
  Sigma = stationary.cov(grid, aRange=aRange)*sigmaSqSpat
  L = t(chol(Sigma))
  truth = L %*% rnorm(nrow(L))
  indepPt = L %*% rnorm(nrow(L))
  seismicEst = corrSeismic^2 * truth + sqrt(1-corrSeismic^2) * indepPt
  
  # sample well point locations at grid centroids based on seismic estimates
  lambdas = exp(seismicPrefPar * seismicEst)
  sampleProbs = lambdas * (1/sum(lambdas))
  sampleI = sample(1:length(lambdas), size = n, replace=T, prob=sampleProbs)
  
  # now sample well point locations uniformly within the cells
  gridWidth = diff(sort(grid[1,])[1:2])
  pts = grid[sampleI,] + matrix(runif(2*n)*gridWidth - gridWidth/2, ncol=2)
  wellSeismicEst = bilinearInterp(pts, cbind(grid, seismicEst))
  
  # transform the truth and seismic data to correspond to fractions
  truthTrans = expit(truth)
  seismicEstTrans = expit(seismicEst)
  wellSeismicEstTrans = expit(wellSeismicEst)
  
  # construct final datasets and truth
  truthDat = cbind(grid, truth=truthTrans)
  seismicDat = cbind(grid, seismicEst=seismicEstTrans)
  wellDat = data.frame(east=pts[,1], north=pts[,2], volFrac=wellSeismicEstTrans)
  
  list(wellDat=wellDat, seismicDat=seismicDat, truthDat=truthDat)
}