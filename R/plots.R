# functions for making plots and illustrations for the manuscript


plotRepulsion = function(seed=123, regenData=FALSE, thisN = c(20, 60, 250), prefPar = c(.5, .75)) {
  thisN = matchArgNumeric(thisN)
  prefPar = matchArgNumeric(prefPar)
  
  # get sampling probabilities and well data locations under high preferentiality, 
  # varying repulsion
  out = load("savedOutput/simStudy/simParListBatch.RData")
  
  
  thisRepI = 2
  
  # parameters
  propVarCase = "diggle"
  
  repI = thisRepI
  sigmaSqErr = 0.01
  sigmaSq = 1
  repelAmount = Inf
  nWells = thisN
  
  # get truth and data
  
  # seismic data
  out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedPred_", repI, ".txt"), force01=TRUE)
  seismicDat = out$surfFrame
  
  # truth
  out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", repI, ".txt"), force01=TRUE)
  truthDat = out$surfFrame
  
  # subsample
  goodCoords = subsampleSimStudyGrid(seismicDat)
  seismicDat = seismicDat[goodCoords,]
  truthDat = truthDat[goodCoords,]
  
  # standardize (center + normalize)
  normDat = getNormFac(seismicDat=seismicDat, truthDat=truthDat, indepDat=truthDat, 
                       subsampled=TRUE, goodCoords=goodCoords)
  truthDatStd = normDat$truthDatStd
  
  # combine them into a realistic mix and convert back to [0,1] scale
  sampleDat = seismicDat
  
  sampleDat[,3] = truthDatStd[,3]
  sampleDat[,3] = exp(sampleDat[,3])
  
  # do batch sampling
  if(thisN == 20) {
    noRepelDist = repAreaToDist(0)
    midRepelDist = repAreaToDist(0.004) # 0.001 in new p_rep parameterization
    highRepelDist = repAreaToDist(0.02) # 0.005 in new p_rep parameterization
  } else if(thisN == 60) {
    noRepelDist = repAreaToDist(0)
    midRepelDist = repAreaToDist(0.004) # 0.001 in new p_rep parameterization
    highRepelDist = repAreaToDist(0.02) # 0.005 in new p_rep parameterization
  } else if (thisN == 250) {
    noRepelDist = repAreaToDist(0)
    midRepelDist = repAreaToDist(0.002)  # 0.001 in new p_rep parameterization
    highRepelDist = repAreaToDist(0.004) # 0.005 in new p_rep parameterization
  } else {
    stop()
  }
  
  
  wellDatFileName = paste0("savedOutput/illustrations/illustrateRepulsion_n", thisN, "_pref", prefPar, ".RData")
  
  if(!file.exists(wellDatFileName) || regenData) {
    # sample well data, making sure to use the same sampling probs and random seed, but 
    # different repulsion levels
    
    wellDatNoRepel = wellSampler(nWells=nWells, wellDat=NULL, seismicDat=seismicDat, 
                                 predGrid=as.matrix(sampleDat[,1:2]), truthDat=truthDat, 
                                 preds=sampleDat[,3], 
                                 transform=logit, invTransform=expit, prefPar=prefPar, 
                                 samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr, 
                                 repelType="none", bwRepel=noRepelDist, 
                                 rbf="uniform", repelAmount=Inf, batchSize=5, isWatson=FALSE, 
                                 seed=seed, int.strategy="eb", strategy="gaussian", getProbsNoRepOnly=FALSE)
    
    # wellDatMidRepel = wellSampler(nWells=nWells, wellDat=NULL, seismicDat=seismicDat, 
    #                               predGrid=as.matrix(sampleDat[,1:2]), truthDat=truthDat, 
    #                               preds=sampleDat[,3], 
    #                               transform=logit, invTransform=expit, prefPar=prefPar, 
    #                               samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr, 
    #                               repelType="rbf", bwRepel=midRepelDist, 
    #                               rbf="uniform", repelAmount=Inf, batchSize=5, isWatson=FALSE, 
    #                               seed=seed, int.strategy="eb", strategy="gaussian", getProbsNoRepOnly=FALSE)
    
    wellDatHighRepel = wellSampler(nWells=nWells, wellDat=NULL, seismicDat=seismicDat, 
                                   predGrid=as.matrix(sampleDat[,1:2]), truthDat=truthDat, 
                                   preds=sampleDat[,3], 
                                   transform=logit, invTransform=expit, prefPar=prefPar, 
                                   samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr, 
                                   repelType="rbf", bwRepel=highRepelDist, 
                                   rbf="uniform", repelAmount=Inf, batchSize=5, isWatson=FALSE, 
                                   seed=seed, int.strategy="eb", strategy="gaussian", getProbsNoRepOnly=FALSE)
    
    # save(wellDatNoRepel, wellDatMidRepel, wellDatHighRepel, file=wellDatFileName)
    save(wellDatNoRepel, wellDatHighRepel, file=wellDatFileName)
  } else {
    out = load(wellDatFileName)
  }
  
  # get well data and sampling probabilities
  noRepelProbs = exp(wellDatNoRepel$logProbsNoRep)
  noRepelProbs = noRepelProbs * (1/sum(noRepelProbs))
  noRepelWellDat = wellDatNoRepel$wellDat[1:thisN,]
  
  # midRepelProbs = wellDatMidRepel$allProbs[,1]
  # midRepelWellDat = wellDatMidRepel$wellDat[1:thisN,]
  
  highRepelProbs = wellDatHighRepel$allProbs[,ncol(wellDatHighRepel$allProbs)]
  highRepelWellDat = wellDatHighRepel$wellDat[1:thisN,]
  highRepelProbs = highRepelProbs * max(noRepelProbs)/max(highRepelProbs)
  
  # plot setup
  seqCols = function(n) {blueSeqCols(n, rev = TRUE)}
  gEast = seismicDat[,1]
  gNorth = seismicDat[,2]
  eastGrid = sort(unique(gEast))
  northGrid = sort(unique(gNorth))
  
  drawCircle <- function(x0, y0, r, n = 200, ...) {
    theta <- seq(0, 2*pi, length.out = n)
    lines(x0 + r * cos(theta), y0 + r * sin(theta), ...)
  }
  
  
  thisCex = ifelse(thisN == 250, .3, 1)
  
  pdf(file=paste0("figures/illustrations/showRepulsion_n", thisN, "_pref", prefPar, ".pdf"), width=12, height=3.4)
  par(mfrow=c(1, 2), oma=c( 0,2.5,0,3.6), mar=c(4.1, 1.1, 1.9, 1.7))
  
  xlab = ifelse(thisN == 250, "Easting", "")
  title1 = ifelse(thisN == 20, "No repulsion", "")
  title2 = ifelse(thisN == 20, "High repulsion", "")
  
  # top row
  squilt(gEast, gNorth, noRepelProbs, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
         xlab=xlab, ylab="", main=title1, scaleFun=log, scaleFunInverse=exp, 
         smallplot=c(.95,.98,.25,.87), addColorBar=FALSE, ticks=getLogScaleTicks(noRepelProbs), 
         tickLabels = as.character(getLogScaleTicks(noRepelProbs)))
  points(noRepelWellDat[,1], noRepelWellDat[,2], cex=thisCex, col="red", pch=19)
  
  mtext("Northing", side=2, line=1.5, outer=TRUE)
  
  
  squilt(gEast, gNorth, noRepelProbs, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
         xlab=xlab, ylab="", main=title2, scaleFun=log, scaleFunInverse=exp, 
         smallplot=c(.95,.98,.25,.87), ticks=getLogScaleTicks(noRepelProbs), 
         tickLabels = as.character(getLogScaleTicks(noRepelProbs)))
  points(highRepelWellDat[,1], highRepelWellDat[,2], cex=thisCex, col="red", pch=19)
  
  # apply(highRepelWellDat[,1:2], 1, function(p) drawCircle(p[1], p[2], r=highRepelDist, lty = 2, lwd=.5, col=rgb(1,0,0,.7)))
  
  dev.off()
  
}

makeAllRepulsionPlots = function(regenData=TRUE) {
  plotRepulsion(regenData=regenData, thisN=20)
  plotRepulsion(regenData=regenData, thisN=60)
  plotRepulsion(regenData=regenData, thisN=250)
  plotRepulsion(regenData=regenData, thisN=20, prefPar=.75)
  plotRepulsion(regenData=regenData, thisN=60, prefPar=.75)
  plotRepulsion(regenData=regenData, thisN=250, prefPar=.75)
}












