



# makes test plots for simulated well data from sim study
testSuccessiveWellDat = function(i=1, adaptScen=c("batch", "adaptPref", "adaptVar"), 
                                 verbose=FALSE, nWells=60, batchSize=nWells/3) {
  adaptScen = match.arg(adaptScen)
  regenData = TRUE
  
  if(verbose) {
    print(paste0("generating well data for i: ", i, ", adapt scenario: ", adaptScen))
  }
  
  adaptScenCap = str_to_title(adaptScen)
  inputListFile = paste0("savedOutput/simStudy/simParList", adaptScenCap, ".RData")
  
  # parameters
  out = load(inputListFile)
  thisPar = wellDatCombsList[[i]]
  propVarCase = thisPar$propVarCase
  prefPar = thisPar$prefPar
  parI = sampleParI = thisPar$sampleParI
  wellDatI = thisPar$wellDatI
  repI = thisPar$repI
  sigmaSqErr = thisPar$nuggetVar
  sigmaSq = thisPar$sigmaSq
  repelAmount = thisPar$repEffect
  # nWells = thisPar$n
  if(adaptScen != "batch") {
    modelFitter = getFitModFuns()[[thisPar$sampModFunI]]
    isWatson = thisPar$sampModFunI == 4
  }
  repelDist = repAreaToDist(thisPar$repelAreaProp)
  seed = thisPar$seed
  
  # if the well data already exists and we don't want to regenerate it, don't
  wellDatFile = paste0("savedOutput/simStudy/wellDat/wellDat_", adaptScen, "_par", sampleParI, "_rep", repI, ".RData")
  if(file.exists(wellDatFile) && !regenData) {
    return(invisible(NULL))
  }
  
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
  
  # set repulsion parameters
  repelType = ifelse(repelDist == 0, "none", "rbf")
  if(repelDist != 0) {
    bwRepel = repelDist
  } else {
    bwRepel = NULL
  }
  
  if(adaptScen != "batch") {
    # sample the wells
    # wellDat = wellSampler(truthDat, seismicDat, modelFitter, nWells=nWells, minN=4,
    #                       predGrid=cbind(east=seismicDat$east, north=seismicDat$north),
    #                       transform=logit, invTransform=expit, prefPar=prefPar,
    #                       samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr,
    #                       repelType=repelType, bwRepel=bwRepel,
    #                       repelAmount=repelAmount, seed=seed,
    #                       int.strategy="eb", strategy="gaussian")
    
    # for testing purposes:
    out = wellSampler(truthDat=truthDat, seismicDat=seismicDat, 
                      modelFitter=modelFitter, nWells=nWells, minN=4,
                      predGrid=cbind(east=seismicDat$east, north=seismicDat$north),
                      transform=logit, invTransform=expit, prefPar=prefPar,
                      samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr,
                      repelType=repelType, bwRepel=bwRepel, isWatson=isWatson, 
                      repelAmount=repelAmount, seed=seed, verbose=verbose, 
                      int.strategy="eb", strategy="gaussian", 
                      batchSize=batchSize, getProbsNoRepOnly=FALSE)
    
    # profvis(wellDat <- wellSampler(truthDat, seismicDat, modelFitter, nWells=5, minN=4,
    #                                predGrid=cbind(east=seismicDat$east, north=seismicDat$north),
    #                                transform=logit, invTransform=expit, prefPar=prefPar,
    #                                samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr,
    #                                repelType=repelType, bwRepel=bwRepel,
    #                                repelAmount=repelAmount, seed=seed,
    #                                int.strategy="eb", strategy="gaussian"))
  } else {
    
    if(propVarCase == "realistic") {
      # indep
      otherRepI = (repI %% 100) + 1
      out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", otherRepI, ".txt"), force01=TRUE)
      indepDat = out$surfFrame
      
      # subsample
      indepDat = indepDat[goodCoords,]
      
      # standardize seismic, truth, and indep data on logit scale
      seismicDatStd = seismicDat
      truthDatStd = truthDat
      indepDatStd = indepDat
      seismicDatStd[,3] = logit(seismicDatStd[,3])
      truthDatStd[,3] = logit(truthDatStd[,3])
      indepDatStd[,3] = logit(indepDatStd[,3])
      seismicDatStd[,3] = (seismicDat[,3] - mean(seismicDat[,3]))/sd(seismicDat[,3])
      truthDatStd[,3] = (truthDat[,3] - mean(truthDat[,3]))/sd(truthDat[,3])
      indepDatStd[,3] = (indepDat[,3] - mean(indepDat[,3]))/sd(indepDat[,3])
      
      # combine them into a realistic mix and convert back to [0,1] scale
      sampleDat = seismicDat
      sampleDat[,3] = 0.5 * seismicDatStd[,3] + 0.25 * truthDatStd[,3] + 0.25 * indepDatStd[,3]
    } else if(propVarCase == "uniform") {
      # in this case, seismic data doesn't matter that much
      sampleDat = seismicDat
      sampleDat[,3] = rep(0, nrow(sampleDat))
    } else {
      stop("unrecognized propVarCase")
    }
    
    sampleDat[,3] = expit(sqrt(sigmaSq) * sampleDat[,3])
    
    # do batch sampling
    out = wellSampler(nWells=nWells, wellDat=NULL, seismicDat=seismicDat, 
                      predGrid=as.matrix(sampleDat[,1:2]), truthDat=truthDat, 
                      preds=sampleDat[,3], 
                      transform=logit, invTransform=expit, prefPar=prefPar, 
                      samplingModel=c("ipp"), sigmaSqErr=sigmaSqErr, 
                      repelType=repelType, bwRepel=bwRepel, 
                      rbf="uniform", repelAmount=Inf, verbose=verbose, 
                      seed=seed, int.strategy="eb", strategy="gaussian", 
                      batchSize=batchSize, getProbsNoRepOnly=FALSE)
    
    nBatches=3
    out$allPreds = matrix(rep(sampleDat[,3], nBatches), ncol=nBatches)
    
  }
  
  wellDat = out$wellDat
  allProbs = out$allProbs
  allLogitProbsNoRep = out$allLogitProbsNoRep
  allPreds = out$allPreds
  nBatches = ncol(allPreds)
  
  # plot things
  
  # plotting setup
  simStudyXlims = c(-12.5, 15012.5)
  simStudyYlims = c(-8.3335, 5008.4336)
  
  batchCols = rainbow(nBatches)
  
  gEast = seismicDat$east
  gNorth = seismicDat$north
  gSeismic = seismicDat$seismicEst
  gTruth = truthDat[,3]
  pEast = wellDat$east
  pNorth = wellDat$north
  pVolFrac = wellDat$volFrac
  
  eastGrid = sort(unique(gEast))
  northGrid = sort(unique(gNorth))
  
  seqCols = function(n) {purpleYellowSeqCols(n)}
  
  thisFigDir = paste0("figures/testSimStudy/")
  
  ticks = seq(0, 1, by=.2)
  tickLabs = as.character(ticks)
  
  # naive estimators ----
  pSeismic = bilinearInterp(cbind(pEast, pNorth), seismicDat, 
                            transform=logit, invTransform = expit)
  pTruth = bilinearInterp(cbind(pEast, pNorth), cbind(gEast, gNorth, gTruth), 
                          transform=logit, invTransform = expit)
  
  # plot one batch at a time with different colors
  
  pdf(file=paste0(thisFigDir, "simWellTruth_par", parI, "_rep", repI, ".pdf"), width=5, height=5)
  # par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
  squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
         zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
         asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
  for(i in 1:nBatches) {
    iStart = (i-1)*batchSize + 1
    iEnd = min(c(iStart + batchSize - 1, nrow(wellDat)))
    thisBatchIs = iStart:iEnd
    
    text(pEast[thisBatchIs], pNorth[thisBatchIs], labels=i, cex=.5)
  }
  dev.off()
  
  pdf(file=paste0(thisFigDir, "simWellSeismic_par", parI, "_rep", repI, ".pdf"), width=5, height=5)
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
  squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
         xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
         zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
         asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
  for(i in 1:nBatches) {
    iStart = (i-1)*batchSize + 1
    iEnd = min(c(iStart + batchSize - 1, nrow(wellDat)))
    thisBatchIs = iStart:iEnd
    
    text(pEast[thisBatchIs], pNorth[thisBatchIs], labels=i, cex=.5)
  }
  dev.off()
  
  pdf(file=paste0(thisFigDir, "simWellWell_par", parI, "_rep", repI, ".pdf"), width=5, height=5)
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
  splot(pEast, pNorth, pVolFrac, 
        xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
        zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
        asp=1, smallplot=c(.83,.87,.25,.8), cex=.3, lwd=.25)
  dev.off()
  
  pdf(file=paste0(thisFigDir, "simWellWellTruth_par", parI, "_rep", repI, ".pdf"), width=5, height=5)
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
  splot(pEast, pNorth, pTruth, 
        xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
        zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
        asp=1, smallplot=c(.83,.87,.25,.8), cex=.3, lwd=.25)
  dev.off()
  
  pdf(file=paste0(thisFigDir, "simWellWellSeismic_par", parI, "_rep", repI, ".pdf"), width=5, height=5)
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
  splot(pEast, pNorth, pSeismic, 
        xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
        zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
        asp=1, smallplot=c(.83,.87,.25,.8), cex=.3, lwd=.25)
  dev.off()
  
  # plot preds/probs
  if(adaptScen != "batch") {
    pdf(file=paste0(thisFigDir, "simWellProbs_par", parI, "_rep", repI, ".pdf"), width=8, height=8)
    par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
    
    # truth
    squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
           zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
           asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
    for(i in 1:nBatches) {
      iStart = (i-1)*batchSize + 1
      iEnd = min(c(iStart + batchSize - 1, nrow(wellDat)))
      thisBatchIs = iStart:iEnd
      
      text(pEast[thisBatchIs], pNorth[thisBatchIs], labels=i, cex=.5)
    }
    
    # plot each batch against probs for that batch (noting that repulsion updates)
    for(i in 1:nBatches) {
      squilt(gEast, gNorth, allProbs[,i], grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
             xlab="Easting", ylab="Northing", main=paste0("Draw probs batch ", i), 
             asp=1, smallplot=c(.83,.87,.25,.8))
      
      iStart = (i-1)*batchSize + 1
      iEnd = min(c(iStart + batchSize - 1, nrow(wellDat)))
      thisBatchIs = iStart:iEnd
      
      text(pEast[thisBatchIs], pNorth[thisBatchIs], labels=1:length(thisBatchIs), cex=.5)
    }
    
    dev.off()
    
    # plot preds/probs
    pdf(file=paste0(thisFigDir, "simWellLogitProbsNoRep_par", parI, "_rep", repI, ".pdf"), width=8, height=8)
    par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
    
    # truth
    squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
           zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
           asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
    for(i in 1:nBatches) {
      iStart = (i-1)*batchSize + 1
      iEnd = min(c(iStart + batchSize - 1, nrow(wellDat)))
      thisBatchIs = iStart:iEnd
      
      text(pEast[thisBatchIs], pNorth[thisBatchIs], labels=i, cex=.5)
    }
    
    # plot each batch against probs for that batch (noting that repulsion updates)
    for(i in 1:nBatches) {
      squilt(gEast, gNorth, allLogitProbsNoRep[,i], grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
             xlab="Easting", ylab="Northing", main=paste0("Logit probs no rep batch ", i), 
             asp=1, smallplot=c(.83,.87,.25,.8))
      
      iStart = (i-1)*batchSize + 1
      iEnd = min(c(iStart + batchSize - 1, nrow(wellDat)))
      thisBatchIs = iStart:iEnd
      
      text(pEast[thisBatchIs], pNorth[thisBatchIs], labels=1:length(thisBatchIs), cex=.5)
    }
    
    dev.off()
  }
  
  pdf(file=paste0(thisFigDir, "simWellPreds_par", parI, "_rep", repI, ".pdf"), width=8, height=8)
  par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
  
  # truth
  squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
         zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
         asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
  for(i in 1:nBatches) {
    iStart = (i-1)*batchSize + 1
    iEnd = min(c(iStart + batchSize - 1, nrow(wellDat)))
    thisBatchIs = iStart:iEnd
    
    text(pEast[thisBatchIs], pNorth[thisBatchIs], labels=i, cex=.5)
  }
  
  # plot each batch against probs for that batch (noting that repulsion updates)
  for(i in 1:nBatches) {
    squilt(gEast, gNorth, allPreds[,i], grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
           zlim=c(0, 1), xlab="Easting", ylab="Northing", main=paste0("Draw preds batch ", i), 
           asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
    
    iStart = (i-1)*batchSize + 1
    iEnd = min(c(iStart + batchSize - 1, nrow(wellDat)))
    thisBatchIs = iStart:iEnd
    
    text(pEast[thisBatchIs], pNorth[thisBatchIs], labels=1:length(thisBatchIs), cex=.5)
  }
  
  dev.off()
  
  browser()
  
  
  
  
}

# for batch uniform repelAreaProp==.01 scenario, make sure repel distances apply
testRepelDists = function() {
  
  inputListFile = paste0("savedOutput/simStudy/simParListBatch.RData")
  
  out = load(inputListFile)
  
  parI = 6
  sampleParCombs[parI,]
  
  
  wellIs = wellDatCombs$wellDatI[wellDatCombs$sampleParI == parI]
  
  repelDist = repAreaToDist(sampleParCombs$repelAreaProp[parI])
  
  minDs = numeric(length(wellIs))
  for(i in 1:length(wellIs)) {
    print(paste0("calculating minum distance for i = ", i, "/", length(wellIs)))
    
    
    wellI = wellIs[i]
    repI = wellDatCombs$repI[wellI]
    
    wellFile = paste0("savedOutput/simStudy/wellDat/wellDat_batch_par", parI, "_rep", repI, ".RData")
    out = load(wellFile)
    
    dMat = rdist(as.matrix(wellDat[,1:2]))
    minDs[i] = min(dMat[dMat != 0])
    
    if(minDs[i] < repelDist) {
      print(paste0("WARNING: minD for i=", i, " is ", minDs[i], " < ", repelDist))
    }
  }
  
  browser()
  
  
  
}











