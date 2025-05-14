# script for testing Successive SPDE model

# _ ----
# load seismic pref data ----
# _ ----
gEast = seismicTestDat$east
gNorth = seismicTestDat$north
gSeismic = seismicTestDat$seismicEst
gTruth = truthTestDat$truth
pEast = wellTestDat$east
pNorth = wellTestDat$north
pVolFrac = wellTestDat$volFrac

eastGrid = sort(unique(gEast))
northGrid = sort(unique(gNorth))

print(cor(gSeismic, gTruth))

seqCols = function(n) {purpleYellowSeqCols(n)}

# fit model (pproc) ----
minN=4
out = fitSuccessiveSPDEsimDat(wellTestDat, seismicTestDat, 
                              successiveMethod="pproc", includeInt=TRUE, 
                              minN=minN)

summary(out$mod)

names(out)

# plot predictions ----

preds = out$predEst
predAggMat = out$predAggMat
predAgg = out$predAggEst
predAggU = out$predAggUpper
predAggL = out$predAggLower
trueAgg = mean(gTruth)
predAgg
trueAgg

allWs = out$allWs

pdf(file=paste0(figDir, "testing/testPreds_SuccessiveSPDE_pproc.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Successive SPDE, pproc)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SuccessiveSPDE_pproc.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (Successive SPDE, pproc)", 
     xlab="Volume Fraction", freq=F, xlim=range(c(predAggMat, trueAgg - .01, trueAgg + .01)))
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
legend("topright", legend=c("Mean", "80% CI", "Truth"), lty=c(1, 2, 1), col=c("purple", "purple", "black"))
dev.off()

pdf(file=paste0(figDir, "testing/testWs_SuccessiveSPDE_pproc.pdf"), width=8, height=12)
par(mfrow=c(3,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Successive SPDE, pproc)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, expit(allWs[,5-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W5) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:5], pNorth[1:5], cex=.5)

squilt(gEast, gNorth, expit(allWs[,10-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W10) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:10], pNorth[1:10], cex=.5)

squilt(gEast, gNorth, expit(allWs[,ncol(allWs)]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main=paste0("expit(W", ncol(allWs)+minN-1, ") (kde)"), 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()








# fit model (SPDE) ----
minN=4
out = fitSuccessiveSPDEsimDat(wellTestDat, seismicTestDat, 
                              successiveMethod="spde", includeInt=FALSE, 
                              minN=minN)

summary(out$mod)

names(out)

# plot predictions ----

preds = out$predEst
predAggMat = out$predAggMat
predAgg = out$predAggEst
predAggU = out$predAggUpper
predAggL = out$predAggLower
trueAgg = mean(gTruth)
predAgg
trueAgg

allWs = out$allWs

pdf(file=paste0(figDir, "testing/testPreds_SuccessiveSPDE_spde.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Successive SPDE, pproc)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SuccessiveSPDE_spde.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (Successive SPDE, pproc)", 
     xlab="Volume Fraction", freq=F, xlim=range(c(predAggMat, trueAgg - .01, trueAgg + .01)))
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
legend("topright", legend=c("Mean", "80% CI", "Truth"), lty=c(1, 2, 1), col=c("purple", "purple", "black"))
dev.off()

pdf(file=paste0(figDir, "testing/testWs_SuccessiveSPDE_spde.pdf"), width=8, height=12)
par(mfrow=c(3,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Successive SPDE, pproc)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, expit(allWs[,5-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W5) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:5], pNorth[1:5], cex=.5)

squilt(gEast, gNorth, expit(allWs[,10-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W10) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:10], pNorth[1:10], cex=.5)

squilt(gEast, gNorth, expit(allWs[,ncol(allWs)]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main=paste0("expit(W", ncol(allWs)+minN-1, ") (kde)"), 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

# _ ----
# load truth pref data ----
# _ ----
gEast = seismicTestDat_truthPref$east
gNorth = seismicTestDat_truthPref$north
gSeismic = seismicTestDat_truthPref$seismicEst
gTruth = truthTestDat_truthPref$truth
pEast = wellTestDat_truthPref$east
pNorth = wellTestDat_truthPref$north
pVolFrac = wellTestDat_truthPref$volFrac

eastGrid = sort(unique(gEast))
northGrid = sort(unique(gNorth))

# fit model (pproc) ----
minN=4
out = fitSuccessiveSPDEsimDat(wellTestDat_truthPref, seismicTestDat_truthPref, 
                              successiveMethod="pproc", includeInt=TRUE, 
                              minN=minN)

summary(out$mod)

names(out)

# plot predictions ----

preds = out$predEst
predAggMat = out$predAggMat
predAgg = out$predAggEst
predAggU = out$predAggUpper
predAggL = out$predAggLower
trueAgg = mean(gTruth)
predAgg
trueAgg

allWs = out$allWs

pdf(file=paste0(figDir, "testing/testPreds_SuccessiveSPDE_pproc_truthPref.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Successive SPDE, pproc)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SuccessiveSPDE_pproc_truthPref.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (Successive SPDE, pproc)", 
     xlab="Volume Fraction", freq=F)
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
legend("topright", legend=c("Mean", "80% CI", "Truth"), lty=c(1, 2, 1), col=c("purple", "purple", "black"))
dev.off()

pdf(file=paste0(figDir, "testing/testWs_SuccessiveSPDE_pproc_truthPref.pdf"), width=8, height=12)
par(mfrow=c(3,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Successive SPDE, pproc)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, expit(allWs[,5-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W5) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:5], pNorth[1:5], cex=.5)

squilt(gEast, gNorth, expit(allWs[,10-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W10) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:10], pNorth[1:10], cex=.5)

squilt(gEast, gNorth, expit(allWs[,ncol(allWs)]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main=paste0("expit(W", ncol(allWs)+minN-1, ") (kde)"), 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()








# fit model (SPDE) ----
minN=4
out = fitSuccessiveSPDEsimDat(wellTestDat, seismicTestDat, 
                              successiveMethod="spde", includeInt=FALSE, 
                              minN=minN)

summary(out$mod)

names(out)

# plot predictions ----

preds = out$predEst
predAggMat = out$predAggMat
predAgg = out$predAggEst
predAggU = out$predAggUpper
predAggL = out$predAggLower
trueAgg = mean(gTruth)
predAgg
trueAgg

allWs = out$allWs

pdf(file=paste0(figDir, "testing/testPreds_SuccessiveSPDE_spde_truthPref.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Successive SPDE, pproc)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SuccessiveSPDE_spde_truthPref.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (Successive SPDE, pproc)", 
     xlab="Volume Fraction", freq=F)
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
legend("topright", legend=c("Mean", "80% CI", "Truth"), lty=c(1, 2, 1), col=c("purple", "purple", "black"))
dev.off()

pdf(file=paste0(figDir, "testing/testWs_SuccessiveSPDE_spde_truthPref.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, expit(allWs[,5-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W5) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:5], pNorth[1:5], cex=.5)

squilt(gEast, gNorth, expit(allWs[,10-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W10) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:5], pNorth[1:5], cex=.5)

squilt(gEast, gNorth, expit(allWs[,15-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W15) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:5], pNorth[1:5], cex=.5)

squilt(gEast, gNorth, expit(allWs[,20-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W20) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:5], pNorth[1:5], cex=.5)
dev.off()

# _ ----
# load successive pref data ----
# _ ----

gEast = seismicTestDat_successiveIPP$east
gNorth = seismicTestDat_successiveIPP$north
gSeismic = seismicTestDat_successiveIPP$seismicEst
gTruth = truthTestDat_successiveIPP$truth
pEast = wellTestDat_successiveIPP$east
pNorth = wellTestDat_successiveIPP$north
pVolFrac = wellTestDat_successiveIPP$volFrac

eastGrid = sort(unique(gEast))
northGrid = sort(unique(gNorth))

print(cor(gSeismic, gTruth))

seqCols = function(n) {purpleYellowSeqCols(n)}

# fit model (pproc) ----
minN=4
out = fitSuccessiveSPDEsimDat(wellTestDat_successiveIPP, seismicTestDat_successiveIPP, 
                              successiveMethod="pproc", includeInt=TRUE, 
                              minN=minN)

summary(out$mod)

out$fixedEffectSummary
out$parameterSummaryTable

names(out)

out$timings
# totalTime modelDefineTime modelFitTime posteriorSamplingTime otherTime modelDefinePct modelFitTimePct
# elapsed      44.2            1.07        10.42                 26.67      6.04     0.02420814       0.2357466
# posteriorSamplingTimePct otherTimePct
# elapsed                0.6033937    0.1366516

# plot predictions ----

preds = out$predEst
obsPreds = out$obsEst
predAggMat = out$predAggMat
predAgg = out$predAggEst
predAggU = out$predAggUpper
predAggL = out$predAggLower
trueAgg = mean(gTruth)
predAgg
trueAgg

allWs = out$allWs

pdf(file=paste0(figDir, "testing/testPreds_SuccessiveSPDE_pproc_succesiveIPP.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Successive SPDE, pproc)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SuccessiveSPDE_pproc_succesiveIPP.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (Successive SPDE, pproc)", 
     xlab="Volume Fraction", freq=F)
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
legend("topright", legend=c("Mean", "80% CI", "Truth"), lty=c(1, 2, 1), col=c("purple", "purple", "black"))
dev.off()

pdf(file=paste0(figDir, "testing/testWs_SuccessiveSPDE_pproc_succesiveIPP.pdf"), width=8, height=12)
par(mfrow=c(3,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Successive SPDE, pproc)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, expit(allWs[,5-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W5) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:5], pNorth[1:5], cex=.5)

squilt(gEast, gNorth, expit(allWs[,10-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W10) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:10], pNorth[1:10], cex=.5)

squilt(gEast, gNorth, expit(allWs[,ncol(allWs)]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main=paste0("expit(W", ncol(allWs)+minN-1, ") (kde)"), 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

print(paste0("MSE (seismic and SPDE): ", mean(gTruth - gSeismic)^2, " and ", mean(gTruth - preds)^2))
print(paste0("logit MSE (seismic and SPDE): ", mean(logit(gTruth) - logit(gSeismic))^2, " and ", mean(logit(gTruth) - logit(preds))^2))

print(paste0("MSE (seismic and SPDE): ", mean(pTruth - pSeismic)^2, " and ", mean(pTruth - obsPreds)^2))
print(paste0("logit MSE (seismic and SPDE): ", mean(logit(pTruth) - logit(pSeismic))^2, " and ", mean(logit(pTruth) - logit(obsPreds))^2))






# fit model (SPDE) ----
minN=4
out = fitSuccessiveSPDEsimDat(wellTestDat_successiveIPP, seismicTestDat_successiveIPP, 
                              successiveMethod="spde", includeInt=FALSE, 
                              minN=minN)

out$fixedEffectSummary
#          mean        sd   0.5quant   0.1quant   0.9quant
# X1 -0.5962212 0.2711696 -0.5963847 -0.9418894 -0.2503176 # w
# X2  0.6452703 0.3387565  0.6401297  0.2190478  1.0790130 # seismic

out$parameterSummaryTable
#                       Est           SD       Qlower          Q50       Qupper
# totalVar     3.813636e-01 7.506408e-02 2.774713e-01 3.692482e-01 4.919691e-01
# spatialVar   3.711340e-01 7.514567e-02 2.668994e-01 3.592318e-01 4.835769e-01
# errorVar     1.022956e-02 1.816972e-03 8.363372e-03 9.935513e-03 1.215749e-02
# totalSD      6.145725e-01 6.056310e-02 5.267554e-01 6.076580e-01 7.014051e-01
# spatialSD    6.061022e-01 6.146500e-02 5.166230e-01 5.993595e-01 6.953970e-01
# errorSD      1.007494e-01 8.899377e-03 9.145147e-02 9.967704e-02 1.102610e-01
# spatialRange 2.816145e+03 7.436888e+02 1.935359e+03 2.741885e+03 3.604644e+03

names(out)

out$timings
# NOTE: these timings are inaccurate. They are only for the final SPDE fit:
# totalTime modelDefineTime modelFitTime posteriorSamplingTime otherTime modelDefinePct modelFitTimePct
# elapsed     52.36             0.5        13.78                 31.66      6.42    0.009549274        0.263178
# posteriorSamplingTimePct otherTimePct
# elapsed                  0.60466    0.1226127

# plot predictions ----

preds = out$predEst
obsPreds = out$obsEst
predAggMat = out$predAggMat
predAgg = out$predAggEst
predAggU = out$predAggUpper
predAggL = out$predAggLower
trueAgg = mean(gTruth)
predAgg
trueAgg

allWs = out$allWs

pdf(file=paste0(figDir, "testing/testPreds_SuccessiveSPDE_spde_noInt_succesiveIPP.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Successive SPDE, pproc)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SuccessiveSPDE_spde_noInt_succesiveIPP.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (Successive SPDE, pproc)", 
     xlab="Volume Fraction", freq=F)
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
legend("topright", legend=c("Mean", "80% CI", "Truth"), lty=c(1, 2, 1), col=c("purple", "purple", "black"))
dev.off()

pdf(file=paste0(figDir, "testing/testWs_SuccessiveSPDE_spde_noInt_succesiveIPP.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, expit(allWs[,5-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W5) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:5], pNorth[1:5], cex=.5)

squilt(gEast, gNorth, expit(allWs[,10-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W10) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:5], pNorth[1:5], cex=.5)

squilt(gEast, gNorth, expit(allWs[,15-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W15) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:5], pNorth[1:5], cex=.5)

squilt(gEast, gNorth, expit(allWs[,20-minN+1]), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="expit(W20) (kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast[1:5], pNorth[1:5], cex=.5)
dev.off()

print(paste0("MSE (seismic and SPDE): ", mean(gTruth - gSeismic)^2, " and ", mean(gTruth - preds)^2))
# [1] "MSE (seismic and SPDE): 0.000575049976194529 and 0.00686684525189851"
print(paste0("logit MSE (seismic and SPDE): ", mean(logit(gTruth) - logit(gSeismic))^2, " and ", mean(logit(gTruth) - logit(preds))^2))
# [1] "logit MSE (seismic and SPDE): 0.00918437079448658 and 0.129870687806384"

print(paste0("MSE (seismic and SPDE): ", mean(pTruth - pSeismic)^2, " and ", mean(pTruth - obsPreds)^2))
# [1] "MSE (seismic and SPDE): 9.36241385035057e-07 and 5.8214704528366e-06"
print(paste0("logit MSE (seismic and SPDE): ", mean(logit(pTruth) - logit(pSeismic))^2, " and ", mean(logit(pTruth) - logit(obsPreds))^2))
# [1] "logit MSE (seismic and SPDE): 0.000706596851712083 and 0.000151513199464639"
