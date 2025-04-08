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
     xlab="Volume Fraction", freq=F)
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
     xlab="Volume Fraction", freq=F)
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
legend("topright", legend=c("Mean", "80% CI", "Truth"), lty=c(1, 2, 1), col=c("purple", "purple", "black"))
dev.off()

pdf(file=paste0(figDir, "testing/testWs_SuccessiveSPDE_spde.pdf"), width=8, height=8)
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








# fit model (SPDE) ----
minN=4
out = fitSuccessiveSPDEsimDat(wellTestDat_successiveIPP, seismicTestDat_successiveIPP, 
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

pdf(file=paste0(figDir, "testing/testPreds_SuccessiveSPDE_spde_succesiveIPP.pdf"), width=8, height=8)
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

pdf(file=paste0(figDir, "testing/testAgg_SuccessiveSPDE_spde_succesiveIPP.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (Successive SPDE, pproc)", 
     xlab="Volume Fraction", freq=F)
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
legend("topright", legend=c("Mean", "80% CI", "Truth"), lty=c(1, 2, 1), col=c("purple", "purple", "black"))
dev.off()

pdf(file=paste0(figDir, "testing/testWs_SuccessiveSPDE_spde_succesiveIPP.pdf"), width=8, height=8)
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