# script for testing simulated test dataset

# seismic data, truth ----
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

# naive estimators ----
pSeismic = bilinearInterp(cbind(pEast, pNorth), seismicTestDat, 
                          transform=logit, invTransform = expit)
pTruth = bilinearInterp(cbind(pEast, pNorth), cbind(gEast, gNorth, gTruth), 
                        transform=logit, invTransform = expit)
lmDat = data.frame(logitVfrac = logit(pVolFrac), logitSeismic=logit(pSeismic))
lmMod = lm(logitVfrac~logitSeismic, data=lmDat)
logitGridPreds = predict(lmMod, data.frame(logitSeismic=logit(gSeismic)))
naiveAggEst = mean(logitNormMean(cbind(logitGridPreds, summary(lmMod)$sigma)))

library(survey)
svyDat = cbind(lmDat, seismic=pSeismic, truth=pTruth)
design = svydesign(~1, weights=~I(1/pSeismic), data=svyDat)
svyMod = svyglm(logitVfrac~logitSeismic, design)
logitGridPreds = predict(svyMod, data.frame(logitSeismic=logit(gSeismic)))
wts = 1/pSeismic * (1/sum(1/pSeismic))
sigmaHat = sqrt(1/sum(wts) * sum(wts * residuals(svyMod)^2))
svyAggEst = mean(logitNormMean(cbind(logitGridPreds, sigmaHat)))

trueAgg = mean(gTruth)

trueAgg
svyAggEst
naiveAggEst

# plot results ----

pdf(file=paste0(figDir, "testing/testDatTruth.pdf"), width=5, height=5)
# par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
# squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), 
#        zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
#        asp=1, legend.args=list(axis.args=list(labels=FALSE, tick=FALSE), 
#                               legend.cex=1, smallplot=c(.825,.86,.19,.85), 
#                               legend.lab="urban", legend.line=0.5))
# quilt.plot(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), 
#            zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
#            asp=1, legend.args=list(axis.args=list(labels=FALSE, tick=FALSE), 
#                                    legend.cex=1, smallplot=c(.825,.86,.19,.85), 
#                                    legend.lab="urban", legend.line=0.5))
dev.off()

pdf(file=paste0(figDir, "testing/testDatSeismic.pdf"), width=5, height=5)
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
dev.off()

pdf(file=paste0(figDir, "testing/testDatWell.pdf"), width=5, height=5)
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
splot(pEast, pNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
      xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))
dev.off()

# plot all together
pdf(file=paste0(figDir, "testing/testData.pdf"), width=8, height=8)
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
      asp=1, smallplot=c(.83,.87,.25,.8), cex=.6)
dev.off()

# Repeat for test data sampled preferentially wrt truth ----
gEast = seismicTestDat_truthPref$east
gNorth = seismicTestDat_truthPref$north
gSeismic = seismicTestDat_truthPref$seismicEst
gTruth = truthTestDat_truthPref$truth
pEast = wellTestDat_truthPref$east
pNorth = wellTestDat_truthPref$north
pVolFrac = wellTestDat_truthPref$volFrac

print(cor(gSeismic, gTruth))

# naive estimators ----
pSeismic = bilinearInterp(cbind(pEast, pNorth), seismicTestDat, 
                          transform=logit, invTransform = expit)
pTruth = bilinearInterp(cbind(pEast, pNorth), cbind(gEast, gNorth, gTruth), 
                        transform=logit, invTransform = expit)
lmDat = data.frame(logitVfrac = logit(pVolFrac), logitSeismic=logit(pSeismic))
lmMod = lm(logitVfrac~logitSeismic, data=lmDat)
logitGridPreds = predict(lmMod, data.frame(logitSeismic=logit(gSeismic)))
naiveAggEst = mean(logitNormMean(cbind(logitGridPreds, summary(lmMod)$sigma)))

library(survey)
svyDat = cbind(lmDat, seismic=pSeismic, truth=pTruth)
design = svydesign(~1, weights=~I(1/truth), data=svyDat)
svyMod = svyglm(logitVfrac~logitSeismic, design)
logitGridPreds = predict(svyMod, data.frame(logitSeismic=logit(gSeismic)))
wts = 1/pTruth * (1/sum(1/pTruth))
sigmaHat = sqrt(1/sum(wts) * sum(wts * residuals(svyMod)^2))
svyAggEst = mean(logitNormMean(cbind(logitGridPreds, sigmaHat)))

trueAgg = mean(gTruth)

trueAgg
svyAggEst
naiveAggEst

# plot results ----

pdf(file=paste0(figDir, "testing/testDatTruth_prefTruth.pdf"), width=5, height=5)
# par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
dev.off()

pdf(file=paste0(figDir, "testing/testDatSeismic_prefTruth.pdf"), width=5, height=5)
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
dev.off()

pdf(file=paste0(figDir, "testing/testDatWell_prefTruth.pdf"), width=5, height=5)
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
splot(pEast, pNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
      xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))
dev.off()

# plot all together
pdf(file=paste0(figDir, "testing/testData_prefTruth.pdf"), width=8, height=8)
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
      asp=1, smallplot=c(.83,.87,.25,.8), cex=.6)
dev.off()


# Repeat for successive test data ----
gEast = seismicTestDat_successiveIPP$east
gNorth = seismicTestDat_successiveIPP$north
gSeismic = seismicTestDat_successiveIPP$seismicEst
gTruth = truthTestDat_successiveIPP$truth
pEast = wellTestDat_successiveIPP$east
pNorth = wellTestDat_successiveIPP$north
pVolFrac = wellTestDat_successiveIPP$volFrac

print(cor(gSeismic, gTruth))

# minimum distance between well points
distMat = rdist(cbind(pEast, pNorth))
minDists = sapply(1:nrow(distMat), function(i) {min(distMat[i, -i])})
min(minDists)

# naive estimators ----
pSeismic = bilinearInterp(cbind(pEast, pNorth), seismicTestDat, 
                          transform=logit, invTransform = expit)
pTruth = bilinearInterp(cbind(pEast, pNorth), cbind(gEast, gNorth, gTruth), 
                        transform=logit, invTransform = expit)
lmDat = data.frame(logitVfrac = logit(pVolFrac), logitSeismic=logit(pSeismic))
lmMod = lm(logitVfrac~logitSeismic, data=lmDat)
logitGridPreds = predict(lmMod, data.frame(logitSeismic=logit(gSeismic)))
naiveAggEst = mean(logitNormMean(cbind(logitGridPreds, summary(lmMod)$sigma)))

library(survey)
svyDat = cbind(lmDat, seismic=pSeismic, truth=pTruth)
design = svydesign(~1, weights=~I(1/truth), data=svyDat)
svyMod = svyglm(logitVfrac~logitSeismic, design)
logitGridPreds = predict(svyMod, data.frame(logitSeismic=logit(gSeismic)))
wts = 1/pTruth * (1/sum(1/pTruth))
sigmaHat = sqrt(1/sum(wts) * sum(wts * residuals(svyMod)^2))
svyAggEst = mean(logitNormMean(cbind(logitGridPreds, sigmaHat)))

trueAgg = mean(gTruth)

trueAgg
svyAggEst
naiveAggEst

# plot results ----

pdf(file=paste0(figDir, "testing/testDatTruth_successiveIPP.pdf"), width=5, height=5)
# par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
dev.off()

pdf(file=paste0(figDir, "testing/testDatSeismic_successiveIPP.pdf"), width=5, height=5)
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
dev.off()

pdf(file=paste0(figDir, "testing/testDatWell_successiveIPP.pdf"), width=5, height=5)
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
splot(pEast, pNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
      xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))
dev.off()

# plot all together
pdf(file=paste0(figDir, "testing/testData_successiveIPP.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
text(pEast, pNorth, labels=as.character(1:length(pEast)), cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
text(pEast, pNorth, labels=as.character(1:length(pEast)), cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8), cex=.6)

plot(1:length(pSeismic), pSeismic, pch=19, cex=.5, main="Sampled values", 
     xlab="Iteration number", ylab="Value", ylim=c(0, 1), col="blue")
points(1:length(pSeismic), pTruth, pch=19, cex=.5, col="green")
legend("bottomright", c("Seismic estimate", "Truth"), pch=19, col=c("blue", "green"))

dev.off()

# plot ratio of truth to seismic estimate at well points through time


