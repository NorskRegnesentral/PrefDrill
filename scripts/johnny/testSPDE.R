# script for testing SPDE model

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

divCols = function(n) {redBlueDivCols(n, rev=TRUE)}

pdf(file=paste0(figDir, "testing/testDatTruth.pdf"), width=5, height=5)
# par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=divCols, 
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
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=divCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
dev.off()

pdf(file=paste0(figDir, "testing/testDatWell.pdf"), width=5, height=5)
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
splot(pEast, pNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
      xlim=simStudyXlims, ylim=simStudyYlims, colScale=divCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
       asp=1, smallplot=c(.83,.87,.25,.8))
dev.off()

# plot all together
pdf(file=paste0(figDir, "testing/testData.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=divCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=divCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=divCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))
dev.off()

# fit model ----
out = fitSPDEsimDat(wellTestDat, seismicTestDat)

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

pdf(file=paste0(figDir, "testing/testSPDE.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=divCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=divCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=divCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=divCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (SPDE)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testSPDE_agg.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (SPDE)", 
     xlab="Volume Fraction", freq=F)
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
legend("topright", legend=c("Mean", "80% CI", "Truth"), lty=c(1, 2, 1), col=c("purple", "purple", "black"))
dev.off()


