# script for testing simulation study data

# setup ----
thisFigDir =  "figures/testSimStudy/"

out = load("C:/Users/jpaige/git/PrefDrill/savedOutput/simStudy/simParListBatch.RData")
sampleParCombs

simStudyXlims <<- c(-12.5, 15012.5)
simStudyYlims <<- c(-8.3335, 5008.4336)

# batch case -----

adaptScen = "batch"

parI = 9
repI = 3
sampleParCombs[parI,]

# seismic data
out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedPred_", repI, ".txt"), force01=TRUE)
seismicDat = out$surfFrame

# truth
out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", repI, ".txt"), force01=TRUE)
truth = out$surfFrame

# subsample
goodCoords = subsampleSimStudyGrid(seismicDat)
seismicDat = seismicDat[goodCoords,]
truth = truth[goodCoords,]

# well data
wellDatFile = paste0("savedOutput/simStudy/wellDat/wellDat_", adaptScen, "_par", parI, "_rep", repI, ".RData")
out = load(wellDatFile)

maxN = 20
wellDat = wellDat[1:maxN,]

# plot some figures for testing purposes

gEast = seismicDat$east
gNorth = seismicDat$north
gSeismic = seismicDat$seismicEst
gTruth = truth[,3]
pEast = wellDat$east
pNorth = wellDat$north
pVolFrac = wellDat$volFrac

eastGrid = sort(unique(gEast))
northGrid = sort(unique(gNorth))

ticks = seq(0, 1, by=.2)
tickLabs = as.character(ticks)

seqCols = function(n) {purpleYellowSeqCols(n)}

print(cor(gSeismic, gTruth))

# minimum distance between well points
distMat = rdist(cbind(pEast, pNorth))
minDists = sapply(1:nrow(distMat), function(i) {min(distMat[i, -i])})
min(minDists)

## naive estimators ----
pSeismic = bilinearInterp(cbind(pEast, pNorth), cbind(gEast, gNorth, gSeismic), 
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

## plot results ----

pdf(file=paste0(thisFigDir, "simDatTruth_par", parI, "_rep", repI, ".pdf"), width=5, height=5)
# par(mar=c(3, 3, 2, 5), mgp=c(1.7, .5, 0))
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
points(pEast, pNorth, cex=.3, lwd=.25)
dev.off()

pdf(file=paste0(thisFigDir, "simDatSeismic_par", parI, "_rep", repI, ".pdf"), width=5, height=5)
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
points(pEast, pNorth, cex=.3, lwd=.25)
dev.off()

pdf(file=paste0(thisFigDir, "simDatWell_par", parI, "_rep", repI, ".pdf"), width=5, height=5)
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
splot(pEast, pNorth, pVolFrac, 
      xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8), cex=.3, lwd=.25)
dev.off()

pdf(file=paste0(thisFigDir, "simDatWellTruth_par", parI, "_rep", repI, ".pdf"), width=5, height=5)
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
splot(pEast, pNorth, pTruth, 
      xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8), cex=.3, lwd=.25)
dev.off()

pdf(file=paste0(thisFigDir, "simDatWellSeismic_par", parI, "_rep", repI, ".pdf"), width=5, height=5)
par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
splot(pEast, pNorth, pSeismic, 
      xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8), cex=.3, lwd=.25)
dev.off()

# plot all together
pdf(file=paste0(thisFigDir, "simDat_par", parI, "_rep", repI, ".pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
text(pEast, pNorth, labels=as.character(1:length(pEast)), cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
text(pEast, pNorth, labels=as.character(1:length(pEast)), cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, legend.args=list(smallplot=c(.83,.87,.25,.8)), cex=.3, lwd=.25)

plot(1:length(pSeismic), pSeismic, pch=19, cex=.5, main="Sampled values", 
     xlab="Iteration number", ylab="Value", ylim=c(0, 1), col="blue")
points(1:length(pSeismic), pTruth, pch=19, cex=.5, col="green")
legend("bottomright", c("Seismic estimate", "Truth"), pch=19, col=c("blue", "green"))

dev.off()

























pdf(file=paste0("figures/testSimStudy/testPreds_simStudy_", adaptScen, "_par", 
                sampleParI, "_rep", repI, ".pdf"), width=8, height=5)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
points(pEast, pNorth, cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0("figures/testSimStudy/testQuants_simStudy_", adaptScen, "_par", 
                sampleParI, "_rep", repI, ".pdf"), width=8, height=5)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
points(pEast, pNorth, cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, legend.args=list(smallplot=c(.83,.87,.25,.8)), ticks=ticks, tickLabels=tickLabs)

squilt(gEast, gNorth, predQuants, grid=list(x=eastGrid, y=northGrid), 
       colScale=quantCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Truth Quantiles", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0("figures/testSimStudy/testLogitSDs_simStudy_", adaptScen, "_par", 
                sampleParI, "_rep", repI, ".pdf"), width=8, height=5)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, logit(gTruth), grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
points(pEast, pNorth, cex=.5)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, legend.args=list(smallplot=c(.83,.87,.25,.8)), ticks=ticks, tickLabels=tickLabs)

squilt(gEast, gNorth, predQuants, grid=list(x=eastGrid, y=northGrid), 
       colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Truth Quantiles", 
       asp=1, smallplot=c(.83,.87,.25,.8), ticks=ticks, tickLabels=tickLabs)
points(pEast, pNorth, cex=.5)
dev.off()