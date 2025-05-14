# script for testing SPDE model

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

# fit model (without kde) ----

out = fitSPDEsimDat(wellTestDat, seismicTestDat)

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

pdf(file=paste0(figDir, "testing/testPreds_SPDE_seismicPref.pdf"), width=8, height=8)
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
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data, seismicPref", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (SPDE)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SPDE_seismicPref.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (SPDE)", 
     xlab="Volume Fraction", freq=F, xlim=c(0,1))
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
abline(v=svyAggEst, col="blue")
legend("topright", legend=c("Mean", "80% CI", "Truth", "HT"), lty=c(1, 2, 1, 1), col=c("purple", "purple", "black", "blue"))
dev.off()

# fit model (with kde covariate) ----
out = fitSPDEsimDat(wellTestDat, seismicTestDat, addKDE=TRUE)

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

pdf(file=paste0(figDir, "testing/testPreds_SPDE_kde_seismicPref.pdf"), width=8, height=8)
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
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data, seismicPref", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (SPDE, kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SPDE_kde_seismicPref.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (SPDE, kde, seismicPref)", 
     xlab="Volume Fraction", freq=F, xlim=c(0,1))
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
abline(v=svyAggEst, col="blue")
legend("topright", legend=c("Mean", "80% CI", "Truth", "HT"), lty=c(1, 2, 1, 1), col=c("purple", "purple", "black", "blue"))
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

# fit model (without kde) ----
out = fitSPDEsimDat(wellTestDat_truthPref, seismicTestDat_truthPref)

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

pdf(file=paste0(figDir, "testing/testPreds_SPDE_truthPref.pdf"), width=8, height=8)
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
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (SPDE)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SPDE_truthPref.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (SPDE)", 
     xlab="Volume Fraction", freq=F, xlim=c(0,1))
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
abline(v=svyAggEst, col="blue")
legend("topright", legend=c("Mean", "80% CI", "Truth", "HT"), lty=c(1, 2, 1, 1), col=c("purple", "purple", "black", "blue"))
dev.off()

# fit model (with kde covariate) ----
out = fitSPDEsimDat(wellTestDat_truthPref, seismicTestDat_truthPref, addKDE=TRUE)

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

pdf(file=paste0(figDir, "testing/testPreds_SPDE_kde_truthPref.pdf"), width=8, height=8)
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
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data, truthPref", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (SPDE, kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SPDE_kde_truthPref.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (SPDE, kde, truthPref)", 
     xlab="Volume Fraction", freq=F, xlim=c(0,1))
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
abline(v=svyAggEst, col="blue")
legend("topright", legend=c("Mean", "80% CI", "Truth", "HT"), lty=c(1, 2, 1, 1), col=c("purple", "purple", "black", "blue"))
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

# fit model (without kde) ----
out = fitSPDEsimDat(wellTestDat_successiveIPP, seismicTestDat_successiveIPP)

summary(out$mod)

out$fixedEffectSummary
# Current:
#         mean        sd  0.5quant   0.1quant  0.9quant
# X1 0.1419464 0.2573685 0.1397473 -0.1748361 0.4621164 # Int
# X2 0.4600541 0.1550950 0.4595754  0.2623176 0.6583665 # Seismic

# Before changing seismic coef prior
#         mean        sd  0.5quant   0.1quant  0.9quant
# X1 0.1354510 0.2482634 0.1338886 -0.1707699 0.4440865 # Int
# X2 0.4000945 0.1634691 0.4003490  0.1913543 0.6084619 # Seismic

out$parameterSummaryTable
#                       Est           SD       Qlower          Q50       Qupper
# totalVar     4.475239e-01 9.556565e-02 3.141938e-01 4.328827e-01 5.981422e-01
# spatialVar   4.375301e-01 9.556078e-02 3.040138e-01 4.228923e-01 5.883281e-01
# errorVar     9.993862e-03 1.732387e-04 9.814147e-03 9.990437e-03 1.017996e-02
# totalSD      6.651788e-01 7.117732e-02 5.605299e-01 6.579382e-01 7.733965e-01
# spatialSD    6.575344e-01 7.199866e-02 5.513745e-01 6.503017e-01 7.670255e-01
# errorSD      9.996556e-02 8.661572e-04 9.906638e-02 9.995217e-02 1.008958e-01
# spatialRange 3.561494e+03 9.032152e+02 2.504507e+03 3.446784e+03 4.569013e+03

names(out)

out$timings
# totalTime modelDefineTime modelFitTime posteriorSamplingTime otherTime modelDefinePct modelFitTimePct
# elapsed     42.59            1.23        11.11                 26.67      3.58     0.02888002       0.2608594
# posteriorSamplingTimePct otherTimePct
# elapsed                0.6262033   0.08405729

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

pdf(file=paste0(figDir, "testing/testPreds_SPDE_successiveIPP.pdf"), width=8, height=8)
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
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (SPDE)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SPDE_successiveIPP.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (SPDE)", 
     xlab="Volume Fraction", freq=F, xlim=c(0,1))
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
abline(v=svyAggEst, col="blue")
legend("topright", legend=c("Mean", "80% CI", "Truth", "HT"), lty=c(1, 2, 1, 1), col=c("purple", "purple", "black", "blue"))
dev.off()

#  Calculate MSEs ----
print(paste0("MSE (seismic and SPDE): ", mean(gTruth - gSeismic)^2, " and ", mean(gTruth - preds)^2))
# [1] "MSE (seismic and SPDE): 0.000575049976194529 and 0.000488413144050024" # current
# [1] "MSE (seismic and SPDE): 0.000575049976194529 and 0.000668483827213597" # before changing seismic prior
print(paste0("logit MSE (seismic and SPDE): ", mean(logit(gTruth) - logit(gSeismic))^2, " and ", mean(logit(gTruth) - logit(preds))^2))
# [1] "logit MSE (seismic and SPDE): 0.00918437079448658 and 0.0118368850592009" # current
# [1] "logit MSE (seismic and SPDE): 0.00918437079448658 and 0.0156641805988305" # before changing seismic prior

print(paste0("MSE (seismic and SPDE): ", mean(pTruth - pSeismic)^2, " and ", mean(pTruth - obsPreds)^2))
# [1] "MSE (seismic and SPDE): 9.36241385035057e-07 and 6.24006759552452e-06" # current
# [1] "MSE (seismic and SPDE): 9.36241385035057e-07 and 8.60959238767509e-06" # before changing seismic prior
print(paste0("logit MSE (seismic and SPDE): ", mean(logit(pTruth) - logit(pSeismic))^2, " and ", mean(logit(pTruth) - logit(obsPreds))^2))
# [1] "logit MSE (seismic and SPDE): 0.000706596851712083 and 0.000160046362674348" # current
# [1] "logit MSE (seismic and SPDE): 0.000706596851712083 and 0.000213128148633414" # before changing seismic prior

tempMod = lm(pVolFrac ~ pSeismic)
summary(tempMod)

# fit model (with kde covariate) ----
out = fitSPDEsimDat(wellTestDat_successiveIPP, seismicTestDat_successiveIPP, addKDE=TRUE)

summary(out$mod)

out$fixedEffectSummary
# current
#          mean        sd   0.5quant    0.1quant   0.9quant
# X1  0.2597330 0.2754897  0.2616151 -0.08341869  0.6002660
# X2  0.6618511 0.1807942  0.6618591  0.43107868  0.8925562
# X3 -0.3695933 0.1949673 -0.3695592 -0.61822008 -0.1209852

out$parameterSummaryTable
#                       Est           SD       Qlower          Q50       Qupper
# totalVar     4.240591e-01 9.187172e-02 2.969014e-01 4.088807e-01 5.724527e-01
# spatialVar   4.140599e-01 9.186522e-02 2.867205e-01 3.988886e-01 5.626378e-01
# errorVar     9.999128e-03 1.744119e-04 9.814889e-03 9.992054e-03 1.018095e-02
# totalSD      6.474111e-01 7.016268e-02 5.448866e-01 6.394378e-01 7.566061e-01
# spatialSD    6.395482e-01 7.101476e-02 5.354629e-01 6.315763e-01 7.500919e-01
# errorSD      9.999184e-02 8.719808e-04 9.907012e-02 9.996026e-02 1.009007e-01
# spatialRange 3.921780e+03 9.807687e+02 2.786811e+03 3.783173e+03 5.045601e+03

names(out)

out$timings
# totalTime modelDefineTime modelFitTime posteriorSamplingTime otherTime modelDefinePct modelFitTimePct
# elapsed     47.34               1        12.22                 31.77      2.35     0.02112379       0.2581327
# posteriorSamplingTimePct otherTimePct
# elapsed                0.6711027    0.0496409

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

pdf(file=paste0(figDir, "testing/testPreds_SPDE_kde_successiveIPP.pdf"), width=8, height=8)
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
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data, Successive IPP", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (SPDE, kde)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SPDE_kde_successiveIPP.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (SPDE, kde, Successive IPP)", 
     xlab="Volume Fraction", freq=F, xlim=c(0,1))
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
abline(v=svyAggEst, col="blue")
legend("topright", legend=c("Mean", "80% CI", "Truth", "HT"), lty=c(1, 2, 1, 1), col=c("purple", "purple", "black", "blue"))
dev.off()

pdf(file=paste0(figDir, "testing/testPair_SPDE_kde_successiveIPP.pdf"), width=5, height=5)
plot(gTruth, gSeismic, main="Predictions vs truth and seismic (SPDE)", 
     xlab="Truth", ylab="Estimate", xlim=c(0,1), ylim=c(0, 1), pch=19, col="black", cex=.1)
points(gTruth, preds, pch=19, col="blue", cex=.1)
abline(0, 1)
legend("topleft", legend=c("Seismic", "SPDE"), pch=19, col=c("black", "blue"))
dev.off()

lims = range(logit(c(gTruth, gSeismic, preds)))
pdf(file=paste0(figDir, "testing/testPairLogit_SPDE_kde_successiveIPP.pdf"), width=5, height=5)
plot(logit(gTruth), logit(gSeismic), main="Logit predictions vs truth and seismic (SPDE)", 
     xlab="Truth", ylab="Estimate", xlim=lims, ylim=lims, pch=19, col="black", cex=.1)
points(logit(gTruth), logit(preds), pch=19, col="blue", cex=.1)
abline(0, 1)
legend("topleft", legend=c("Seismic", "SPDE"), pch=19, col=c("black", "blue"))
dev.off()

#  Calculate MSEs ----
print(paste0("MSE (seismic and SPDE): ", mean(gTruth - gSeismic)^2, " and ", mean(gTruth - preds)^2))
# [1] "MSE (seismic and SPDE): 0.000575049976194529 and 0.00246729193698229" # current
# [1] "MSE (seismic and SPDE): 0.000575049976194529 and 0.00243152077758097" # before adding seismic prior
print(paste0("logit MSE (seismic and SPDE): ", mean(logit(gTruth) - logit(gSeismic))^2, " and ", mean(logit(gTruth) - logit(preds))^2))
# [1] "logit MSE (seismic and SPDE): 0.00918437079448658 and 0.0509517964478948" # current
# [1] "logit MSE (seismic and SPDE): 0.00918437079448658 and 0.0501667280153215" # before adding seismic prior

print(paste0("MSE (seismic and SPDE): ", mean(pTruth - pSeismic)^2, " and ", mean(pTruth - obsPreds)^2))
# [1] "MSE (seismic and SPDE): 9.36241385035057e-07 and 7.35800194044614e-06" # current
# [1] "MSE (seismic and SPDE): 9.36241385035057e-07 and 6.24781947966474e-06" # before adding seismic prior
print(paste0("logit MSE (seismic and SPDE): ", mean(logit(pTruth) - logit(pSeismic))^2, " and ", mean(logit(pTruth) - logit(obsPreds))^2))
# [1] "logit MSE (seismic and SPDE): 0.000706596851712083 and 0.000192303145691416" # current
# [1] "logit MSE (seismic and SPDE): 0.000706596851712083 and 0.000164310870152474" # before adding seismic prior

# fit model (without kde, strict spatial prior) ----
mesh=getSPDEmeshSimStudy()
prior=getSPDEprior(mesh, U=1, alpha=0.01)
out = fitSPDEsimDat(wellTestDat_successiveIPP, seismicTestDat_successiveIPP, mesh=mesh, prior=prior)

summary(out$mod)

out$fixedEffectSummary
#         mean        sd  0.5quant   0.1quant  0.9quant
# X1 0.1339131 0.2291812 0.1325140 -0.1503416 0.4202935 # Int
# X2 0.3952025 0.1599935 0.3954386  0.1908322 0.5992302 # Seismic

out$parameterSummaryTable
#                       Est           SD       Qlower          Q50       Qupper
# totalVar     4.044620e-01 7.934443e-02 2.931786e-01 3.937815e-01 5.268410e-01
# spatialVar   3.942288e-01 7.942296e-02 2.838826e-01 3.833713e-01 5.184387e-01
# errorVar     1.023315e-02 1.794580e-03 8.376822e-03 9.971021e-03 1.218642e-02
# totalSD      6.329278e-01 6.219521e-02 5.414597e-01 6.275201e-01 7.258381e-01
# spatialSD    6.247048e-01 6.306107e-02 5.328063e-01 6.191699e-01 7.200269e-01
# errorSD      1.007758e-01 8.801881e-03 9.152498e-02 9.985500e-02 1.103921e-01
# spatialRange 3.240931e+03 8.168992e+02 2.262528e+03 3.157367e+03 4.130542e+03

names(out)

out$timings
# totalTime modelDefineTime modelFitTime posteriorSamplingTime otherTime modelDefinePct modelFitTimePct
# elapsed     42.59            1.23        11.11                 26.67      3.58     0.02888002       0.2608594
# posteriorSamplingTimePct otherTimePct
# elapsed                0.6262033   0.08405729

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

pdf(file=paste0(figDir, "testing/testPreds_SPDE_alpha01_successiveIPP.pdf"), width=8, height=8)
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
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (SPDE)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_SPDE_alpha01_successiveIPP.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (SPDE)", 
     xlab="Volume Fraction", freq=F, xlim=c(0,1))
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
abline(v=svyAggEst, col="blue")
legend("topright", legend=c("Mean", "80% CI", "Truth", "HT"), lty=c(1, 2, 1, 1), col=c("purple", "purple", "black", "blue"))
dev.off()

#  Calculate MSEs ----
print(paste0("MSE (seismic and SPDE): ", mean(gTruth - gSeismic)^2, " and ", mean(gTruth - preds)^2))
# [1] "MSE (seismic and SPDE): 0.000575049976194529 and 0.000619012685497106"
print(paste0("logit MSE (seismic and SPDE): ", mean(logit(gTruth) - logit(gSeismic))^2, " and ", mean(logit(gTruth) - logit(preds))^2))
# [1] "logit MSE (seismic and SPDE): 0.00918437079448658 and 0.014654008303109"

print(paste0("MSE (seismic and SPDE): ", mean(pTruth - pSeismic)^2, " and ", mean(pTruth - obsPreds)^2))
# [1] "MSE (seismic and SPDE): 9.36241385035057e-07 and 7.38644495813154e-06"
print(paste0("logit MSE (seismic and SPDE): ", mean(logit(pTruth) - logit(pSeismic))^2, " and ", mean(logit(pTruth) - logit(obsPreds))^2))
# [1] "logit MSE (seismic and SPDE): 0.000706596851712083 and 0.000188186827925881"
