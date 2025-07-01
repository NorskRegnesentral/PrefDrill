# script for testing Watson et al. model

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

pdf(file=paste0(figDir, "testing/testPreds_Watson_seismicPref.pdf"), width=8, height=8)
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

pdf(file=paste0(figDir, "testing/testAgg_Watson_seismicPref.pdf"), width=5, height=5)
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

pdf(file=paste0(figDir, "testing/testPreds_Watson_kde_seismicPref.pdf"), width=8, height=8)
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

pdf(file=paste0(figDir, "testing/testAgg_Watson_kde_seismicPref.pdf"), width=5, height=5)
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

pdf(file=paste0(figDir, "testing/testPreds_Watson_truthPref.pdf"), width=8, height=8)
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

pdf(file=paste0(figDir, "testing/testAgg_Watson_truthPref.pdf"), width=5, height=5)
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

pdf(file=paste0(figDir, "testing/testPreds_Watson_kde_truthPref.pdf"), width=8, height=8)
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

pdf(file=paste0(figDir, "testing/testAgg_Watson_kde_truthPref.pdf"), width=5, height=5)
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
out = fitWatsonSimDat(wellTestDat_successiveIPP, seismicTestDat_successiveIPP)

summary(out$mod)
out$fixedEffectSummary

# current
#              mean        sd   0.5quant   0.1quant   0.9quant
# X.y1    0.1459610 0.2568172  0.1450624 -0.1737417  0.4669562
# X.y2    0.4604633 0.1543527  0.4601099  0.2634151  0.6579391
# X.pp1  -8.3717763 1.0029317 -8.3717705 -9.6576995 -7.0858605
# X.pp2  -8.3661212 1.0029793 -8.3661155 -9.6521053 -7.0801445
# X.pp3  -8.3519629 1.0028375 -8.3519577 -9.6377648 -7.0661676
# X.pp4  -8.3491821 1.0028942 -8.3491769 -9.6350568 -7.0633140
# X.pp5  -8.3451006 1.0029776 -8.3450954 -9.6310823 -7.0591256
# X.pp6  -8.3289746 1.0027121 -8.3289694 -9.6146159 -7.0433402
# X.pp7  -8.3220262 1.0027526 -8.3220209 -9.6077194 -7.0363398
# X.pp8  -8.3064680 1.0025854 -8.3064631 -9.5919466 -7.0209958
# X.pp9  -8.2880496 1.0023142 -8.2880445 -9.5731805 -7.0029251
# X.pp10 -8.2817135 1.0023603 -8.2817085 -9.5669035 -6.9965299
# X.pp11 -8.2738898 1.0023808 -8.2738849 -9.5591061 -6.9886799
# X.pp12 -8.2647045 1.0023952 -8.2646995 -9.5499393 -6.9794761
# X.pp13 -8.2537026 1.0023598 -8.2536976 -9.5388920 -6.9685196
# X.pp14 -8.2432482 1.0023389 -8.2432434 -9.5284107 -6.9580920
# X.pp15 -8.2378450 1.0024108 -8.2378402 -9.5230996 -6.9525964
# X.pp16 -8.2239422 1.0022171 -8.2239374 -9.5089485 -6.9389420
# X.pp17 -8.2151223 1.0021704 -8.2151176 -9.5000688 -6.9301820
# X.pp18 -8.2019628 1.0020211 -8.2019583 -9.4867176 -6.9172137
# X.pp19 -8.1958252 1.0020851 -8.1958208 -9.4806620 -6.9109939
# X.pp20 -8.1800085 1.0019942 -8.1800039 -9.4647287 -6.8952940
# X.pp21  0.5466815 0.2083730  0.5466438  0.2795451  0.8138659

# before seismic prior
#              mean        sd   0.5quant   0.1quant   0.9quant
# X.y1    0.1436738 0.2475690  0.1432831 -0.1651981  0.4530783
# X.y2    0.4002405 0.1634461  0.4004527  0.1913370  0.6088452
# X.pp1  -8.3537734 1.0020151 -8.3537651 -9.6385231 -7.0690345
# X.pp2  -8.3478585 1.0020533 -8.3478501 -9.6326572 -7.0630705
# X.pp3  -8.3344883 1.0019331 -8.3344805 -9.6191326 -7.0498541
# X.pp4  -8.3312719 1.0019846 -8.3312641 -9.6159821 -7.0465717
# X.pp5  -8.3262610 1.0020581 -8.3262535 -9.6110654 -7.0414663
# X.pp6  -8.3120621 1.0018524 -8.3120544 -9.5966029 -7.0275312
# X.pp7  -8.3048878 1.0018806 -8.3048800 -9.5894647 -7.0203208
# X.pp8  -8.2902996 1.0017406 -8.2902924 -9.5746968 -7.0059118
# X.pp9  -8.2740463 1.0015417 -8.2740390 -9.5581885 -6.9899136
# X.pp10 -8.2674967 1.0015745 -8.2674892 -9.5516810 -6.9833220
# X.pp11 -8.2592856 1.0015858 -8.2592784 -9.5434841 -6.9750962
# X.pp12 -8.2502633 1.0015926 -8.2502559 -9.5344708 -6.9660653
# X.pp13 -8.2396788 1.0015586 -8.2396713 -9.5238426 -6.9555245
# X.pp14 -8.2291040 1.0015355 -8.2290970 -9.5132380 -6.9449790
# X.pp15 -8.2232677 1.0015916 -8.2232606 -9.5074738 -6.9390708
# X.pp16 -8.2109078 1.0014528 -8.2109008 -9.4949358 -6.9268889
# X.pp17 -8.2025541 1.0014142 -8.2025471 -9.4865325 -6.9185848
# X.pp18 -8.1904547 1.0013058 -8.1904481 -9.4742940 -6.9066239
# X.pp19 -8.1839631 1.0013508 -8.1839564 -9.4678601 -6.9000748
# X.pp20 -8.1696060 1.0012875 -8.1695987 -9.4534222 -6.8857992
# X.pp21  0.4499548 0.2287569  0.4499411  0.1566673  0.7432595 # seismic estimate

out$parameterSummaryTable
#                        Est           SD        Qlower           Q50       Qupper
# totalVar        0.44404892 8.563846e-02  3.275257e-01  4.331413e-01 5.875535e-01
# spatialVar      0.43404790 8.564066e-02  3.175274e-01  4.231433e-01 5.777115e-01
# errorVar        0.01000102 1.523966e-04  9.842362e-03  9.998005e-03 1.015937e-02
# totalSD         0.66329526 6.397197e-02  5.722986e-01  6.581347e-01 7.665204e-01
# spatialSD       0.65563995 6.471741e-02  5.634957e-01  6.504946e-01 7.600733e-01
# errorSD         0.10000222 7.618549e-04  9.920868e-02  9.999002e-02 1.007937e-01
# spatialRange 3544.98016576 7.891418e+02  2.576741e+03  3.542685e+03 4.647723e+03
# prefPar        -0.01598913 1.789553e-01 -2.571784e-01 -2.147669e-02 2.187120e-01

out$timings
#         totalTime modelDefineTime modelFitTime posteriorSamplingTime otherTime modelDefinePct modelFitTimePct
# elapsed    192.09            3.95        79.67                102.22      6.25     0.02056328       0.4147535
# posteriorSamplingTimePct otherTimePct
# elapsed                0.5321464   0.03253683

# plot predictions ----

fixedPreds = rowMeans(out$fixedPredMat.y)
spatialPreds = rowMeans(out$spatialPredMat.y)
preds = out$pred.yEst
obsPreds = out$obs.yEst
preds.pp = out$pred.ppEst
predAggMat = out$pred.yAggMat
predAgg = out$pred.yAggEst
predAggU = out$pred.yAggUpper
predAggL = out$pred.yAggLower
trueAgg = mean(gTruth)
predAgg
trueAgg

pdf(file=paste0(figDir, "testing/testPreds_Watson_successiveIPP.pdf"), width=8, height=8)
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
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Watson)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testPredsEps_Watson_successiveIPP.pdf"), width=8, height=8)
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

splot(pEast, pNorth, pVolFrac-pSeismic, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Watson)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testPredsPPEps_Watson_successiveIPP.pdf"), width=8, height=8)
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

splot(pEast, pNorth, pVolFrac-pSeismic, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))

cellArea = diff(eastGrid[1:2]) * diff(northGrid[1:2])
squilt(gEast, gNorth, exp(preds.pp)/cellArea, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       xlab="Easting", ylab="Northing", main="Next point intensity estimate (Watson)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testPredParts_Watson_successiveIPP.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, fixedPreds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       xlab="Easting", ylab="Northing", main="Fixed Effect", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, logit(gSeismic), grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       xlab="Easting", ylab="Northing", main="Logit Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, spatialPreds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       xlab="Easting", ylab="Northing", main="Spatial Effect", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, logit(preds), grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       xlab="Easting", ylab="Northing", main="Logit Estimate (Watson)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_Watson_successiveIPP.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (Watson)", 
     xlab="Volume Fraction", freq=F, xlim=c(0,1))
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
abline(v=svyAggEst, col="blue")
legend("topright", legend=c("Mean", "80% CI", "Truth", "HT"), lty=c(1, 2, 1, 1), col=c("purple", "purple", "black", "blue"))
dev.off()

pdf(file=paste0(figDir, "testing/testPair_Watson_successiveIPP.pdf"), width=5, height=5)
plot(gTruth, gSeismic, main="Predictions vs truth and seismic (Watson)", 
     xlab="Truth", ylab="Estimate", xlim=c(0,1), ylim=c(0, 1), pch=19, col="black", cex=.1)
points(gTruth, preds, pch=19, col="blue", cex=.1)
abline(0, 1)
legend("topleft", legend=c("Seismic", "Watson"), pch=19, col=c("black", "blue"))
dev.off()

lims = range(logit(c(gTruth, gSeismic, preds)))
pdf(file=paste0(figDir, "testing/testPairLogit_Watson_successiveIPP.pdf"), width=5, height=5)
plot(logit(gTruth), logit(gSeismic), main="Logit predictions vs truth and seismic (Watson)", 
     xlab="Truth", ylab="Estimate", xlim=lims, ylim=lims, pch=19, col="black", cex=.1)
points(logit(gTruth), logit(preds), pch=19, col="blue", cex=.1)
abline(0, 1)
legend("topleft", legend=c("Seismic", "Watson"), pch=19, col=c("black", "blue"))
dev.off()

print(paste0("MSE (seismic and Watson): ", mean(gTruth - gSeismic)^2, " and ", mean(gTruth - preds)^2))
# [1] "MSE (seismic and Watson): 0.000575049976194529 and 0.000457862891827518" # current
# [1] "MSE (seismic and Watson): 0.000575049976194529 and 0.000743311172475913" # before seismic prior
print(paste0("logit MSE (seismic and Watson): ", mean(logit(gTruth) - logit(gSeismic))^2, " and ", mean(logit(gTruth) - logit(preds))^2))
# [1] "logit MSE (seismic and Watson): 0.00918437079448658 and 0.0111454503420435" # current
# [1] "logit MSE (seismic and Watson): 0.00918437079448658 and 0.0171151574649274" # before seismic prior

print(paste0("MSE (seismic and Watson): ", mean(pTruth - pSeismic)^2, " and ", mean(pTruth - obsPreds)^2))
# [1] "logit MSE (seismic and Watson): 0.000706596851712083 and 0.00017709516220975" # current
# [1] "MSE (seismic and Watson): 9.36241385035057e-07 and 6.15114141922741e-06" # before seismic prior
print(paste0("logit MSE (seismic and Watson): ", mean(logit(pTruth) - logit(pSeismic))^2, " and ", mean(logit(pTruth) - logit(obsPreds))^2))
# [1] "logit MSE (seismic and Watson): 0.000706596851712083 and 0.00017709516220975" # current # before seismic prior
# [1] "logit MSE (seismic and Watson): 0.000706596851712083 and 0.000159484950630608"

temp= logit(gTruth)
loess_fit <- loess(I(logit(preds)) ~ temp)

# Predict fitted values
x_vals <- seq(0, 1, length.out = 100)
x_vals <- seq(-3, 3, length.out = 100)
y_pred <- predict(loess_fit, newdata = data.frame(temp = x_vals))

lines(x_vals, y_pred, col="green")

library(ggplot2)
data = data.frame(gTruth=gTruth, preds=preds)
ggplot(data, aes(x = gTruth, y = preds)) +
  geom_point() +                              # Scatter plot
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, color = "blue") +
  theme_minimal()


var(logit(pTruth) - logit(pVolFrac))
# 0.01149867
var(logit(pTruth) - logit(pSeismic))
# 0.6878443
var(logit(gTruth) - logit(gSeismic))
# 0.5432175

# fit model with no repulsion ----
out = fitWatsonSimDat(wellTestDat_successiveIPP, seismicTestDat_successiveIPP, repelDist = 10)

summary(out$mod)
out$fixedEffectSummary
out$parameterSummaryTable

names(out)

out$timings
# totalTime modelDefineTime modelFitTime posteriorSamplingTime otherTime modelDefinePct modelFitTimePct
# elapsed    188.99            6.11        67.63                106.35       8.9     0.03232975       0.3578496
# posteriorSamplingTimePct otherTimePct
# elapsed                0.5627282   0.04709244

# plot predictions ----

fixedPreds = rowMeans(out$fixedPredMat.y)
spatialPreds = rowMeans(out$spatialPredMat.y)
preds = out$pred.yEst
obsPreds = out$obs.yEst
preds.pp = out$pred.ppEst
predAggMat = out$pred.yAggMat
predAgg = out$pred.yAggEst
predAggU = out$pred.yAggUpper
predAggL = out$pred.yAggLower
trueAgg = mean(gTruth)
predAgg
trueAgg

pdf(file=paste0(figDir, "testing/testPreds_Watson_rep10_successiveIPP.pdf"), width=8, height=8)
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
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Watson)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testPredsEps_Watson_rep10_successiveIPP.pdf"), width=8, height=8)
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

splot(pEast, pNorth, pVolFrac-pSeismic, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Watson)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testPredsPPEps_Watson_rep10_successiveIPP.pdf"), width=8, height=8)
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

splot(pEast, pNorth, pVolFrac-pSeismic, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      xlab="Easting", ylab="Northing", main="Well Data", 
      asp=1, smallplot=c(.83,.87,.25,.8))

cellArea = diff(eastGrid[1:2]) * diff(northGrid[1:2])
squilt(gEast, gNorth, exp(preds.pp)/cellArea, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       xlab="Easting", ylab="Northing", main="Next point intensity estimate (Watson)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testPredParts_Watson_rep10_successiveIPP.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, fixedPreds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       xlab="Easting", ylab="Northing", main="Fixed Effect", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, logit(gSeismic), grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       xlab="Easting", ylab="Northing", main="Logit Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, spatialPreds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       xlab="Easting", ylab="Northing", main="Spatial Effect", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)

squilt(gEast, gNorth, logit(preds), grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       xlab="Easting", ylab="Northing", main="Logit Estimate (Watson)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_Watson_rep10_successiveIPP.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (Watson)", 
     xlab="Volume Fraction", freq=F, xlim=c(0,1))
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
abline(v=svyAggEst, col="blue")
legend("topright", legend=c("Mean", "80% CI", "Truth", "HT"), lty=c(1, 2, 1, 1), col=c("purple", "purple", "black", "blue"))
dev.off()

pdf(file=paste0(figDir, "testing/testPair_Watson_rep10_successiveIPP.pdf"), width=5, height=5)
plot(gTruth, gSeismic, main="Predictions vs truth and seismic (Watson)", 
     xlab="Truth", ylab="Estimate", xlim=c(0,1), ylim=c(0, 1), pch=19, col="black", cex=.1)
points(gTruth, preds, pch=19, col="blue", cex=.1)
abline(0, 1)
legend("topleft", legend=c("Seismic", "Watson"), pch=19, col=c("black", "blue"))
dev.off()

lims = range(logit(c(gTruth, gSeismic, preds)))
pdf(file=paste0(figDir, "testing/testPairLogit_Watson_rep10_successiveIPP.pdf"), width=5, height=5)
plot(logit(gTruth), logit(gSeismic), main="Logit predictions vs truth and seismic (Watson)", 
     xlab="Truth", ylab="Estimate", xlim=lims, ylim=lims, pch=19, col="black", cex=.1)
points(logit(gTruth), logit(preds), pch=19, col="blue", cex=.1)
abline(0, 1)
legend("topleft", legend=c("Seismic", "Watson"), pch=19, col=c("black", "blue"))
dev.off()

print(paste0("MSE (seismic and Watson): ", mean(gTruth - gSeismic)^2, " and ", mean(gTruth - preds)^2))
# [1] "MSE (seismic and Watson): 0.000575049976194529 and 0.000655688160261606"
print(paste0("logit MSE (seismic and Watson): ", mean(logit(gTruth) - logit(gSeismic))^2, " and ", mean(logit(gTruth) - logit(preds))^2))
# [1] "logit MSE (seismic and Watson): 0.00918437079448658 and 0.0153906195214439"

print(paste0("MSE (seismic and Watson): ", mean(pTruth - pSeismic)^2, " and ", mean(pTruth - obsPreds)^2))
# [1] "MSE (seismic and Watson): 9.36241385035057e-07 and 7.44658347577478e-06"
print(paste0("logit MSE (seismic and Watson): ", mean(logit(pTruth) - logit(pSeismic))^2, " and ", mean(logit(pTruth) - logit(obsPreds))^2))
# [1] "logit MSE (seismic and Watson): 0.000706596851712083 and 0.000182783198580526"

# _ ----
# Test pseudosites ----
# _ ----
gEast = seismicTestDat_successiveIPP$east
gNorth = seismicTestDat_successiveIPP$north

pdf(file=paste0(figDir, "testing/testPseudo_Watson.pdf"), width=5, height=5)
plot(gEast, gNorth, main="Pseudosites", 
     xlab="Easting", ylab="Northing", xlim=c(0, 10000), ylim=c(0, 10000), pch=19, col="black", cex=.1)
dev.off()
