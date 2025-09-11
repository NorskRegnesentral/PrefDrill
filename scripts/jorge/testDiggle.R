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

# fit model ----
out = fitDigglesimDat(wellTestDat, seismicTestDat)

# plot predictions ----
preds = out$predEst
predAggMat = out$predAggMat
predAgg = out$predAggEst
predAggU = out$predAggUpper
predAggL = out$predAggLower


pdf(file=paste0(figDir, "testing/testPreds_Diggle_seismicPref.pdf"), width=8, height=8)
par(mfrow=c(2,2), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
squilt(gEast, gNorth, gTruth, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
       xlim=simStudyXlims, ylim=simStudyYlims, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="True Sand Volume Frac", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.8)

squilt(gEast, gNorth, gSeismic, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Seismic Estimate", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.8)

splot(pEast, pNorth, pVolFrac, grid=list(x=eastGrid, y=northGrid), colScale=seqCols, 
      xlim=simStudyXlims, ylim=simStudyYlims, resetGraphics=FALSE, 
      zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Well Data, seismicPref", 
      asp=1, smallplot=c(.83,.87,.25,.8))

squilt(gEast, gNorth, preds, grid=list(x=eastGrid, y=northGrid), 
       xlim=simStudyXlims, ylim=simStudyYlims, colScale=seqCols, 
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Diggle)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.8)

dev.off()

pdf(file=paste0(figDir, "testing/testAgg_Diggle_seismicPref.pdf"), width=5, height=5)
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (Diggle)", 
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

# fit model  ----
# goodCoords = subsampleSimStudyGrid(seismicTestDat_truthPref)
# thisSeismic = seismicTestDat_truthPref[goodCoords,]

out = fitDigglesimDat(wellTestDat_truthPref, seismicTestDat_truthPref, 
                      mesh=getSPDEmeshTestDat(), prefMean=2.5)

summary(out$mod)
# prefMean = 0:
# Fixed effects:
#                 mean    sd 0.025quant 0.5quant 0.975quant    mode kld
# Intercept_pp -16.903 0.633    -18.233  -16.874    -15.735 -16.777   0
# X_pp           1.923 0.388      1.185    1.915      2.707   1.915   0
# Intercept_y   -0.093 0.294     -0.699   -0.083      0.462  -0.083   0
# X_y            0.832 0.165      0.512    0.831      1.159   0.831   0
# 
# Random effects:
#   Name	  Model
# field_y SPDE2 model
# field_pp Copy
# 
# Model hyperparameters:
#                                                mean       sd 0.025quant 0.5quant 0.975quant     mode
# Precision for the Gaussian observations[2]  100.005    3.160     93.857   99.979    106.299   99.979
# Range for field_y                          3049.729 1050.036   1418.359 2906.014   5497.542 2642.207
# Stdev for field_y                             0.671    0.148      0.418    0.658      0.999    0.636
# Beta for field_pp                             1.968    0.524      0.939    1.968      3.001    1.965

# prefMean = 2.5:
# Fixed effects:
#                 mean    sd 0.025quant 0.5quant 0.975quant    mode kld
# Intercept_pp -16.984 0.643    -18.334  -16.955    -15.799 -16.857   0
# X_pp           1.963 0.395      1.210    1.954      2.761   1.954   0
# Intercept_y   -0.081 0.280     -0.659   -0.073      0.448  -0.073   0
# X_y            0.823 0.161      0.510    0.821      1.142   0.821   0
# 
# Random effects:
#   Name	  Model
# field_y SPDE2 model
# field_pp Copy
# 
# Model hyperparameters:
#                                                mean       sd 0.025quant 0.5quant 0.975quant     mode
# Precision for the Gaussian observations[2]  100.020    3.164     93.945   99.967    106.401   99.854
# Range for field_y                          3070.845 1051.629   1506.376 2905.781   5594.439 2601.978
# Stdev for field_y                             0.661    0.144      0.424    0.645      0.987    0.615
# Beta for field_pp                             2.145    0.544      1.115    2.132      3.253    2.074


# plot predictions ----

preds = out$predEst
predAggMat = out$predAggMat
predAgg = out$predAggEst
predAggU = out$predAggUpper
predAggL = out$predAggLower
trueAgg = mean(gTruth)
predAgg
trueAgg

pdf(file=paste0(figDir, "testing/testPreds_Diggle_truthPref.pdf"), width=8, height=8)
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
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Diggle)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_Diggle_truthPref.pdf"), width=5, height=5)
par(mfrow=c(1,1), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (Diggle)", 
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

# fit model ----
profvis({out <- fitDigglesimDat(wellTestDat_successiveIPP, seismicTestDat_successiveIPP)})
system.time(out <- fitDigglesimDat(wellTestDat_successiveIPP, seismicTestDat_successiveIPP))
# n=3500, max.n=5000: 288.56 
# n=1500, max.n=2500: 71.86  

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

pdf(file=paste0(figDir, "testing/testPreds_Diggle_successiveIPP.pdf"), width=8, height=8)
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
       zlim=c(0, 1), xlab="Easting", ylab="Northing", main="Estimate (Diggle)", 
       asp=1, smallplot=c(.83,.87,.25,.8))
points(pEast, pNorth, cex=.5)
dev.off()

pdf(file=paste0(figDir, "testing/testAgg_Diggle_successiveIPP.pdf"), width=5, height=5)
par(mfrow=c(1,1), oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 5.5))
hist(predAggMat, breaks=30, col="skyblue", main="Posterior for Volume Fraction (Diggle)", 
     xlab="Volume Fraction", freq=F, xlim=c(0,1))
abline(v=c(predAggL, predAggU), col="purple", lty=2)
abline(v=predAgg, col="purple")
abline(v=trueAgg, col="black")
abline(v=svyAggEst, col="blue")
legend("topright", legend=c("Mean", "80% CI", "Truth", "HT"), lty=c(1, 2, 1, 1), col=c("purple", "purple", "black", "blue"))
dev.off()






# ~~~~~~~~~~~~~~ ----
# n=250 truth pref ----
# ~~~~~~~~~~~~~~ ----

testDat = getTestDat(n=250, prefType="truth", saveDat=FALSE)

wellDat = testDat$wellDat
seismicDat = testDat$seismicDat
truthDat = testDat$truthDat

out = fitDigglesimDat(wellDat, seismicDat, 
                      mesh=getSPDEmeshTestDat(), prefMean=2.5)

summary(out$mod)

# prefMean = 2.5
# Fixed effects:
#                 mean    sd 0.025quant 0.5quant 0.975quant    mode kld
# Intercept_pp -14.979 0.534    -16.062  -14.969    -13.957 -14.969   0
# X_pp           2.373 0.146      2.089    2.372      2.662   2.372   0
# Intercept_y   -0.354 0.241     -0.841   -0.350      0.107  -0.350   0
# X_y            1.123 0.058      1.010    1.122      1.236   1.122   0
# 
# Random effects:
#   Name	  Model
# field_y SPDE2 model
# field_pp Copy
# 
# Model hyperparameters:
#                                                mean      sd 0.025quant 0.5quant 0.975quant     mode
# Precision for the Gaussian observations[2]   97.969   3.091     92.062   97.907     104.23   97.757
# Range for field_y                          2736.558 418.489   2027.538 2697.846    3671.08 2609.116
# Stdev for field_y                             0.814   0.106      0.632    0.805       1.05    0.785
# Beta for field_pp                             2.205   0.138      1.935    2.205       2.48    2.203





















