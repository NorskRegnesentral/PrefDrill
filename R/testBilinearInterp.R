# test script for the bilinearInterp function

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
splot(seismicTestDat$east, seismicTestDat$north, seismicTestDat$seismicEst)
pSeismic = bilinearInterp(cbind(pEast, pNorth), seismicTestDat, 
                          transform=logit, invTransform = expit)
splot(gEast, gNorth, gSeismic, zlim=c(0,1))
splot(pEast, pNorth, pSeismic, xlim=c(0,10000), ylim=c(0,10000), zlim=c(0,1))

gEast = seq(0, 2, l=100)
gNorth = seq(0, 1, l=30)
pEast = runif(100)*2
pNorth = runif(100)
gridCoords = make.surface.grid(list(east=gEast, north=gNorth))
zs = sin(gridCoords[,1]*pi)
gridDat = data.frame(east=gridCoords[,1], north=gridCoords[,2], z=zs)
pSeismic = bilinearInterp(cbind(pEast, pNorth), gridDat)
splot(gridDat$east, gridDat$north, gridDat$z, xlim=c(0,2), ylim=c(0,1), zlim=c(-1,1))
splot(pEast, pNorth, pSeismic, xlim=c(0,2), ylim=c(0,1), zlim=c(-1,1))




