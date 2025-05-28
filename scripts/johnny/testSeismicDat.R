# script for checking seismic data and truth

out = readSurfaceRMS("K:/Internal/190693_GP_PrefDrill/07_Work/synthetic_model/RegularizedPred.txt")
regPred = out$surfMat
image(regPred)
dim(regPred)

out = readSurfaceRMS("K:/Internal/190693_GP_PrefDrill/07_Work/synthetic_model/RegularizedSand.txt")
regSand = out$surfMat
image(regSand)
dim(regSand)

cor(c(regSand), c(regPred))

out = readSurfaceRMS("K:/Internal/190693_GP_PrefDrill/07_Work/synthetic_model/RegularizedTrue.txt")
regTrue = out$surfMat
image(regTrue)
dim(regTrue)

# Sand and true seem identical and not very relative to pred...

names(out)
range(out$xgrid)
range(out$ygrid)
dx = diff(out$xgrid)[1]
dy = diff(out$ygrid)[1]
xr = c(min(out$xgrid) - dx/2, max(out$xgrid) + dx/2)
yr = c(min(out$ygrid) - dy/2, max(out$ygrid) + dy/2)
xr
yr

totArea = diff(xr) * diff(yr)
repAreas = c(.0025, .005, .01) * totArea
repDists = sqrt(repAreas/pi)
repDists
repDists^2 * pi
repAreas

simStudyMesh = getSPDEmeshSimStudy()
plot(simStudyMesh)
tempPrior = getSPDEprior(simStudyMesh)

plot(regSand, regPred, pch=".")
abline(0, 1, col="blue")

# 





