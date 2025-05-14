# script for checking seismic data and truth

out = readSurfaceRMS("K:/Internal/190693_GP_PrefDrill/07_Work/synthetic_model/RegularizedPred.txt")
regPred = out$surfMat
image(regPred)
dim(regPred)

out = readSurfaceRMS("K:/Internal/190693_GP_PrefDrill/07_Work/synthetic_model/RegularizedSand.txt")
regSand = out$surfMat
image(regSand)
dim(regSand)

out = readSurfaceRMS("K:/Internal/190693_GP_PrefDrill/07_Work/synthetic_model/RegularizedTrue.txt")
regTrue = out$surfMat
image(regTrue)
dim(regTrue)

# Sand and true seem identical and not very relative to pred...
