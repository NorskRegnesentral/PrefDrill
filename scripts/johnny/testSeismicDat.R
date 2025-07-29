# script for checking seismic data and truth

# Ragnar's datasets ----
out = readSurfaceRMS("K:/Internal/190693_GP_PrefDrill/07_Work/synthetic_model/RegularizedPred.txt")
regPred = out$surfMat
regPredFrame = out$surfFrame
image(regPred)
dim(regPred)

squilt(x=regPredFrame$east, y=regPredFrame$north, z=regPredFrame$seismicEst)

out = readSurfaceRMS("K:/Internal/190693_GP_PrefDrill/07_Work/synthetic_model/RegularizedSand.txt")
regSand = out$surfMat
image(regSand)
dim(regSand)

cor(c(regSand), c(regPred))

out = readSurfaceRMS("K:/Internal/190693_GP_PrefDrill/07_Work/synthetic_model/RegularizedTrue.txt")
regTrue = out$surfMat
image(regTrue)
dim(regTrue)

cor(c(regSand), c(regTrue))

# Dataset I generated
out = readSurfaceRMS("C:/Users/jpaige/Desktop/synthetic_model/RegularizedSand.txt")
regSand = out$surfMat
image(regSand)
dim(regSand)

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



# Sim study datasets ----

# Ragnar's datasets
out = readSurfaceRMS("C:/Users/jpaige/Desktop/synthetic_model/RegularizedPred_2.txt")
regPred = out$surfMat
regPredFrame = out$surfFrame
image(regPred)
dim(regPred)

# squilt(x=regPredFrame$east, y=regPredFrame$north, z=regPredFrame$seismicEst)

out = readSurfaceRMS("C:/Users/jpaige/Desktop/synthetic_model/RegularizedSand_2.txt")
regSand = out$surfMat
image(regSand)
dim(regSand)

cor(c(regPred), c(regSand))

out = readSurfaceRMS("C:/Users/jpaige/Desktop/synthetic_model/RegularizedPred_1.txt")
regPred = out$surfMat
image(regPred)
dim(regPred)

cor(c(regSand), c(regPred))

# Dataset I generated
out = readSurfaceRMS("C:/Users/jpaige/Desktop/synthetic_model/RegularizedSand_1.txt")
regSand = out$surfMat
image(regSand)
dim(regSand)

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


root = "C:/Users/jpaige/Desktop/synthetic_model/"
root = "~/synthetic_model/"
allCors = numeric(140)
for(i in 1:140) {
  print(paste0("i: ", i, "/140"))
  out = readSurfaceRMS(paste0(root, "RegularizedPred_", i, ".txt"))
  regPred = out$surfMat
  regPredFrame = out$surfFrame
  image(regPred)
  dim(regPred)
  
  # squilt(x=regPredFrame$east, y=regPredFrame$north, z=regPredFrame$seismicEst)
  
  out = readSurfaceRMS(paste0(root, "RegularizedSand_", i, ".txt"))
  regSand = out$surfMat
  image(regSand)
  dim(regSand)
  
  allCors[i] = cor(c(regPred), c(regSand))
}

fixedAllCors = allCors[-c(53, 28, 24, 78, 22)]

fnameTab = data.frame(oldI=1:length(allCors), goodL=allCors>0.5)
newI = rep(NA, nrow(fnameTab))
newI[fnameTab$goodL] = 1:sum(fnameTab$goodL)
newI[newI > 100] = NA
fnameTab$newI = newI

fnameTab$goodL = !is.na(fnameTab$newI)

# Loop through files, move them to the new directory only if their correlation 
# is above 0.5
for (r in 1:nrow(fnameTab)) {
  old_i <- fnameTab$oldI[r]
  new_i <- fnameTab$newI[r]
  
  if (!is.na(new_i)) {
    old_name <- paste0("~/synthetic_model/RegularizedSand_", old_i, ".txt")
    new_name <- paste0("~/git/PrefDrill/data/seisTruthReplicates/RegularizedSand_", new_i, ".txt")
    
    file.copy(old_name, new_name)
    
    old_name <- paste0("~/synthetic_model/RegularizedPred_", old_i, ".txt")
    new_name <- paste0("~/git/PrefDrill/data/seisTruthReplicates/RegularizedPred_", new_i, ".txt")
    
    file.copy(old_name, new_name)
  } else {
    message(sprintf("Skipping oldI = %d (marked as bad)", old_i))
  }
}

pdf("figures/testSimStudy/corHist.pdf", width=5, height=5)
hist(allCors, main="Cor(regPred, regSand)", breaks=20, xlab="Correlation")
dev.off()

pdf("figures/testSimStudy/corHistFixed.pdf", width=5, height=5)
hist(fixedAllCors, main="Cor(regPred, regSand)", breaks=20, xlab="Correlation")
dev.off()

pdf("figures/testSimStudy/corECDF.pdf", width=5, height=5)
plot(ecdf(allCors), main="Cor(regPred, regSand)", xlab="Correlation", ylab="ECDF")
dev.off()

pdf("figures/testSimStudy/corECDFfixed.pdf", width=5, height=5)
plot(ecdf(fixedAllCors), main="Cor(regPred, regSand)", xlab="Correlation", ylab="ECDF")
dev.off()

tail(sort(allCors))
head(allCors)

image(regPred)
image(regSand)

# CRAVA error runs:
# jpaige@loginsamba:~/synthetic_model$ ls -lt $(grep -l "Error with" *.log)
# -rw-rw-r-- 1 jpaige u_sand 15451 Jun 26 14:58 invseisLog_53.log
# -rw-rw-r-- 1 jpaige u_sand 15451 Jun 26 01:32 invseisLog_28.log
# -rw-rw-r-- 1 jpaige u_sand 15530 Jun 25 23:27 invseisLog_24.log

# CRAVA crashes:
# jpaige@loginsamba:~/synthetic_model$ ls -lt $(grep -L "CRAVA finished" invseis*.log)
# -rw-rw-r-- 1 jpaige u_sand 18947 Jun 27 04:36 invseisLog_78.log
# -rw-rw-r-- 1 jpaige u_sand 15451 Jun 26 14:58 invseisLog_53.log
# -rw-rw-r-- 1 jpaige u_sand 18946 Jun 26 04:12 invseisLog_33.log
# -rw-rw-r-- 1 jpaige u_sand 15451 Jun 26 01:32 invseisLog_28.log
# -rw-rw-r-- 1 jpaige u_sand 15530 Jun 25 23:27 invseisLog_24.log






