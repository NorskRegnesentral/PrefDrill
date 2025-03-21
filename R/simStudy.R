# functions for the simulation study

# construct a rectangular prediction grid over a domain
# xLims: a 2-vector with the min and max x value
# yLims: a 2-vector with the min and max y value
getSimStudyPredGrid = function(xLims = simStudyXlims, yLims = simStudyXlims) {
  make.surface.grid(list(x=xLims, y=yLims))
}

# samples the next drilling location.
# 
# If betaP == gammaP == 0, sampling is uniform, 
# otherwise with prob proportional to the exponential of seismic and spatial data 
# times betaP and gammaP respectively. Regardless, samples uniformly within 
# selected grid cell. Calculates seismic estimate at sampled location based on 
# bilinear interpolation.
# 
# Inputs:
# seismicDat: data.frame with columns:
#   east
#   north
#   seismicEst: estimate of sand volume fraction
# spatDat: data.frame with columns:
#   east
#   north
#   spatEst: estimate of the spatial component of the model. If 0, only seismic 
#            data is used to select the next drill location
# betaP: preferentiality parameter with respect to seismic estimates
# gammaP: preferentiality parameter with respect to spatial estimates
# seed: if not NULL, sets random number seed
# 
# Output: 
# a 3-vector with values: east, north, seismicEst
getPrefDrillLoc = function(seismicDat, 
                           spatDat=cbind(seismicDat$east, seismicDat$north, spatEst=0), 
                           betaP=0, gammaP=0, seed=NULL) {
  
  # set random number seed
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # calculate grid width, delta (assume north and east grids are equal along with delta)
  eastGrid = sort(unique(seismicDat$east))
  delta = eastGrid[2] - eastGrid[1]
  
  # calculate linear predictor, sampling probabilities
  etaP = logit(seismicDat$seismicEst) * betaP + spatDat$spatEst * gammaP
  sampleProbs = exp(etaP)
  sampleProbs = sampleProbs/sum(sampleProbs)
  
  # sample centroid location
  sampleI = sample(1:nrow(seismicDat), 1, replace=FALSE, prob = sampleProbs)
  centEast = seismicDat$east[sampleI]
  centNorth = seismicDat$north[sampleI]
  
  # sample point location uniformly within grid cell
  east = centEast = runif(1) * delta - delta/2
  north = centNorth = runif(1) * delta - delta/2
  
  # use bilinear interpolation to get seismic estimate at that point
  seismicEst = bilinearInterp(rbind(c(east, north)), seismicDat)
  
  # return sampled drill location and associated seismic estimate
  c(east=east, north=north, seismicEst=seismicEst)
}

# Inputs:
# wellDat: A data.frame with columns:
#   east: easting
#   north: northing
#   volFrac: sand volume fraction (nett og gross)
# seismicDat: a data.frame containing a grid of of seismic estimates with columns:
#   east: easting
#   north: northing
#   est: central estimate
# fitModFun: Function for fitting a model and generating predictions. Contains the 
#            following inputs, plus potentially more: wellDat, seismicDat.
# 
# Outputs:
# preds: prediction grid with central est plus SDs or grid of draws or just pred 
#        of total nett og gross estimate over entire domain???
#   
fitModel = function(wellDat, seismicDat, fitModFun, ...) {
  
}






