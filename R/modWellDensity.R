# function for fitting densities and intensity functions (with no response 
# model), and obtaining component due to well data alone

# estimate log density of well points after adjusting/controlling for seismic 
# data. 
# 
# Inputs:
# seismicDat: seismic data on a grid (see modSPDE)
# wellDat: well data at point locations (see modSPDE)
# predPts: points where final predictions will be made/final estimates are 
#          desired
# method: if kde, estimates will be based on residuals of regression of log 
#         densities as estimated by kernel density estimator against logit 
#         seismic estimates. If inlabru (not currently supported), estimates 
#         will be based on point process model including logit seismic estimates 
#         as a covariate and an additional spatial effect. The spatial effect 
#         represents effects due to wells but not seismic data.
# kde.args: arguments to ks::kde function
# lm.args: arguments to lm
# centerScale: whether or not to centerScale the final set of estimates to have 
#              mean 0 and sd 1 on seismic data grid
getDensityDueToWells = function(seismicDat, wellDat, predPts=NULL, 
                                method=c("kde", "inlabru"), 
                                kde.args=NULL, lm.args=NULL, 
                                centerScale=TRUE) {
  require(ks)
  require(inlabru)
  
  method = match.arg(method)
  
  predGrid = seismicDat[,1:2]
  pts = wellDat[,1:2]
  
  if(method == "kde") {
    # first fit the kernel density smoother
    densEst = do.call("kde", c(list(pts, eval.points=predGrid, density=TRUE), kde.args))$estimate
    
    # normalize to integrate to 1 over domain
    densEst = densEst * (1/sum(densEst))
    log1pDens = log1p(densEst) # would do log, but sometimes density estimates are numerically 0
    
    # now subtract out the part of the density due to seismic data
    # logitSeismicEst = logit(seismicDat$seismicEst)
    # mod = do.call("lm", c(list(logDens~logitSeismicEst), lm.args))
    # resids = residuals(mod)
    # 
    # wellEffect = resids
    
    wellEffect = log1pDens
  } else if(method == "inlabru") {
    stop("inlabru not yet supported")
    
  }
  
  if(centerScale) {
    # if requested by user, center and scale the estimate of well preferential 
    # effect
    wellEffect = wellEffect - mean(wellEffect)
    wellEffect = wellEffect * (1/sd(wellEffect))
  }
  
  if(!is.null(predPts)) {
    # use bilinear interpolation to estimate wellEffect at custom points
    wellEffectPreds = bilinearInterp(predPts, cbind(predGrid, wellEffect))
  } else {
    # just return predictions over the grid
    wellEffectPreds = NULL
  }
  
  
  list(wellEffectGrid=wellEffect, wellEffectPredPts=wellEffectPreds)
}








