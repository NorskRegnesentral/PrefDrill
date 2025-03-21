# function for fitting densities and intensity functions (with no response 
# model), and obtaining component due to well data alone

getDensityDueToWells = function(seismicDat, wellDat, 
                                method=c("kde", "inlabru"), 
                                kde.args=NULL, lm.args=NULL, 
                                centerScale=TRUE, 
                                ...) {
  require(ks)
  require(inlabru)
  
  method = match.arg(method)
  
  predGrid = seismicDat[,1:2]
  pts = wellDat[,1:2]
  
  if(method == "kde") {
    # first fit the kernel density smoother
    densEst = do.call("kde", c(list(pts, eval.points=predGrid), kde.args))$estimate
    
    browser() # make sure densEst is only densities not coords
    
    # normalize to integrate to 1 over domain
    densEst = densEst * (1/sum(densEst))
    
    # now subtract out the part of the density due to seismic data
    logDens = log(densEst)
    logitSeismicEst = logit(seismicDat$seismicEst)
    mod = do.call("lm", c(list(logDens~logitSeismicEst), lm.args))
    resids = residuals(mod)
    
    wellEffect = resids
  } else if(method == "inlabru") {
    stop("inlabru not yet supported")
    
    
  }
  
  if(centerScale) {
    # if requested by user, center and scale the estimate of well preferential 
    # effect
    wellEffect = wellEffect - mean(wellEffect)
    wellEffect = wellEffect * (1/sd(wellEffect))
  }
  
  # use bilinear interpolation to estimate wellEffect at well points
  wellEffectWells = bilinearInterp(pts, cbind(predGrid, wellEffect))
  
  list(wellEffectGrid=wellEffect, wellEffectPts=wellEffectWells)
}








