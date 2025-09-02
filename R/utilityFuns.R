# Script with miscellaneous useful functions

logit <- function(x) {
  log(x/(1-x))
}

# aka the inverse logit function
expit <- function(x) {
  res = 1/(1+exp(-x))
  res[x >= 37] = 1
  res[x <= -101] = 0
  res
}

# Use bilinear interpolation on a grid, where gridDat has long format
# 
# inputs:
# pts: easting/northing coordinates of pts to extract bilinearly interpolated 
#      value for based on gridDat
# gridDat: 3 columns, with first two being names "east" and "north", third being 
#          what to interpolate
# transform: transformation to third column of gridDat prior to interpolation
# invTransform: inverse of transform (i.e. back transformation)
# 
# outputs: vector of interpolated values at input pts
bilinearInterp = function(pts, gridDat, transform=I, invTransform=I, safeExtrap=FALSE) {
  require(akima)
  
  # use bilinear interpolation to get seismic estimate at that point
  # First must convert to akima format into matrix grid form for z
  eastGrid = sort(unique(gridDat[,1]))
  northGrid = sort(unique(gridDat[,2]))
  sortI1 = order(gridDat[,1]) # east
  gridDat = gridDat[sortI1,]
  sortI2 = order(gridDat[,2]) # north
  gridDat = gridDat[sortI2,]
  
  z = matrix(transform(gridDat[,3]), nrow=length(eastGrid), ncol=length(northGrid))
  
  if(safeExtrap) {
    # if necessary, make sure extrapolated pts get closest value
    lowEastL = pts[,1] < eastGrid[1]
    highEastL = pts[,1] > eastGrid[length(eastGrid)]
    lowNorthL = pts[,2] < northGrid[1]
    highNorthL = pts[,2] > northGrid[length(northGrid)]
    
    if(any(lowEastL)) {
      
    }
    if(any(highEastL)) {
      
    }
    if(any(lowNorthL)) {
      
    }
    if(any(highNorthL)) {
      
    }
  }
  
  seismicEst = invTransform(akima::bilinear(x=eastGrid, y=northGrid, z=z, x0 = pts[,1], y0=pts[,2])$z)
  seismicEst
}

# adapted from logitnorm package.  Calculates the mean of a distribution whose 
# logit is Gaussian. Each row of muSigmaMat has a mean and standard deviation 
# on the logit scale
logitNormMean = function(muSigmaMat, parClust=NULL, logisticApproximation=FALSE, splineApproximation=FALSE, ...) {
  if(length(muSigmaMat) > 2) {
    if(is.null(parClust) && !splineApproximation) {
      apply(muSigmaMat, 1, logitNormMean, logisticApproximation=logisticApproximation, ...)
    }
    else if(splineApproximation) {
      # not parallel and using spline approximation
      uniqueSigmas = sort(unique(muSigmaMat[,2]))
      muSigmaIndsList = lapply(uniqueSigmas, function(s) {
        inds = which(muSigmaMat[,2] == s)
        list(mu=muSigmaMat[inds,1], inds=inds, sigma=s)
      })
      outList = lapply(muSigmaIndsList, function(l) {
        mus = l$mu
        sigma = l$sigma
        logitNormMeanSplineApprox(mus, sigma, ...)
      })
      outVals = 1:nrow(muSigmaMat)
      for(i in 1:length(muSigmaIndsList)) {
        inds = muSigmaIndsList[[i]]$inds
        outVals[inds] = outList[[i]]$vals
      }
      outVals
    } else {
      parApply(parClust, muSigmaMat, 1, logitNormMean, logisticApproximation=logisticApproximation, ...)
    }
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    if(sigma == 0)
      expit(mu)
    else {
      if(any(is.na(c(mu, sigma))))
        NA
      else if(!logisticApproximation) {
        # numerically calculate the mean
        fExp <- function(x) exp(plogis(x, log.p=TRUE) + dnorm(x, mean = mu, sd = sigma, log=TRUE))
        outVal = integrate(fExp, mu-10*sigma, mu+10*sigma, abs.tol = 0, ...)$value
        
        if(is.na(outVal)) {
          # if the result is NA, such as in the case of numerical instability, just use logistic approximation
          k = 16 * sqrt(3) / (15 * pi)
          expit(mu / sqrt(1 + k^2 * sigma^2))
          outVal = expit(mu)
        }
        
        outVal
      } else {
        # use logistic approximation
        k = 16 * sqrt(3) / (15 * pi)
        expit(mu / sqrt(1 + k^2 * sigma^2))
      }
    }
  }
}

# Similar to logitNormMean, except assumes one sigma at the top of each column
# of muSigmaMat
logitNormMeanGrouped = function(muSigmaMat, logisticApproximation=TRUE, splineApproximation=FALSE, parClust=NULL, ...) {
  if(is.matrix(muSigmaMat) && ncol(muSigmaMat) > 1) {
    if(is.null(parClust)) {
      apply(muSigmaMat, 2, logitNormMeanGrouped, 
            logisticApproximation=logisticApproximation, 
            splineApproximation=splineApproximation, ...)
    } else {
      parApply(parClust, muSigmaMat, 2, logitNormMeanGrouped, 
               logisticApproximation=logisticApproximation, 
               splineApproximation=splineApproximation, ...)
    }
  } else {
    sigma = muSigmaMat[1]
    mus = muSigmaMat[-1]
    
    if(is.na(sigma)) {
      rep(NA, length(mus))
    } else if(sigma == 0) {
      expit(mus)
    } else if (splineApproximation) {
      logitNormMeanSplineApprox(mus, sigma, ...)$vals
    } else {
      logitNormMean(muSigmaMat=cbind(mus, sigma), parClust=parClust, 
                    logisticApproximation=logisticApproximation, 
                    splineApproximation=splineApproximation, ...)
    }
  }
}

# Approximates logitNormMean at a single value of sigma and many mus using spline 
# on a logit scale. npts determines number of values of mu in the range of mu 
# over which the monotonic spline function is generated.
# Note: Uses a monotonic cubic spline. See ?splinefun for method="Hyman", and:
# Hyman, J. M. (1983). Accurate monotonicity preserving cubic interpolation. 
# SIAM Journal on Scientific and Statistical Computing, 4, 645–654. 
# doi:10.1137/0904045.
logitNormMeanSplineApprox = function(mus, sigma, npts=250, ...) {
  
  rangeMu = range(mus)
  rangeExpitMu = expit(rangeMu)
  if(!any(rangeExpitMu %in% c(0, 1))) {
    # all values of mu are in a good range for numerical evaluation. 
    # Construct evaluation grid on probability scale
    seqMus = logit(seq(rangeExpitMu[1], rangeExpitMu[2], l=npts))
  } else {
    # some values of mu outside of a good range for numerical evaluation. 
    # Stitch together evaluation grids on probability and logit scales. 
    if(all(rangeExpitMu %in% c(0, 1))) {
      if(all(rangeExpitMu == 0)) {
        seqMus = seq(rangeMu[1], rangeMu[2], l=npts)
      } else if(all(rangeExpitMu == 1)) {
        seqMus = seq(rangeMu[1], rangeMu[2], l=npts)
      } else if(all(rangeExpitMu %in% c(0, 1))) {
        nPtsBelow = nPtsAbove = ceiling(npts/10)
        nPtsGood = round(npts*.8)
        goodRange = c(-709, 36)
        seqMus = c(seq(rangeMu[1], goodRange[1], l=nPtsBelow), 
                   logit(seq(expit(goodRange[1]), expit(goodRange[2]), l=nPtsGood)), 
                   seq(goodRange[2], rangeMu[2], l=nPtsAbove))
      }
    } else if(rangeExpitMu[1] == 0) {
      nPtsBelow = ceiling(npts/10)
      nPtsGood = round(npts*.9)
      goodRange = c(-709, 36)
      seqMus = c(seq(rangeMu[1], goodRange[1], l=nPtsBelow), 
                 logit(seq(expit(goodRange[1]), expit(goodRange[2]), l=nPtsGood)))
    } else if(rangeExpitMu[2] == 1) {
      nPtsAbove = ceiling(npts/10)
      nPtsGood = round(npts*.9)
      goodRange = c(-709, 36)
      seqMus = c(logit(seq(expit(goodRange[1]), expit(goodRange[2]), l=nPtsGood)), 
                 seq(goodRange[2], rangeMu[2], l=nPtsAbove))
    }
  }
  
  
  muSigmaMat = cbind(seqMus, sigma)
  vals = logit(logitNormMean(muSigmaMat, logisticApproximation=FALSE, splineApproximation=FALSE))
  
  spFun = splinefun(seqMus, vals)
  
  outVals = expit(spFun(mus))
  
  list(vals = outVals, fun=spFun, range=rangeMu)
}

# converts data.frame to a list of lists
dfToListOfLists = function(df) {
  list_of_lists <- split(df, seq(nrow(df)))
  lapply(list_of_lists, as.list)
}

# converts the first letter to uppercase, keeps the case of all other letters
myTitleCase <- function(x) {
  if (nchar(x) == 0) return(x)
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}


cppFunction('
NumericMatrix removeOnePerRow(NumericMatrix mat, IntegerVector tIs) {
  int nr = mat.nrow();
  int nc = mat.ncol();
  NumericMatrix out(nr, nc - 1);

  for (int i = 0; i < nr; ++i) {
    int skip_col = tIs[i] - 1; // R is 1-based, C++ is 0-based
    int out_col = 0;

    for (int j = 0; j < nc; ++j) {
      if (j != skip_col) {
        out(i, out_col) = mat(i, j);
        out_col++;
      }
    }
  }

  return out;
}
')

cppFunction('
List insertSortedColumnWithIndices(NumericMatrix mat, NumericVector colvec) {
  int nr = mat.nrow();
  int nc = mat.ncol();
  NumericMatrix out(nr, nc + 1);
  IntegerVector insert_pos(nr);

  for (int i = 0; i < nr; ++i) {
    double new_val = colvec[i];
    int j = 0;

    // Find insert position in sorted row
    while (j < nc && mat(i, j) < new_val) {
      out(i, j) = mat(i, j);
      j++;
    }

    out(i, j) = new_val;
    insert_pos[i] = j + 1; // R is 1-based

    // Copy the rest of the row
    for (int k = j; k < nc; ++k) {
      out(i, k + 1) = mat(i, k);
    }
  }

  return List::create(Named("newMat") = out,
                      Named("insertIndices") = insert_pos);
}
')

cppFunction('
NumericMatrix sweep_col_mult(const NumericMatrix& mat, const NumericVector& vec) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  
  if (vec.size() != ncol) {
    stop("Length of vec must be equal to number of columns in mat");
  }
  
  NumericMatrix out(nrow, ncol);
  
  for (int j = 0; j < ncol; j++) {
    double multiplier = vec[j];
    for (int i = 0; i < nrow; i++) {
      out(i, j) = mat(i, j) * multiplier;
    }
  }
  
  return out;
}
')
