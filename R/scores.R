# this script contains all used scoring rules used to evaluate the models for a given data set

# this function computes all scoring rules
# truth: the true values
# est: the estimates, of the same length as truth. By default calculated from estMat
# var: the estimates, of the same length as truth. By default calculated from estMat
# lower: the lower end of the credible interval. By default calculated from estMat
# upper: the upper end of the credible interval. By default calculated from estMat
# estMat: a matrix of joint estimate draws, with number of rows equal to the length of truth, a number of 
#          columns equal to the number of draws. If not included, a gaussian distribution is assumed.
# significance: the significance level of the credible interval. By default 80%
# distances: the distances to the nearest observation if not NULL, scores are broken up 
#            as a function of nearest neighbor distances
# breaks: the number of equal spaced bins to break the scores into as a function of distance
# NOTE: Discrete, count level credible intervals are estimated based on the input estMat along with coverage and CRPS
getScores = function(truth, est=NULL, var=NULL, lower=NULL, upper=NULL, estMat=NULL, 
                     significance=c(.8, .95), distances=NULL, breaks=30, aggScores=TRUE, aggFun=robustWeightedMean, aggFunArgs=NULL, 
                     ns=NULL, anyNAisNA=FALSE, returnNAs=FALSE, na.rm=FALSE, setInfToNA=FALSE, throwOutAllNAs=FALSE) {
  
  if(setInfToNA) {
    naRows = !is.finite(truth)
    if(!is.null(est)) {
      naRows = naRows | !is.finite(est)
    }
    if(!is.null(var)) {
      naRows = naRows | !is.finite(var)
    }
    if(!is.null(lower)) {
      naRows = naRows | !is.finite(lower)
    }
    if(!is.null(upper)) {
      naRows = naRows | !is.finite(upper)
    }
    if(!is.null(estMat) && throwOutAllNAs) {
      naRows = naRows | apply(estMat, 1, function(x) {any(!is.finite(x))})
    }
    
    truth[naRows] = NA
    if(!is.null(est)) {
      est[naRows] = NA
    }
    if(!is.null(var)) {
      var[naRows] = NA
    }
    if(!is.null(lower)) {
      lower[naRows] = NA
    }
    if(!is.null(upper)) {
      upper[naRows] = NA
    }
    if(!is.null(estMat)) {
      estMat[!is.finite(estMat)] = NA
      estMat[naRows,] = NA
    }
  }
  
  # compute central estimates if estMat is not null
  if(!is.null(estMat)) {
    if(is.null(est))
      est = rowMeans(estMat, na.rm=na.rm)
  }
  
  # first calculate bias, variance, and MSE
  out = mse(truth, est, aggScores=aggScores, aggFun=aggFun, aggFunArgs=aggFunArgs)
  thisMSE = out$MSE
  thisBias = out$bias
  thisVar = out$var
  
  # calculate coverage and credible interval width with and without binomial variation
  intScore = intervalScore(truth, est, var, lower, upper, estMat=estMat, 
                           significance=significance, returnIntervalWidth=TRUE, 
                           returnCoverage=TRUE, aggFun=aggFun, aggFunArgs=aggFunArgs, ns=ns, 
                           doFuzzyReject=doFuzzyReject, aggScores=aggFunArgs, na.rm=na.rm)
  if(aggScores) {
    thisIntScore = intScore[grepl("intScore", names(intScore))]
    thisCoverage = intScore[grepl("coverage", names(intScore))]
    thisWidth = intScore[grepl("width", names(intScore))]
  } else {
    thisIntScore = intScore[,grepl("intScore", colnames(intScore))]
    thisCoverage = intScore[,grepl("coverage", colnames(intScore))]
    thisWidth = intScore[,grepl("width", colnames(intScore))]
  }
  
  # calculate CRPS
  thisCRPS = crps(truth, est, var, estMat=estMat, aggScores=aggScores, aggFun=aggFun, aggFunArgs=aggFunArgs, na.rm=na.rm)
  
  # collect the results in a data frame
  results = matrix(c(thisBias, thisVar, thisMSE, sqrt(thisMSE), thisCRPS, thisIntScore, thisCoverage, 
                     thisWidth), ncol=5 + 3*length(significance))
  colnames(results) = c("Bias", "Var", "MSE", "RMSE", "CRPS", 
                        paste("IntervalScore", 100*significance, sep=""), 
                        paste("Coverage", 100*significance, sep=""), 
                        paste("Width", 100*significance, sep=""))
  results = as.data.frame(results)
  
  results
}

# calculate bias, variance, and MSE
mse <- function(truth, est, aggScores=TRUE, aggFun=robustWeightedMean, aggFunArgs=NULL){
  res = est - truth
  
  if(aggScores) {
    MSE = do.call("aggFun", c(list(res^2), resFunArgs))
    bias = do.call("aggFun", c(list(res), resFunArgs))
    thisVar = do.call("aggFun", c(list((res - mean(res, na.rm=TRUE))^2), resFunArgs))
  } else {
    MSE = res^2
    bias = res
    thisVar = (res - mean(res, na.rm=TRUE))^2
  }
  out = list(MSE=MSE, bias=bias, var=thisVar)
  
  out
}

# either include both lower and upper, or include either: 
#    - the joint estimate draw matrix
#    - estimates and variances (assumes gaussian)
# truth: the true empirical proportions of mortality rates within the regions or enumeration areas of interest
# lower: the lower end of the credible interval
# upper: the upper end of the credible interval
# estMat: a matrix of joint draws of estimates, with number of rows equal to the length of truth, a number of 
#         columns equal to the number of draws. If not included, a lgaussian distribution is assumed. Can be 
#         Gaussian or discrete values such as empirical proportions
# significance: the significance level of the credible interval. By default 80%
# doFuzzyReject: based on https://www.jstor.org/stable/pdf/20061193.pdf
# ns: a vector of maximum possible counts (denominators) for each observation. Used only for random/fuzzy reject. 
#     Can be left out, in which case it will be inferred from the minimum draw difference in each row of estMat.
coverage = function(truth, est=NULL, var=NULL, lower=NULL, upper=NULL, 
                    estMat=NULL, significance=.8, returnIntervalWidth=FALSE, 
                    doFuzzyReject=TRUE, aggScores=TRUE, ns=NULL, aggFun=robustWeightedMean, 
                    aggFunArgs=NULL){
  
  # if more than 1 significance level, return results for each
  if(length(significance) > 1) {
    res = lapply(significance, coverage, truth=truth, est=est, var=var, lower=lower, upper=upper, 
                 estMat=estMat, returnIntervalWidth=returnIntervalWidth, 
                 doFuzzyReject=doFuzzyReject, aggScores=aggScores, ns=ns, 
                 aggFun=aggFun, aggFunArgs=aggFunArgs)
    out = do.call("cbind", res)
  }
  
  if(any(is.null(lower)) || any(is.null(upper))) {
    # if the user did not supply their own credible intervals, we must get them ourselves given the other information
    
    if(is.null(estMat) && (is.null(est) || is.null(var)))
      stop("either include both lower and upper, est and var, or estMat")
    
    if(!is.null(est) && !is.null(var) && is.null(estMat)) {
      # in this case, we must calculate lower and upper assuming gaussianity
      lower = qnorm((1 - significance) / 2, est, sqrt(var))
      upper = qnorm(1 - (1 - significance) / 2, est, sqrt(var))
    }
    else {
      # we don't have information about the predictive distribution, and don't assume normality. 
      # Instead, use the user supplied to probability matrix estMat
      
      # take the quantiles of the probability draws
      CIs = apply(estMat, 1, function(ps) {quantile(ps, probs=c((1 - significance) / 2, 1 - (1 - significance) / 2))})
      lower = CIs[1,]
      upper = CIs[2,]
    }
  }
  
  if(any(lower > upper, na.rm=na.rm)) {
    warning("lower > upper, reordering")
    tmp = lower
    wrongOrder = lower > upper
    lower[wrongOrder] = upper[wrongOrder]
    upper[wrongOrder] = tmp[wrongOrder]
  }
  
  res = lower <= truth & upper >= truth
  
  width = upper - lower
  
  if(doFuzzyReject) {
    # in this case, we fuzzy reject if the truth is at the edge of the coverage interval. First 
    # determine what values are at the edge of the intervals, then determine the probability of rejection 
    # for each, then subtract fuzzy rejection probability from the coverage for that observation at the edge
    atLowerEdge = lower == truth
    atUpperEdge = upper == truth
    lowerEdgeInds = which(atLowerEdge)
    upperEdgeInds = which(atUpperEdge)
    
    # Fuzzy rejection probabilities at the interval boundaries. Nonzero only at 
    # respective interval boundaries
    rejectLower = rep(0, length(truth))
    rejectUpper = rep(0, length(truth))
    if(length(lowerEdgeInds) != 0) {
      probRejectLower = sapply(lowerEdgeInds, function(i) {
        if(mean(estMat[i,] == upper[i]) != 0) {
          ((1 - significance) / 2 - mean(estMat[i,] > lower[i], na.rm=na.rm)) / mean(estMat[i,] == lower[i], na.rm=na.rm)
        } else {
          warning("lower end of CI not equal to any samples from estMat. Setting reject probability to 0.5")
          0.5
        }
      })
      rejectLower[lowerEdgeInds] = probRejectLower
    }
    if(length(upperEdgeInds) != 0) {
      probRejectUpper = sapply(upperEdgeInds, function(i) {
        if(mean(estMat[i,] == upper[i]) != 0) {
          ((1 - significance) / 2 - mean(estMat[i,] > upper[i], na.rm=na.rm)) / mean(estMat[i,] == upper[i], na.rm=na.rm)
        } else {
          warning("upper end of CI not equal to any samples from estMat. Setting reject probability to 0.5")
          0.5
        }
      })
      rejectUpper[upperEdgeInds] = probRejectUpper
    }
    
    # determine minimum differences between probabilities
    if(is.null(ns))
      deltas = apply(estMat, 1, function(x) {min(diff(sort(unique(x))), na.rm=na.rm)})
    else
      deltas = 1 / ns
    
    # reduce CI width based on fuzzy boundaries
    width = width - deltas*rejectLower - deltas*rejectUpper
    upper = upper - deltas*rejectUpper
    lower = lower + deltas*rejectLower
    
    if(length(lowerEdgeInds) != 0) {
      res[lowerEdgeInds] = 1-rejectLower[lowerEdgeInds]
    }
    if(length(upperEdgeInds) != 0) {
      res[upperEdgeInds] = 1-rejectUpper[upperEdgeInds]
    }
  }
  
  if(aggScores)
    allResults = c(coverage=do.call("aggFun", c(list(res, significance=significance), aggFunArgs)))
  else
    allResults = c(coverage=res)
  
  if(returnIntervalWidth) {
    if(getAverage)
      allResults = c(allResults, width=do.call("aggFun", c(list(width), aggFunArgs)))
    else
      allResults = cbind(allResults, width=width)
  }
  
  # adjust names to match the significance in case there are multiple significances
  if(is.matrix(allResults)) {
    colnames(allResults) = paste(colnames(allResults), 100*significance, sep="")
  } else {
    names(allResults) = paste(names(allResults), 100*significance, sep="")
  }
  
  allResults
}

# truth: a vector of observations on the desired scale
# est: a vector of logit-scale predictions of the same length as truth 
# my.var: a vector of logit-scale predictive variances of the same length as truth
# estMat: if available, use these probability draws in the integration. Use this argument 
#         when a gaussian approximation to the (possibly transformed) posterior is unreasonable
# aggScores: if FALSE, returns score for individual observations. Otherwise for all observations
crps <- function(truth, est=NULL, my.var=NULL, estMat=NULL, aggScores=TRUE, na.rm=FALSE, 
                 aggFun=robustWeightedMean, aggFunArgs=NULL){
  if(!is.null(est) && !is.null(my.var) && is.null(estMat)) {
    sig = sqrt(my.var)
    x0 <- (truth - est) / sig
    res <- sig * (1 / sqrt(pi) -  2 * dnorm(x0) - x0 * (2 * pnorm(x0) - 1))
    
    ## sign as in Held (2008)
    res <- -res
  }
  else {
    # Integrate numerically using estMat
    if(is.null(estMat))
      stop("must include either or both est and my.var, or estMat")
    
    crpsRow = function(rowI) {
      thisTruth = truth[rowI]
      
      # build the predictive cdf assuming from the empirical distribution given by 
      # estMat
      
      # thisCdf = ecdf(estMat[rowI,])
      sorted = estMat[rowI,] # already sorted
      
      if(na.rm) {
        sorted = sorted[!is.na(sorted)]
      }
      
      # since we are using the empirical distribution, there is a closed form for the integral
      allPoints = sort(c(sorted, thisTruth))
      deltas = diff(allPoints)
      firstGreater = match(TRUE, sorted >= thisTruth)
      vals = (1:length(sorted))/length(sorted)
      if(is.na(firstGreater))
        return(sum((vals)^2 * deltas, na.rm=na.rm))
      else if(firstGreater == 1)
        return(deltas[1] + sum((1-vals[1:(length(sorted)-1)])^2 * deltas[2:length(deltas)], na.rm=na.rm))
      else {
        left = sum(vals[1:(firstGreater-1)]^2 * deltas[1:(firstGreater-1)], na.rm=na.rm)
        mid = sum((1 - vals[firstGreater-1])^2 * deltas[firstGreater], na.rm=na.rm)
        right = ifelse(firstGreater == length(vals), 0, sum((1 - vals[firstGreater:(length(vals)-1)])^2 * deltas[(firstGreater+1):length(deltas)], na.rm=na.rm))
        return(left+mid+right)
      }
      
      # intFun = function(ws) {
      #   (thisCdf(ws) - (ws >= thisTruth))^2
      # }
    }
    
    if(!is.null(estMat))
      estMat = t(apply(estMat, 1, sortWithNAs))
    res = sapply(1:length(truth), crpsRow)
  }
  
  if(aggScores) {
    do.call("aggFun", c(list(res), aggFunArgs))
  } else {
    res
  }
}

# either include both lower and upper, or include either: 
#    - the joint estimate draw matrix
#    - estimates and variances (assumes gaussian)
# truth: the true empirical proportions of mortality rates within the regions or enumeration areas of interest
# lower: the lower end of the credible interval
# upper: the upper end of the credible interval
# estMat: a matrix of joint draws of estimates, with number of rows equal to the length of truth, a number of 
#         columns equal to the number of draws. If not included, a lgaussian distribution is assumed. Can be 
#         Gaussian or discrete values such as empirical proportions
# significance: the significance level of the credible interval. By default 80%
# doFuzzyReject: based on https://www.jstor.org/stable/pdf/20061193.pdf
# ns: a vector of maximum possible counts (denominators) for each observation. Used only for random/fuzzy reject. 
#     Can be left out, in which case it will be inferred from the minimum draw difference in each row of estMat.
# getAverage: if FALSE, returns score for individual observations. Otherwise for all observations
# NOTE: this does not account for fuzzy CIs for discrete data. Defines on p 13 of: 
# https://www.tandfonline.com/doi/pdf/10.1198/016214506000001437?casa_token=0vXXqMZ3M2IAAAAA:BYmw_z2zaASEcAvFrNDf6PQ157vq6FAQuDuI9depRZp44RJ_M8zbY47CN_KGXHMXP9CHJL02bTDT
intervalScore = function(truth, est=NULL, var=NULL, lower=NULL, upper=NULL, 
                         estMat=NULL, significance=.8, returnIntervalWidth=FALSE, 
                         returnCoverage=FALSE, doFuzzyReject=TRUE, aggScores=TRUE, ns=NULL, 
                         na.rm=FALSE, aggFun=robustWeightedMean, aggFunArgs=NULL){
  
  # if more than 1 significance level, return results for each
  if(length(significance) > 1) {
    res = lapply(significance, intervalScore, truth=truth, est=est, var=var, lower=lower, upper=upper, 
                 estMat=estMat, returnIntervalWidth=returnIntervalWidth, returnCoverage=returnCoverage, 
                 doFuzzyReject=doFuzzyReject, aggScores=aggScores, ns=ns, na.rm=na.rm, aggFun=aggFun, aggFunArgs=aggFunArgs)
    if(aggScores) {
      return(unlist(res))
    } else {
      return(do.call("cbind", res))
    }
  }
  
  if(any(is.null(lower)) || any(is.null(upper))) {
    # if the user did not supply their own credible intervals, we must get them ourselves given the other information
    
    if(is.null(estMat) && (is.null(est) || is.null(var)))
      stop("either include both lower and upper, est and var, or estMat")
    
    if(!is.null(est) && !is.null(var) && is.null(estMat)) {
      # in this case, we must calculate lower and upper assuming gaussianity
      lower = qnorm((1 - significance) / 2, est, sqrt(var))
      upper = qnorm(1 - (1 - significance) / 2, est, sqrt(var))
    }
    else {
      # we don't have information about the predictive distribution, and don't assume normality. 
      # Instead, use the user supplied to probability matrix estMat
      
      # take the quantiles of the probability draws
      CIs = apply(estMat, 1, function(ps) {quantile(ps, probs=c((1 - significance) / 2, 1 - (1 - significance) / 2), na.rm=na.rm)})
      lower = CIs[1,]
      upper = CIs[2,]
    }
  }
  
  if(any(lower > upper, na.rm=na.rm)) {
    warning("lower > upper, reordering")
    tmp = lower
    wrongOrder = lower > upper
    lower[wrongOrder] = upper[wrongOrder]
    upper[wrongOrder] = tmp[wrongOrder]
  }
  
  greaterThanLower = lower <= truth
  lessThanUpper = upper >= truth
  if(returnCoverage) {
    cvg = greaterThanLower & lessThanUpper
  }
  
  width = upper - lower
  
  if(doFuzzyReject) {
    # in this case, we fuzzy reject if the truth is at the edge of the coverage interval. First 
    # determine what values are at the edge of the intervals, then determine the probability of rejection 
    # for each, then subtract fuzzy rejection probability from the coverage for that observation at the edge
    atLowerEdge = lower == truth
    atUpperEdge = upper == truth
    lowerEdgeInds = which(atLowerEdge)
    upperEdgeInds = which(atUpperEdge)
    
    # Fuzzy rejection probabilities at the interval boundaries. Nonzero only at 
    # respective interval boundaries
    rejectLower = rep(0, length(truth))
    rejectUpper = rep(0, length(truth))
    if(length(lowerEdgeInds) != 0) {
      probRejectLower = sapply(lowerEdgeInds, function(i) {
        if(mean(estMat[i,] == lower[i]) != 0) {
          ((1 - significance) / 2 - mean(estMat[i,] < lower[i], na.rm=na.rm)) / mean(estMat[i,] == lower[i], na.rm=na.rm)
        } else {
          warning("lower end of CI not equal to any samples from estMat. Setting reject probability to 0.5")
          0.5
        }
      })
      rejectLower[lowerEdgeInds] = probRejectLower
    }
    if(length(upperEdgeInds) != 0) {
      probRejectUpper = sapply(upperEdgeInds, function(i) {
        if(mean(estMat[i,] == upper[i]) != 0) {
          ((1 - significance) / 2 - mean(estMat[i,] > upper[i], na.rm=na.rm)) / mean(estMat[i,] == upper[i], na.rm=na.rm)
        } else {
          warning("upper end of CI not equal to any samples from estMat. Setting reject probability to 0.5")
          0.5
        }
      })
      rejectUpper[upperEdgeInds] = probRejectUpper
    }
    
    # determine minimum differences between probabilities
    if(is.null(ns))
      deltas = apply(estMat, 1, function(x) {min(c(0, diff(sort(unique(x)))), na.rm=na.rm)})
    else
      deltas = 1 / ns
    
    # reduce CI width based on fuzzy boundaries
    width = width - deltas*rejectLower - deltas*rejectUpper
    upper = upper - deltas*rejectUpper
    lower = lower + deltas*rejectLower
    # width = upper - lower (this should be the same as above)
    
    # cvg at the boundaries was 1, but the fuzzy coverage is 1 minus the rejection probability
    if(returnCoverage) {
      if(length(lowerEdgeInds) != 0) {
        # cvg[lowerEdgeInds] = sapply(lowerEdgeInds, function(i) {min(cvg[i], (1-rejectLower[i]))})
        cvg[lowerEdgeInds] = 1-rejectLower[lowerEdgeInds]
      }
      if(length(upperEdgeInds) != 0) {
        # cvg[upperEdgeInds] = sapply(upperEdgeInds, function(i) {min(cvg[i], (1-rejectUpper[i]))})
        cvg[upperEdgeInds] = 1-rejectUpper[upperEdgeInds]
      }
    }
  }
  
  # calculate interval score
  alpha = 1 - significance
  theseScores = upper - lower + 
    2/alpha * (lower - truth) * as.numeric(!greaterThanLower) + 
    2/alpha * (truth - upper) * as.numeric(!lessThanUpper)
  
  # aggregate individual scores if need be, and concatenate all relevant scores
  if(aggScores) {
    allResults = c(intScore=do.call("aggFun", c(list(theseScores), aggFunArgs)))
  } else {
    allResults = c(intScore=theseScores)
  }
  
  if(returnCoverage) {
    if(aggScores) {
      allResults = c(allResults, coverage=do.call("aggFun", c(list(cvg, significance=significance), aggFunArgs)))
    } else {
      allResults = cbind(allResults, coverage=cvg)
    }
  }
  
  if(returnIntervalWidth) {
    if(aggScores) {
      allResults = c(allResults, width=do.call("aggFun", c(list(width), aggFunArgs)))
    } else {
      allResults = cbind(allResults, width=width)
    }
  }
  
  # adjust names to match the significance in case there are multiple significances
  if(is.matrix(allResults)) {
    colnames(allResults) = paste(colnames(allResults), 100*significance, sep="")
  } else {
    names(allResults) = paste(names(allResults), 100*significance, sep="")
  }
  
  allResults
}

# averages a list of many tables, each returned from the getScores function with distanceBreaks set by user
averageBinnedScores = function(tableList) {
  
  if(length(tableList) == 1) {
    # base case
    
    return(as.data.frame(tableList[[1]]))
  } else if(is.null(tableList[[2]])) {
    # minor case (some of the tables might be NULL)
    
    return(averageBinnedScores(tableList[-2]))
  } else {
    # recursive case
    
    firstTable = tableList[[1]]
    secondTable = tableList[[2]]
    
    # make sure all tables are matrices
    firstTable = as.matrix(firstTable)
    secondTable = as.matrix(secondTable)
    
    # make sure tables have matched bins
    uniqueDists = sort(unique(c(firstTable[,1], secondTable[,1])))
    firstMatch = match(firstTable[,1], uniqueDists)
    secondMatch = match(secondTable[,1], uniqueDists)
    newFirstTable = matrix(0, nrow=length(uniqueDists), ncol=ncol(firstTable))
    newFirstTable[,1] = uniqueDists
    newFirstTable[firstMatch,] = firstTable
    colnames(newFirstTable) = colnames(firstTable)
    newSecondTable = matrix(0, nrow=length(uniqueDists), ncol=ncol(secondTable))
    newSecondTable[,1] = uniqueDists
    newSecondTable[secondMatch,] = secondTable
    colnames(newSecondTable) = colnames(secondTable)
    
    firstTable = newFirstTable
    secondTable = newSecondTable
    
    # calculate weights for averaging
    # ns1 = firstTable$nPerBin
    # ns2 = secondTable$nPerBin
    ns1 = firstTable[,2]
    ns2 = secondTable[,2]
    nsTotal = ns1 + ns2
    ws1 = ns1 / nsTotal
    ws2 = ns2 / nsTotal
    
    # perform weighted averaging
    newTable = firstTable
    newTable[,2] = ns1 + ns2
    newTable[,3:ncol(firstTable)] = sweep(firstTable[,3:ncol(firstTable)], 1, ws1, "*") + sweep(secondTable[,3:ncol(secondTable)], 1, ws2, "*")
    
    # return results recursively
    return(averageBinnedScores(c(list(newTable), tableList[-(1:2)])))
  }
}

aggregateScoresByDistance = function(singleScores, breaks=30, observationType=c("All", "Urban", "Rural"), predictionType=c("All", "Urban", "Rural"), 
                                     dat=NULL, targetPop=c("women", "children"), distanceBreaksType=c("quantiles", "even"), nPerBin=NULL, maxDist=Inf) {
  # NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA, dataI
  
  targetPop = match.arg(targetPop)
  observationType = match.arg(observationType)
  predictionType = match.arg(predictionType)
  distanceBreaksType = match.arg(distanceBreaksType)
  
  # subset prediction points by urbanicity if necessary. First determine whether prediction points are urban, then filter
  if(targetPop == "women") {
    resultNameRoot="Ed"
    if(is.null(dat)) {
      out = load("../U5MR/kenyaDataEd.RData")
      dat = ed
    }
    load("../U5MR/popGridAdjustedWomen.RData")
  } else if(targetPop == "children") {
    resultNameRoot="Mort"
    if(is.null(dat)) {
      out = load("../U5MR/kenyaData.RData")
      dat = mort
    }
    load("../U5MR/popGridAdjusted.RData")
  }
  predictionUrban = dat$urban[singleScores$dataI]
  
  if(predictionType == "Urban") {
    singleScores = singleScores[predictionUrban,]
  } else if(predictionType == "Rural") {
    singleScores = singleScores[!predictionUrban,]
  }
  
  # Now determine distance to which type of point we use
  distanceType = ""
  if(observationType=="Urban") {
    distanceType = "U"
  } else if(observationType == "Rural"){
    distanceType = "u"
  }
  distanceVar = paste0("NNDist", distanceType)
  distances = singleScores[[distanceVar]]
  
  # remove distances beyond the maximum of breaks
  if(length(breaks) != 1) {
    maxDist = min(maxDist, max(breaks))
  }
  badDistances = distances >= maxDist
  singleScores = singleScores[!badDistances,]
  distances = distances[!badDistances]
  
  # sort table by distances
  sortI = sort(distances, index.return=TRUE)$ix
  singleScores = singleScores[sortI,]
  distances = distances[sortI]
  
  # calculate default breaks for the bin limits if necessary
  if(length(breaks) == 1) {
    nBreaks = breaks
    if(distanceBreaksType == "even" && is.null(nPerBin)) {
      breaks = seq(0, max(distances), l=nBreaks)
    } else {
      if(is.null(nPerBin))
        nPerBin = ceiling(nrow(singleScores)/nBreaks)
      
      # get endpoints of the bins, average their values when calculating breaks
      startI = seq(nPerBin+1, nrow(singleScores), by=nPerBin)
      endI = startI - 1
      breaks = c(0, c(rowMeans(cbind(distances[startI], distances[endI]))), distances[length(distances)]+1e-6)
    }
  }
  
  # construct the distance bins with which to group the data and compute scores within
  binsI = cut(distances, breaks, labels=1:(length(breaks)-1), include.lowest=TRUE)
  centers = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  uniqueBinsI = sort(unique(binsI))
  
  # determine the number of observations per bin
  nPerBin = as.numeric(table(binsI))
  
  # helper function to compute the scoring rules for a given bin
  getSubScores = function(uniqueBinI) {
    thisDatI = binsI == uniqueBinI
    
    # thisSingleScoresBinomial = data.frame(c(list(Region=thisRegion, dataI=which(thisSampleI), NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA), getScores(truth, est, vars, lower, upper, estMatBinomial, getAverage=FALSE), Time=time[3]))
    colMeans(singleScores[thisDatI,-c(1, 2)])
  }
  
  # calculate scores for each bin individually
  binnedScores = t(sapply(uniqueBinsI, getSubScores))
  
  # make sure each variable in binnedScores is a numeric, not a list...
  temp = matrix(unlist(binnedScores), nrow=length(uniqueBinsI))
  theseNames = colnames(binnedScores)
  binnedScores = data.frame(temp)
  names(binnedScores) = theseNames
  out = as.data.frame(cbind(nPerBin=nPerBin[uniqueBinsI], binnedScores))
  out[[distanceVar]] = centers[uniqueBinsI]
  out
}

aggregateScoresByDistanceBasic = function(singleScores, breaks=30, distanceVar="NNDist", 
                                          distanceBreaksType=c("quantiles", "even"), 
                                          nPerBin=NULL, maxDist=Inf) {
  
  distanceBreaksType = match.arg(distanceBreaksType)
  
  # Now determine distance to which type of point we use
  distanceVarI = grep(distanceVar, colnames(singleScores))
  distances = singleScores[,distanceVarI]
  
  # remove distances beyond the maximum of breaks
  if(length(breaks) != 1) {
    maxDist = min(maxDist, max(breaks))
  }
  badDistances = distances >= maxDist
  singleScores = singleScores[!badDistances,]
  distances = distances[!badDistances]
  
  # sort table by distances
  sortI = sort(distances, index.return=TRUE)$ix
  singleScores = singleScores[sortI,]
  distances = distances[sortI]
  
  # calculate default breaks for the bin limits if necessary
  if(length(breaks) == 1) {
    nBreaks = breaks
    if(distanceBreaksType == "even" && is.null(nPerBin)) {
      breaks = seq(0, max(distances), l=nBreaks)
    } else {
      if(is.null(nPerBin))
        nPerBin = ceiling(nrow(singleScores)/nBreaks)
      
      # get endpoints of the bins, average their values when calculating breaks
      startI = seq(nPerBin+1, nrow(singleScores), by=nPerBin)
      endI = startI - 1
      breaks = c(0, c(rowMeans(cbind(distances[startI], distances[endI]))), distances[length(distances)]+1e-6)
    }
  }
  
  # construct the distance bins with which to group the data and compute scores within
  binsI = cut(distances, breaks, labels=1:(length(breaks)-1), include.lowest=TRUE)
  centers = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  uniqueBinsI = sort(unique(binsI))
  
  # determine the number of observations per bin
  nPerBin = as.numeric(table(binsI))
  
  # helper function to compute the scoring rules for a given bin
  getSubScores = function(uniqueBinI) {
    thisDatI = binsI == uniqueBinI
    
    colMeans(singleScores[thisDatI,-distanceVarI], na.rm=TRUE)
  }
  
  # calculate scores for each bin individually
  binnedScores = t(sapply(uniqueBinsI, getSubScores))
  
  # make sure each variable in binnedScores is a numeric, not a list...
  temp = matrix(unlist(binnedScores), nrow=length(uniqueBinsI))
  theseNames = colnames(binnedScores)
  binnedScores = data.frame(temp)
  names(binnedScores) = theseNames
  out = as.data.frame(cbind(nPerBin=nPerBin[uniqueBinsI], binnedScores))
  out[[distanceVar]] = centers[uniqueBinsI]
  out
}

# sort a vector containing NAs, which are placed at the end
sortWithNAs = function(x, ...) {
  numNAs = sum(is.na(x))
  out = sort(x, ...)
  c(out, rep(NA, numNAs))
}

# Score aggregation functions ----

# weighted mean that removes NAs
robustWeightedMean = function(x, weights=rep(1, length(x))) {
  nas = !is.finite(x)
  x = x[!nas]
  weights = weights[!nas]
  weights = weights * (1/sum(weights))
  sum(x*weights)
}

# get most nonzero value in the vector
getWorst = function(x, significance=NULL) {
  if(is.null(significance)) {
    res = abs(x)
  } else {
    res = x - significance
  }
  whichI = which.max(abs(res))
  sign(x[whichI]) * max(abs(x[whichI]))
}













