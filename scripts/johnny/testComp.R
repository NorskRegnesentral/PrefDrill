library(profvis)

# get seismic data
i = 1
out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedPred_", i, ".txt"), force01=TRUE)
seismicDat = out$surfFrame

# truth
out = readSurfaceRMS(paste0("data/seisTruthReplicates/RegularizedSand_", i, ".txt"), force01=TRUE)
truth = out$surfFrame

# generate prediction "distribution"
predMat = cbind(seismicDat[,3], seismicDat[,3])
meanSeis = mean(seismicDat[,3])
predAggMat = matrix(c(meanSeis, meanSeis), nrow=1)

# calculate scoring rules and metrics based on predictions
profvis(pwScoresMean <- getScores(truth[,3], estMat=predMat, significance=c(.8, .95)))

system.time(pwScoresMean <- getScores(truth[,3], estMat=predMat, significance=c(.8, .95)))[3]

system.time(
CIs <- apply(predMat, 1, function(ps) 
  {
    quantile(ps, probs=c((1 - .8) / 2, 1 - (1 - .8) / 2), na.rm=TRUE)
  }
  )
)[3]

system.time(
  {
    low = (1 - .8) / 2
    up = 1 - (1 - .8) / 2
    CIs = apply(predMat, 1, function(ps) 
    {
      quantile(ps, probs=c(low, up), na.rm=TRUE, names=FALSE)
    }
    )
  }
)[3] # 5.32

system.time(
  {
    low = (1 - .8) / 2
    up = 1 - (1 - .8) / 2
    CIs = apply(predMat, 1, function(ps) 
    {
      quantile(ps, probs=c(low, up), na.rm=FALSE, names=FALSE, type=7)
    }
    )
  }
)[3]
# 

install.packages("matrixStats")
library(matrixStats)

system.time(
  {
    low = (1 - .8) / 2
    up = 1 - (1 - .8) / 2
    CIs = row_q <- rowQuantiles(predMat, probs = c(low, up), useNames=TRUE)
  }
)[3] # 1.36

temp = predMat
temp[,1] = temp[,1] - .1

system.time({
  test = t(apply(temp, 1, sortWithNAs))
})[3] # 4.72 


system.time({
  test = rowRanks(temp)
  test2 = t(sapply(1:nrow(temp), function(i) {
    temp[i,][test[i,]]
  }))
  # test2 = temp[test]
})[3] # 0.17

finalPredMat = matrix(rep(seismicDat[,3], 1000), ncol=1000)
system.time(pwScoresMean <- getScores(truth[,3], estMat=finalPredMat, significance=c(.8, .95)))[3]
# with improvements: 39
# without improvements: 51.34

profvis(pwScoresMean <- getScores(truth[,3], estMat=finalPredMat, significance=c(.8, .95)))

system.time(pwScoresMean <- getScores(truth[,3], estMat=finalPredMat, significance=c(.8, .95)))[3]
# with improvements: 53
# without improvements: 41.97 
# with final improvements (including scoringRules::crps_sample): 25.96 

pwScoresMean <- getScores(truth[,3], estMat=finalPredMat, significance=c(.8, .95))

testTIs = sample(1:1000, nrow(finalPredMat), replace=TRUE)
remove_one_per_row
system.time(test0 <- remove_one_per_row(finalPredMat, testTIs)) # 3.98  
system.time(test <- remove_one_per_row_fast(finalPredMat, testTIs)) # 6.52 
system.time(test2 <- removeOnePerRow(finalPredMat, testTIs)) # 1.70 


finalPredMat = matrix(1:length(finalPredMat), byrow = TRUE, ncol=ncol(finalPredMat))
system.time(test0 <- remove_one_per_row(finalPredMat, testTIs)) # 3.19   
system.time(test <- remove_one_per_row_fast(finalPredMat, testTIs)) # 5.11  
system.time(test2 <- removeOnePerRow(finalPredMat, testTIs)) # 1.92 
head(testTIs)
finalPredMat[1:5,1:5]
test0[1:5,1:5]
test[1:5,1:5]
test2[1:5,1:5]
test2[2,1:64]

test = insertSortedColumnWithIndices(finalPredMat, rep(0, nrow(finalPredMat)))
test$newMat[1:5,1:5]
test$insertIndices[1:5]

test = insertSortedColumnWithIndices(finalPredMat, rep(length(finalPredMat)+1, nrow(finalPredMat)))
test$newMat[1:5,995:1001]
test$insertIndices[1:5]

system.time(getSeismicEsts(1, regenData=TRUE))[3]
# with final improvements: 22.81 

