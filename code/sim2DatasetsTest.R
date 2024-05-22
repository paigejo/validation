# Test cross validation methods in the case of two datasets with correlated 
# sampling rates. Fix a sampling distribution that is variable on unit square 
# based on GRF, and for a number of different grid resolutions. Error will vary 
# between LOO-CV error and asymptotic extrapolation error, the true predictive 
# error being somewhere in between
#   1.  for iter in 1:number of simulations:
#   2.    Simulate 1 GRF on unit square for responses
#   3.    Simulate 2 GRFs on unit square for sample distributions
#   4.    Simulate observations in domain based on the sampling distn and the GRF 
#         for responses include nugget
#   5.    Calculate full covariate matrix for sample and crossCov to all locs
#   6.    for i in 1:number of block resolutions:
#   7.      Group observations with blocks
#   8.      Leave out one cell at a time, calculate gridded CV MSE
#   9.    Leave out one observation at a time, calculate LOO-CV MSE
#   10.   Calculate importance weighted LOO-CV MSE, or LOOIS-CV MSE
#   11.   Calculate actual predictive MSE
#   12.   Save result
#   13. Combine and save all results
# Note that the below function is a single iteration of the for loop 
# for easy parallelization.
griddedResTestIter2Datasets = function(rGRFargsTruth=NULL, rGRFargsSample=NULL, 
                                       n1=50, n2=n1, gridNs=2^(1:6), iter=1, rGRFargsWrong=NULL, 
                                       nx=100, ny=100, sigmaEpsSq1=.1^2, sigmaEpsSq2=1^2, allSeeds=123, 
                                       alpha=0, beta=1, rho=-.8, printProgress=FALSE) {
  
  print(paste0("iteration ", iter, "/", length(allSeeds)))
  
  if(!is.null(allSeeds)) {
    set.seed(allSeeds[iter])
  }
  
  # set default GRF parameters
  if(is.null(rGRFargsTruth)) {
    rGRFargsTruth = list(mu=0, sigma=1, 
                         cov.args=list(Covariance="Matern", range=0.2, smoothness=0.5), 
                         delta=3, sigmaEpsSq=0, nx=nx, ny=ny)
  }
  if(is.null(rGRFargsWrong)) {
    rGRFargsWrong = rGRFargsTruth
    rGRFargsWrong$cov.args$range = 0.25 * rGRFargsTruth$cov.args$range
    rGRFargsWrong$cov.args$smoothness = 1.5
    if(n1 == 50) {
      # rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt(4)
      rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt(1.2)
    } else if(n1 == 500) {
      # rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt(2)
      rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt(1.1)
    } else if(n1 == 1000) {
      # rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt(1.75)
      rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt(1.075)
    } else {
      rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt((1.2-1)/sqrt(n1/50) + 1)
    }
  }
  if(is.null(rGRFargsSample)) {
    rGRFargsSample = list(mu1=rGRFargsTruth$mu, mu2=rGRFargsTruth$mu, 
                          sigma1=rGRFargsTruth$sigma, sigma2=rGRFargsTruth$sigma, rho=rho, 
                          cov.args=rGRFargsTruth$cov.args, delta=rGRFargsTruth$delta, 
                          sigmaEpsSq1=0, 
                          sigmaEpsSq2=0, 
                          nx=nx, ny=ny)
  }
  
  # 2. Simulate 1 GRF for responses ----
  #    on unit square
  truthGRF = do.call("rGRF", rGRFargsTruth)
  truth = truthGRF$truth
  locs = truthGRF$locs
  
  # 3. Simulate sample distribution ----
  #    as 2 GRF on unit square
  sampleGRFs = do.call("r2GRFs", rGRFargsSample)
  sampleGRF1 = sampleGRFs$Y1
  sampleGRF2 = sampleGRFs$Y2
  sampleRates1 = exp(sampleGRF1$truth)
  sampleProbs1 = sampleRates1 * (1/sum(sampleRates1))
  sampleRates2 = exp(sampleGRF2$truth)
  sampleProbs2 = sampleRates2 * (1/sum(sampleRates2))
  
  # relevant sampleRates and sampleProbs for validation are those of the first dataset
  sampleRates = sampleRates1
  sampleProbs = sampleProbs1
  
  # 4. Simulate observations ----
  #    in domain based on the sampling distn and the GRF 
  #    for responses (don't include nugget?)
  gridResX = 1/nx
  gridResY = 1/ny
  sampleI1 = sample(1:nrow(locs), n1, prob=sampleProbs1, replace=TRUE)
  xs1 = locs[sampleI1,] + cbind(runif(n1, max=gridResX)-gridResX/2, 
                                runif(n1, max=gridResY)-gridResY/2)
  ys1 = truth[sampleI1] + rnorm(n1, sd=sqrt(sigmaEpsSq1))
  sampleI2 = sample(1:nrow(locs), n2, prob=sampleProbs2, replace=TRUE)
  xs2 = locs[sampleI2,] + cbind(runif(n2, max=gridResX)-gridResX/2, 
                                runif(n2, max=gridResY)-gridResY/2)
  ys2 = truth[sampleI2] + rnorm(n2, sd=sqrt(sigmaEpsSq2))
  
  xs = rbind(xs1, xs2)
  ys = c(ys1, ys2)
  
  # 5. Get covariance matrices ----
  #    for sample and crossCov to all locs
  distMatXs = rdist(xs, xs)
  distXs1ToXs1 = distMatXs[1:n1, 1:n1]
  distXs1ToXs2 = distMatXs[1:n1, (n1+1):(n1+n2)]
  SigmaSample = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsTruth$cov.args$range, 
                               smoothness=rGRFargsTruth$cov.args$smoothness, distMat=distMatXs) * rGRFargsTruth$sigma^2
  SigmaSampleWrong = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsWrong$cov.args$range, 
                                    smoothness=rGRFargsWrong$cov.args$smoothness, distMat=distMatXs) * rGRFargsWrong$sigma^2
  SigmaSample = SigmaSample + diag(c(rep(sigmaEpsSq1, n1), rep(sigmaEpsSq2, n2)))
  SigmaSampleWrong = SigmaSampleWrong + diag(c(rep(sigmaEpsSq1, n1), rep(sigmaEpsSq2, n2)))
  distLocsToXs1 = rdist(locs, xs1)
  distLocsToXs2 = rdist(locs, xs2)
  distLocsToSample = cbind(distLocsToXs1, distLocsToXs2)
  SigmaGridToSample = stationary.cov(locs, xs, Covariance="Matern", aRange=rGRFargsTruth$cov.args$range, 
                                     smoothness=rGRFargsTruth$cov.args$smoothness, distMat=distLocsToSample) * rGRFargsTruth$sigma^2
  SigmaGridToSampleWrong = stationary.cov(locs, xs, Covariance="Matern", aRange=rGRFargsWrong$cov.args$range, 
                                          smoothness=rGRFargsWrong$cov.args$smoothness, distMat=distLocsToSample) * rGRFargsWrong$sigma^2
  
  # 6. for i in 1:number block resolutions: ----
  
  # determine grid cells associated with observations
  getCellI = function(thisLoc) {
    xI = match(FALSE, thisLoc[1] > highs)
    yI = match(FALSE, thisLoc[2] > highs)
    
    indMat[yI, xI]
  }
  
  griddedCVs = numeric(length(gridNs))
  griddedCVsWrong = numeric(length(gridNs))
  griddedRCVs = numeric(length(gridNs))
  griddedRCVsWrong = numeric(length(gridNs))
  griddedR2CVs = numeric(length(gridNs))
  griddedR2CVsWrong = numeric(length(gridNs))
  for(i in 1:length(gridNs)) {
    gridN = gridNs[i]
    if(printProgress) {
      print(paste0("grid resolution ", gridN, " (", i, "/", length(gridNs), ")"))
    }
    
    # 7. Group data by block ----
    # index cells and select test/train cells:
    # cell i is in:
    # row floor(i) + 1 (row counts upward along unit square)
    # col floor(i/row) + i
    nCellsTot = gridN^2
    indMat = matrix(1:nCellsTot, nrow=gridN)
    
    highs = seq(0, 1, l=gridN+1)[-1]
    cellIsSample = apply(xs1, 1, getCellI) # only leave out observations from first dataset
    uniqueCellIs = sort(unique(cellIsSample))
    blockCVs = numeric(length(uniqueCellIs))
    blockCVsWrong = numeric(length(uniqueCellIs))
    n1s = numeric(length(uniqueCellIs))
    n2s = numeric(length(uniqueCellIs))
    
    cellIsSample2 = apply(xs2, 1, getCellI) # Also get observations from second dataset in each block
    
    # 8. Get gridded CV MSE ----
    for(j in 1:length(uniqueCellIs)) {
      if(printProgress) {
        print(paste0("Leaving out cell ", j, "/", length(uniqueCellIs)))
      }
      
      cellI = uniqueCellIs[j]
      
      # separate data in test/train
      isTest = c(cellIsSample == cellI, rep(FALSE, n2))
      testYs = ys[isTest]
      trainYs = ys[!isTest]
      
      # check number of observations in block from each dataset
      n1s[j] = sum(isTest)
      n2s[j] = sum(cellIsSample2 == cellI)
      
      # predict test data for both true and wrong model
      SigmaAB = matrix(SigmaSample[isTest, !isTest], nrow=n1s[j])
      SigmaBB = SigmaSample[!isTest, !isTest]
      SigmaAA = SigmaSample[isTest, isTest]
      
      SigmaABWrong = matrix(SigmaSampleWrong[isTest, !isTest], nrow=n1s[j])
      SigmaBBWrong = SigmaSampleWrong[!isTest, !isTest]
      SigmaAAWrong = SigmaSampleWrong[isTest, isTest]
      
      if((n1s[j] > 0) && (sum(!isTest) > 0)) {
        condDistn = condMeanMVN(SigmaAA=SigmaAA, SigmaAB=SigmaAB, SigmaBB=SigmaBB, 
                                ysB=trainYs, getFullCov=FALSE, getCondVar=FALSE)
        muAcondB = condDistn$muAcondB
        
        condDistnWrong = condMeanMVN(SigmaAA=SigmaAAWrong, SigmaAB=SigmaABWrong, SigmaBB=SigmaBBWrong, 
                                     ysB=trainYs, getFullCov=FALSE, getCondVar=FALSE)
        muAcondBWrong = condDistnWrong$muAcondB
        
        # calculate MSE for the grid cell
        blockCVs[j] = mean((testYs - muAcondB)^2)
        blockCVsWrong[j] = mean((testYs - muAcondBWrong)^2)
      } else {
        blockCVs[j] = NA
        blockCVsWrong[j] = NA
      }
    }
    
    # average over MSEs of each grid cell
    griddedCVs[i] = mean(blockCVs, na.rm=TRUE)
    griddedCVsWrong[i] = mean(blockCVsWrong, na.rm=TRUE)
    
    # when calculating controlled/robust gridded CV, notes that each block is an 
    # areal error estimate, not point level. Hence, IWs are equal to 1, since 
    # grid cell areas are all equal.
    thisGridCellArea = 1/gridNs[i]^2
    rate1Ests = (n1s/thisGridCellArea)/n1
    rate2Ests = (n2s/thisGridCellArea)/n2
    griddedRCVs[i] = getRobustIWCV(scores=blockCVs, IWs=1, controlVarMat=1/rate1Ests, controlVarMeans=1)
    griddedRCVsWrong[i] = getRobustIWCV(scores=blockCVsWrong, IWs=1, controlVarMat=1/rate1Ests, controlVarMeans=1)
    griddedR2CVs[i] = getRobustIWCV(scores=blockCVs, IWs=1, 
                                    controlVarMat=cbind(1/rate1Ests, rate2Ests/rate1Ests), controlVarMeans=c(1, 1))
    griddedR2CVsWrong[i] = getRobustIWCV(scores=blockCVsWrong, IWs=1, 
                                         controlVarMat=cbind(1/rate1Ests, rate2Ests/rate1Ests), controlVarMeans=c(1, 1))
  }
  
  # 9. get LOO-CV MSE ----
  LOOCVs = numeric(n1)
  LOOCVsWrong = numeric(n1)
  for(i in 1:n1) {
    if(printProgress) {
      print(paste0("Leaving out obs ", i, "/", n1))
    }
    
    # separate data in test/train
    isTest = (1:(n1+n2)) == i
    testYs = ys[isTest]
    trainYs = ys[!isTest]
    
    # predict test data
    SigmaAB = matrix(SigmaSample[isTest, !isTest], nrow=1)
    SigmaBB = SigmaSample[!isTest, !isTest]
    SigmaAA = SigmaSample[isTest, isTest]
    
    SigmaABWrong = matrix(SigmaSampleWrong[isTest, !isTest], nrow=1)
    SigmaBBWrong = SigmaSampleWrong[!isTest, !isTest]
    SigmaAAWrong = SigmaSampleWrong[isTest, isTest]
    
    condDistn = condMeanMVN(SigmaAA=SigmaAA, SigmaAB=SigmaAB, SigmaBB=SigmaBB, 
                            ysB=trainYs, getFullCov=FALSE, getCondVar=FALSE)
    muAcondB = condDistn$muAcondB
    
    condDistnWrong = condMeanMVN(SigmaAA=SigmaAAWrong, SigmaAB=SigmaABWrong, SigmaBB=SigmaBBWrong, 
                                 ysB=trainYs, getFullCov=FALSE, getCondVar=FALSE)
    muAcondBWrong = condDistnWrong$muAcondB
    
    # calculate MSE for the grid cell
    LOOCVs[i] = mean((testYs - muAcondB)^2)
    LOOCVsWrong[i] = mean((testYs - muAcondBWrong)^2)
  }
  LOOCV = mean(LOOCVs)
  LOOCVWrong = mean(LOOCVsWrong)
  
  # get control variates
  minDists1Grid = apply(distLocsToXs1, 1, min)
  minDist1Mu = mean(minDists1Grid)
  minDists2Grid = apply(distLocsToXs2, 1, min)
  minDist2Mu = mean(minDists2Grid)
  minDists1 = apply(distXs1ToXs1, 1, function(x) {min(x[x != 0])})
  minDists2 = apply(distXs1ToXs2, 1, min)
  
  # LOORCV = getRobustIWCV(scores=LOOCVs, IWs=rep(1, n1), controlVarMat=minDists1, 
  #                        controlVarMeans=minDist1Mu, wtdReg=FALSE)
  LOORCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=minDists1, muf=minDist1Mu)[1]
  # LOORCVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=rep(1, n1), controlVarMat=minDists1, 
  #                             controlVarMeans=minDist1Mu, wtdReg=FALSE)
  LOORCVWrong = getWeightedControlVarEst(Ys=LOOCVsWrong, Zf=minDists1, muf=minDist1Mu)[1]
  
  # LOOR2CV = getRobustIWCV(scores=LOOCVs, IWs=rep(1, n1), controlVarMat=cbind(minDists1, minDists2), 
  #                        controlVarMeans=c(minDist1Mu, minDist2Mu), wtdReg=FALSE)
  LOOR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=cbind(minDists1, minDists2), muf=c(minDist1Mu, minDist2Mu))[1]
  # LOOR2CVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=rep(1, n1), controlVarMat=cbind(minDists1, minDists2), 
  #                              controlVarMeans=c(minDist1Mu, minDist2Mu), wtdReg=FALSE)
  LOOR2CVWrong = getWeightedControlVarEst(Ys=LOOCVsWrong, Zf=cbind(minDists1, minDists2), muf=c(minDist1Mu, minDist2Mu))[1]
  
  # 10. Calculate LOOIS-CV MSE ----
  cellArea = 1/(nx*ny)
  LOOISrates = (sampleProbs[sampleI1]/sum(sampleProbs)) / cellArea # convert from areal rate to rate density
  # LOOISCV = weighted.mean(LOOCVs, w=cellArea/LOOISrates)
  LOOISCV = sum(LOOCVs * (1/LOOISrates))/n1
  LOOISCVWrong = sum(LOOCVsWrong * (1/LOOISrates))/n1
  
  # estimate ratio of MSEs
  out = getP(LOOCVs, 1/LOOISrates)
  p = out$p
  LOOISPCV = sum(LOOCVs * (1/LOOISrates^p))/sum(1/LOOISrates^p)
  LOOISPCVWrong = sum(LOOCVsWrong * (1/LOOISrates^p))/sum(1/LOOISrates^p)
  
  # LOOISRCV = getRobustIWCV(scores=LOOCVs, IWs=1/LOOISrates, printVarRatio=FALSE, printRoot="LOOISR", controlVarMeans=1)
  # LOOISRCV = getRobustIWCV(scores=LOOCVs, IWs=1/LOOISrates, controlVarMat=LOOISrates, printVarRatio=FALSE, printRoot="LOOISR", controlVarMeans=1, wtdReg=TRUE)
  # LOOISRCV = getRobustIWCV(scores=LOOCVs, IWs=1/LOOISrates, controlVarMat=minDists1, printVarRatio=FALSE, printRoot="LOOISR", controlVarMeans=minDist1Mu, wtdReg=FALSE)
  # LOOISRCV = getGREGCV(scores=LOOCVs, IWs=1/LOOISrates, controlVarMat=cbind(1, minDists1), controlVarMeans=c(1, minDist1Mu), normalizeWeights=FALSE, 
  #                      shrinkWeights=FALSE)
  
  LOOISRCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=minDists1, muf=minDist1Mu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates)[1]
  # pGregInfo = getPGREG(LOOCVs, IWs=1/LOOISrates, controlVarMat=minDists1, controlVarMeans=minDist1Mu)
  
  # LOOISRCVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMeans=1)
  # LOOISRCVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMat=LOOISrates, controlVarMeans=1, wtdReg=TRUE)
  # LOOISRCVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMat=minDists1, controlVarMeans=minDist1Mu, wtdReg=FALSE)
  # LOOISRCVWrong = getGREGCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMat=cbind(1, minDists1), controlVarMeans=c(1, minDist1Mu), normalizeWeights=FALSE, 
  #                           shrinkWeights=FALSE)
  LOOISRCVWrong = getWeightedControlVarEst(Ys=LOOCVsWrong, Zf=minDists1, muf=minDist1Mu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates)[1]
  # we know that both 1/rateEsts and rateEsts2/rateEsts have expected value (over 
  # the sampling distn of dataset 1) equal to 1. They can therefore both improve 
  # the estimation. Get the sampling rate of dataset 2 at each location in dataset 1
  LOOISrates2 = (sampleProbs2[sampleI1]/sum(sampleProbs2)) / cellArea
  # LOOISR2CV = getRobustIWCV(scores=LOOCVs, IWs=1/LOOISrates, controlVarMat=cbind(1/LOOISrates, LOOISrates2/LOOISrates), controlVarMeans=c(1, NA), 
  #                           printVarRatio=FALSE, printRoot="LOOISR2")
  # LOOISR2CVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMat=cbind(1/LOOISrates, LOOISrates2/LOOISrates), controlVarMeans=c(1, NA))
  # LOOISR2CV = getRobustIWCV(scores=LOOCVs, IWs=1/LOOISrates, controlVarMat=cbind(1/LOOISrates, LOOISrates2/LOOISrates), controlVarMeans=c(1, 1), 
  #                           printVarRatio=FALSE, printRoot="LOOISR2")
  # LOOISR2CV = getRobustIWCV(scores=LOOCVs, IWs=1/LOOISrates, controlVarMat=cbind(LOOISrates, LOOISrates2), controlVarMeans=c(1, 1), 
  #                           printVarRatio=FALSE, printRoot="LOOISR2", wtdReg=TRUE)
  # LOOISR2CV = getRobustIWCV(scores=LOOCVs, IWs=1/LOOISrates, controlVarMat=cbind(minDists1, minDists2), controlVarMeans=c(minDist1Mu, minDist2Mu), 
  #                           printVarRatio=FALSE, printRoot="LOOISR2", wtdReg=FALSE)
  # LOOISR2CV = getGREGCV(scores=LOOCVs, IWs=1/LOOISrates, controlVarMat=cbind(1, minDists1, minDists2), controlVarMeans=c(1, minDist1Mu, minDist2Mu), normalizeWeights=FALSE, 
  #                       shrinkWeights = FALSE)
  # LOOISR2CV = getGREGCV(scores=LOOCVs, IWs=1/LOOISrates, controlVarMat=cbind(rep(1, length(minDists1))), controlVarMeans=1, normalizeWeights=FALSE, 
  #                       shrinkWeights = FALSE)
  LOOISR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=minDists1, muf=minDist1Mu, Zg=cbind(1/LOOISrates, LOOISrates2/LOOISrates), mug=c(1, 1), ws=1/LOOISrates)[1]
  # LOOISR2CVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMat=cbind(1/LOOISrates, LOOISrates2/LOOISrates), controlVarMeans=c(1))
  # LOOISR2CVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMat=cbind(LOOISrates, LOOISrates2), controlVarMeans=c(1, 1), wtdReg=TRUE)
  # LOOISR2CVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMat=cbind(minDists1, minDists2), controlVarMeans=c(minDist1Mu, minDist2Mu), wtdReg=FALSE)
  # LOOISR2CVWrong = getGREGCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMat=cbind(1, minDists1, minDists2), controlVarMeans=c(1, minDist1Mu, minDist2Mu), normalizeWeights=FALSE, 
  #                            shrinkWeights = FALSE)
  # LOOISR2CVWrong = getGREGCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMat=cbind(rep(1, length(minDists1))), controlVarMeans=1, normalizeWeights=FALSE, 
  #                            shrinkWeights = FALSE)
  LOOISR2CVWrong = getWeightedControlVarEst(Ys=LOOCVsWrong, Zf=minDists1, muf=minDist1Mu, Zg=cbind(1/LOOISrates, LOOISrates2/LOOISrates), mug=c(1, 1), ws=1/LOOISrates)[1]
  
  # LOOISPRCV = getGREGCV(scores=LOOCVs, IWs=1/LOOISrates^p, controlVarMat=cbind(1, minDists1), controlVarMeans=c(1, minDist1Mu), normalizeWeights=TRUE)
  LOOISPRCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=minDists1, muf=minDist1Mu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, shrinkWeights=TRUE)[1]
  # LOOISPRCVWrong = getGREGCV(scores=LOOCVsWrong, IWs=1/LOOISrates^p, controlVarMat=cbind(1, minDists1), controlVarMeans=c(1, minDist1Mu), normalizeWeights=TRUE)
  LOOISPRCVWrong = getWeightedControlVarEst(Ys=LOOCVsWrong, Zf=minDists1, muf=minDist1Mu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, shrinkWeights=TRUE)[1]
  # LOOISPR2CV = getGREGCV(scores=LOOCVs, IWs=1/LOOISrates^p, 
  #                        controlVarMat=cbind(1, minDists1, minDists2), controlVarMeans=c(1, minDist1Mu, minDist2Mu), normalizeWeights=TRUE)
  LOOISPR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=c(minDists1), muf=c(minDist1Mu), Zg=c(1/LOOISrates, LOOISrates2/LOOISrates), mug=c(1, 1), ws=1/LOOISrates, shrinkWeights=TRUE)[1]
  # LOOISPR2CVWrong = getGREGCV(scores=LOOCVsWrong, IWs=1/LOOISrates^p, 
  #                        controlVarMat=cbind(1, minDists1, minDists2), controlVarMeans=c(1, minDist1Mu, minDist2Mu), normalizeWeights=TRUE)
  LOOISPR2CVWrong = getWeightedControlVarEst(Ys=LOOCVsWrong, Zf=c(minDists1), muf=c(minDist1Mu), Zg=c(1/LOOISrates, LOOISrates2/LOOISrates), mug=c(1, 1), ws=1/LOOISrates, shrinkWeights=TRUE)[1]
  
  meanIW1 = mean(1/LOOISrates)
  meanIW2 = mean(LOOISrates2/LOOISrates)
  print(paste0("Empirical mean IW1: ", meanIW1, " IW2: ", meanIW2))
  cors1IS = cor(1/LOOISrates, LOOCVs/LOOISrates)
  cors2IS = cor(LOOISrates2/LOOISrates, LOOCVs/LOOISrates)
  
  # 11. Calculate LOOVC-CV MSE ----
  domainPoly = rbind(c(0, 0), 
                     c(0, 1), 
                     c(1, 1), 
                     c(1, 0), 
                     c(0, 0))
  vcellInfo = getVCellAreas(xs1, domainPoly=domainPoly)
  vcellArea = vcellInfo$area
  rateEsts = 1/vcellArea
  rateEsts = rateEsts/n1 # divide by n to get rate for a single pt draw instead of n pts
  
  LOOVCCV = sum(LOOCVs * (1/rateEsts))/n1
  LOOVCCVWrong = sum(LOOCVsWrong * (1/rateEsts))/n1
  # LOOVCRCV = getRobustIWCV(scores=LOOCVs, IWs=1/rateEsts, printVarRatio=FALSE, printRoot="LOOVCR", controlVarMeans=1)
  # LOOVCRCV = getRobustIWCV(scores=LOOCVs, IWs=1/rateEsts, controlVarMat=rateEsts, 
  #                          printVarRatio=FALSE, printRoot="LOOVCR", controlVarMeans=1, wtdReg=TRUE)
  # browser()
  # LOOVCRCV = getRobustIWCV(scores=LOOCVs, IWs=1/rateEsts, controlVarMat=minDists1, controlVarMeans=minDist1Mu, 
  #                          printVarRatio=FALSE, printRoot="LOOVCR", wtdReg=FALSE)
  LOOVCRCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=c(minDists1), muf=c(minDist1Mu), ws=1/rateEsts, shrinkWeights=FALSE)[1]
  # LOOVCRCVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/rateEsts, controlVarMat=rateEsts, 
  #                          printVarRatio=FALSE, printRoot="LOOVCR", controlVarMeans=1, wtdReg=TRUE)
  # LOOVCRCVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/rateEsts, controlVarMat=minDists1, controlVarMeans=minDist1Mu, 
  #                               printVarRatio=FALSE, printRoot="LOOVCR", wtdReg=FALSE)
  LOOVCRCVWrong = getWeightedControlVarEst(Ys=LOOCVsWrong, Zf=c(minDists1), muf=c(minDist1Mu), ws=1/rateEsts, shrinkWeights=FALSE)[1]
  
  # try using density estimate of second dataset as a control variate
  vcellInfo2 = getVCellAreas(xs2, domainPoly=domainPoly)
  vcellArea2 = vcellInfo2$area # area of each centroidal Voronoi cell
  rateEsts2 = 1/vcellArea2
  rateEsts2 = rateEsts2/n2 # divide by n to get rate for a single pt draw instead of n pts
  
  # figure out density estimates of second dataset at first dataset locations
  dists = rdist(xs1, vcellInfo2$centroids)
  minIs = apply(dists, 1, which.min)
  rateEsts2 = rateEsts2[minIs]
  
  # LOOVCR2CV = getRobustIWCV(scores=LOOCVs, IWs=1/rateEsts, controlVarMat=cbind(rateEsts, rateEsts2), 
  #                          printVarRatio=FALSE, printRoot="LOOVCR", controlVarMeans=1, wtdReg=TRUE)
  # LOOVCR2CV = getRobustIWCV(scores=LOOCVs, IWs=1/rateEsts, controlVarMat=cbind(minDists1, minDists2), 
  #                           printVarRatio=FALSE, printRoot="LOOVCR", controlVarMeans=c(minDist1Mu, minDist2Mu), wtdReg=FALSE)
  LOOVCR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=c(minDists1), muf=c(minDist1Mu), Zg=c(rateEsts2/rateEsts), mug=1, ws=1/rateEsts, shrinkWeights=FALSE)[1]
  # LOOVCR2CVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/rateEsts, controlVarMat=cbind(rateEsts, rateEsts2), 
  #                               printVarRatio=FALSE, printRoot="LOOVCR", controlVarMeans=1, wtdReg=TRUE)
  # LOOVCR2CVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/rateEsts, controlVarMat=cbind(minDists1, minDists2), 
  #                                printVarRatio=FALSE, printRoot="LOOVCR", controlVarMeans=c(minDist1Mu, minDist2Mu), wtdReg=FALSE)
  LOOVCR2CVWrong = getWeightedControlVarEst(Ys=LOOCVsWrong, Zf=c(minDists1), muf=c(minDist1Mu), Zg=c(rateEsts2/rateEsts), mug=1, ws=1/rateEsts, shrinkWeights=FALSE)[1]
  
  cvcellInfo = getVCellAreas(xs1, domainPoly=domainPoly, useCVC = TRUE)
  vcellArea = cvcellInfo$ptArea
  ptsPerArea = cvcellInfo$ptNPerArea
  rateEsts = ptsPerArea/vcellArea
  rateEsts = rateEsts/n1 # divide by n to get rate for a single pt draw instead of n pts
  
  LOOCVCCV = sum(LOOCVs * (1/rateEsts))/n1
  # LOOCVCRCV = getRobustIWCV(scores=LOOCVs, IWs=1/rateEsts, controlVarMat=rateEsts, 
  #                           printVarRatio=FALSE, printRoot="LOOVCR", controlVarMeans=1)
  # LOOCVCRCV = getRobustIWCV(scores=LOOCVs, IWs=1/rateEsts, controlVarMat=rateEsts, 
  #                           printVarRatio=FALSE, printRoot="LOOVCR", controlVarMeans=1, wtdReg=TRUE)
  # LOOCVCRCV = getRobustIWCV(scores=LOOCVs, IWs=1/rateEsts, controlVarMat=minDists1, controlVarMeans=minDist1Mu, 
  #                           printVarRatio=FALSE, printRoot="LOOVCR", wtdReg=FALSE, wtdMeanCentered=TRUE)
  LOOCVCRCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=c(minDists1), muf=c(minDist1Mu), ws=1/rateEsts, shrinkWeights=FALSE)[1]
  LOOCVCCVWrong = sum(LOOCVsWrong * (1/rateEsts))/n1
  # LOOCVCRCVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/rateEsts, 
  #                                controlVarMat=rateEsts, controlVarMeans=1, wtdReg=TRUE)
  # LOOCVCRCVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/rateEsts, 
  #                                controlVarMat=minDists1, controlVarMeans=minDist1Mu, wtdReg=FALSE, wtdMeanCentered=TRUE)
  LOOCVCRCVWrong = getWeightedControlVarEst(Ys=LOOCVsWrong, Zf=c(minDists1), muf=c(minDist1Mu), ws=1/rateEsts, shrinkWeights=FALSE)[1]
  
  # try using density estimate of second dataset as a control variate
  cvcellInfo2 = getVCellAreas(xs2, domainPoly=domainPoly, useCVC = TRUE)
  vcellArea2 = cvcellInfo2$area # area of each centroidal Voronoi cell
  ptsPerArea2 = cvcellInfo2$ptNPerArea
  rateEsts2 = ptsPerArea2/vcellArea2
  rateEsts2 = rateEsts2/n2 # divide by n to get rate for a single pt draw instead of n pts
  
  # figure out density estimates of second dataset at first dataset locations
  dists = rdist(xs1, cvcellInfo2$centroids)
  minIs = apply(dists, 1, which.min)
  rateEsts2 = rateEsts2[minIs]
  
  # we know that both 1/rateEsts and rateEsts2/rateEsts have expected value (over 
  # the sampling distn of dataset 1) equal to 1. They can therefore both improve 
  # the estimation.
  cors1CVC = cor(1/rateEsts, LOOCVs/rateEsts)
  cors2CVC = cor(rateEsts2/rateEsts, LOOCVs/rateEsts)
  # LOOCVCR2CV = getRobustIWCV(scores=LOOCVs, IWs=1/rateEsts, controlVarMat=cbind(1/rateEsts, rateEsts2/rateEsts), printVarRatio=FALSE, printRoot="LOOCVCR2", controlVarMeans=c(1, 1))
  # LOOCVCR2CVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/rateEsts, controlVarMat=cbind(1/rateEsts, rateEsts2/rateEsts), controlVarMeans=c(1, 1))
  # LOOCVCR2CV = getRobustIWCV(scores=LOOCVs, IWs=1/rateEsts, controlVarMat=cbind(rateEsts, rateEsts2), 
  #                            printVarRatio=FALSE, printRoot="LOOCVCR2", controlVarMeans=c(1, 1), wtdReg=TRUE)
  # LOOCVCR2CV = getRobustIWCV(scores=LOOCVs, IWs=1/rateEsts, 
  #                            controlVarMat=cbind(minDists1, minDists2), controlVarMeans=c(minDist1Mu, minDist2Mu), 
  #                            printVarRatio=FALSE, printRoot="LOOCVCR2", wtdReg=FALSE)
  LOOCVCR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=c(minDists1), muf=c(minDist1Mu), Zg=rateEsts2/rateEsts, mug=1, ws=1/rateEsts, shrinkWeights=FALSE)[1]
  # LOOCVCR2CVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/rateEsts, controlVarMat=cbind(rateEsts, rateEsts2), 
  #                                 controlVarMeans=c(1, 1), wtdReg=TRUE)
  # LOOCVCR2CVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/rateEsts, 
  #                                 controlVarMat=cbind(minDists1, minDists2), controlVarMeans=c(minDist1Mu, minDist2Mu), 
  #                                 wtdReg=FALSE)
  LOOCVCR2CVWrong = getWeightedControlVarEst(Ys=LOOCVsWrong, Zf=c(minDists1), muf=c(minDist1Mu), Zg=rateEsts2/rateEsts, mug=1, ws=1/rateEsts, shrinkWeights=FALSE)[1]
  
  meanCVCIW1 = mean(1/rateEsts)
  meanCVCIW2 = mean(rateEsts2/rateEsts)
  print(paste0("Empirical mean CVCIW1: ", meanCVCIW1, " CVCIW2: ", meanCVCIW2))
  
  # 12. Calculate model MSEs ----
  condDistn = condMeanMVN(SigmaAA=rep(rGRFargsTruth$sigma^2, nrow(locs)), SigmaAB=SigmaGridToSample, SigmaBB=SigmaSample, 
                          ysB=ys, getFullCov=FALSE, getCondVar=FALSE)
  muAcondB = condDistn$muAcondB
  condDistn = condMeanMVN(SigmaAA=rep(rGRFargsWrong$sigma^2, nrow(locs)), SigmaAB=SigmaGridToSampleWrong, SigmaBB=SigmaSampleWrong, 
                          ysB=ys, getFullCov=FALSE, getCondVar=FALSE)
  muAcondBwrong = condDistn$muAcondB
  
  # calculate MSE for the grid cell
  trueMSE = mean((truth - muAcondB)^2) + sigmaEpsSq1
  wrongMSE = mean((truth - muAcondBwrong)^2) + sigmaEpsSq1 # add sigmaEpsSq1, the true error variance, not sigmaEpsSq2
  
  # rGRFargsTruth=NULL, rGRFargsSample=NULL, 
  # n=50, gridNs=2^(1:6), iter=1, seed=123, 
  # nx=100, ny=100, sigmaEpsSq=0
  
  if(FALSE) {
    # out = estPredMSErange(sigma=1, cov.args=rGRFargsTruth$cov.args, doPlot=TRUE, 
    #                       nPts=n)
    
    
    plot(gridNs, griddedCVs, type="n", log="x", axes=FALSE, 
         ylim=c(0, max(griddedCVs)), xlab="Blocks per side", 
         ylab="MSE", main="Gridded CV vs resolution")
    axis(side=1, at=gridNs)
    axis(side=2)
    box()
    lines(gridNs, griddedCVs, type="o", pch=19, col="blue")
    abline(h=LOOCV, lty=2, col="purple")
    abline(h=trueMSE, lty=2, col="green")
    legend("topright", c("Gridded", "LOO", "Truth"), col=c("blue", "purple", "green"), 
           pch=c(19, NA, NA), lty=c(1, 2, 2))
    
    plot(gridNs, griddedCVs, type="n", log="x", axes=FALSE, 
         ylim=c(0, 1+sigmaEpsSq), xlab="Blocks per side", 
         ylab="MSE", main="Gridded CV vs resolution")
    axis(side=1, at=gridNs)
    axis(side=2)
    box()
    lines(gridNs, griddedCVs, type="o", pch=19, col="blue")
    abline(h=LOOCV, lty=2, col="purple")
    abline(h=trueMSE, lty=2, col="green")
    abline(h=1+sigmaEpsSq, lty=2, col="red")
    legend("right", c("Gridded", "LOO", "Extrapolation", "Interpolation"), 
           col=c("blue", "purple", "red", "green"), 
           pch=c(19, NA, NA, NA), lty=c(1, 2, 2, 2))
  }
  
  # calculate some theoretical properties of the sampling weights
  # mean(1/sampleGRF1$truth)
  # mean(1/sampleGRF2$truth)
  trueVarW1 = -1 + mean(1/sampleProbs1)
  trueVarW2 = -1 + mean(1/sampleProbs2)
  
  # 13. Save result ----
  save(trueMSE, wrongMSE, LOOCVs, LOOCVsWrong, 
       LOOCV, LOORCV, LOOR2CV, 
       LOOCVWrong, LOORCVWrong, LOOR2CVWrong, 
       LOOISCV, LOOISRCV, LOOISR2CV, 
       LOOISPCV, LOOISPRCV, LOOISPR2CV, 
       LOOISPCVWrong, LOOISPRCVWrong, LOOISPR2CVWrong, 
       LOOISCVWrong, LOOISRCVWrong, LOOISR2CVWrong, 
       LOOVCCV, LOOVCRCV, LOOVCR2CV, 
       LOOCVCCV, LOOCVCRCV, LOOCVCR2CV, 
       LOOVCCVWrong, LOOVCRCVWrong, LOOVCR2CVWrong, 
       LOOCVCCVWrong, LOOCVCRCVWrong, LOOCVCR2CVWrong, 
       griddedCVs, griddedRCVs, griddedR2CVs, 
       griddedCVsWrong, griddedRCVsWrong, griddedR2CVsWrong, gridNs, 
       iter, rGRFargsTruth, rGRFargsSample, rGRFargsWrong, rho, 
       n1, n2, nx, ny, sigmaEpsSq1, sigmaEpsSq2, 
       alpha, beta, allSeeds, meanIW1, meanIW2, meanCVCIW1, meanCVCIW2, cors1IS, 
       cors2IS, cors1CVC, cors2CVC, trueVarW1, trueVarW2, 
       file=paste0("savedOutput/griddedCVtest2Datasets/n1", n1, "_n2", n2, "_rho", rho, "_iter", iter, ".RData"))
  
  list(trueMSE=trueMSE, wrongMSE=wrongMSE, LOOCVs=LOOCVs, LOOCVsWrong=LOOCVsWrong, 
       LOOCV=LOOCV, LOORCV=LOORCV, LOORCV=LOORCV, 
       LOOCVWrong=LOOCVWrong, LOOR2CVWrong=LOOR2CVWrong, LOOR2CVWrong=LOOR2CVWrong, 
       LOOISCV=LOOISCV, LOOISRCV=LOOISRCV, LOOISR2CV=LOOISR2CV, 
       LOOISPCV=LOOISPCV, LOOISPRCV=LOOISPRCV, LOOISPR2CV=LOOISPR2CV, 
       LOOISPCVWrong=LOOISPCVWrong, LOOISPRCVWrong=LOOISPRCVWrong, LOOISPR2CVWrong=LOOISPR2CVWrong, 
       LOOISCVWrong=LOOISCVWrong, LOOISRCVWrong=LOOISRCVWrong, LOOISR2CVWrong=LOOISR2CVWrong, 
       LOOVCCV=LOOVCCV, LOOVCRCV=LOOVCRCV, LOOVCR2CV=LOOVCR2CV, 
       LOOCVCCV=LOOCVCCV, LOOCVCRCV=LOOCVCRCV, LOOCVCR2CV=LOOCVCR2CV, 
       LOOVCCVWrong=LOOVCCVWrong, LOOVCRCVWrong=LOOVCRCVWrong, LOOVCR2CVWrong=LOOVCR2CVWrong, 
       LOOCVCCVWrong=LOOCVCCVWrong, LOOCVCRCVWrong=LOOCVCRCVWrong, LOOCVCR2CVWrong=LOOCVCR2CVWrong, 
       griddedCVs=griddedCVs, griddedRCVs=griddedRCVs, griddedR2CVs=griddedR2CVs, 
       griddedCVsWrong=griddedCVsWrong, griddedRCVsWrong=griddedRCVsWrong, griddedR2CVsWrong=griddedR2CVsWrong, 
       gridNs=gridNs, iter=iter, rGRFargsTruth=rGRFargsTruth, 
       rGRFargsSample=rGRFargsSample, rGRFargsWrong=rGRFargsWrong, 
       n1=n1, n2=n2, nx=nx, ny=ny, rho=rho, sigmaEpsSq1=sigmaEpsSq1, 
       sigmaEpsSq2=sigmaEpsSq2, allSeeds=allSeeds, 
       meanIW1=meanIW1, meanIW2=meanIW2, meanCVCIW1=meanCVCIW1, meanCVCIW2=meanCVCIW2, 
       cors1IS=cors1IS, cors2IS=cors2IS, cors1CVC=cors1CVC, cors2CVC=cors2CVC, 
       trueVarW1=trueVarW1, trueVarW2=trueVarW2)
}