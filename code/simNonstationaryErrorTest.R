# Test cross validation methods in the case of nonstationary error variance. 
# Fix a sampling distribution that is variable on unit square 
# based on GRF. Error variance will vary depending on a GRF potentially 
# correlated with the GRF determining the sampling rate.
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
griddedResTestIterNonstatError = function(rGRFargsTruth=NULL, rGRFargsMount=NULL, 
                                          propMount=.3, propSamplesMount=propMount/5, 
                                          n1=50, gridNs=2^(1:6), iter=1, 
                                          rGRFargsWrong1=NULL, rGRFargsWrong2=NULL, 
                                          nx=100, ny=100, sigmaEpsSqNonMount=.1^2, sigmaEpsSqMount=1^2, 
                                          sigmaEpsSqNonMountWrong1=.1^2, sigmaEpsSqMountWrong1=.1^2, 
                                          sigmaEpsSqNonMountWrong2=.15^2, sigmaEpsSqMountWrong2=.9^2, 
                                          allSeeds=123, printProgress=FALSE) {
  
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
  if(is.null(rGRFargsWrong1)) {
    rGRFargsWrong1 = rGRFargsTruth
    rGRFargsWrong1$cov.args$range = 1.25 * rGRFargsTruth$cov.args$range
    rGRFargsWrong1$cov.args$smoothness = 1.5
  }
  if(is.null(rGRFargsWrong2)) {
    rGRFargsWrong2 = rGRFargsWrong1
  }
  if(is.null(rGRFargsMount)) {
    rGRFargsMount = rGRFargsTruth
  }
  
  # 2. Simulate 1 GRF for responses ----
  #    on unit square
  truthGRF = do.call("rGRF", rGRFargsTruth)
  truth = truthGRF$truth
  locs = truthGRF$locs
  
  # 2. Simulate 1 GRF for mountains ----
  #    on unit square
  mountGRF = do.call("rGRF", rGRFargsMount)
  mountCov = mountGRF$truth
  isMount = mountCov > quantile(mountCov, .7, type=1) # exactly 30% of domain is mountainous
  
  # 3. Simulate sample distribution ----
  #    on mountainous and non-mountainous part of domain
  sampleRates = rep(1, nrow(locs))
  nMount = round(n1*propMount)
  nNonMount = n1 - nMount
  sampleRates[isMount] = nMount/propMount
  sampleRates[!isMount] = nNonMount/(1-propMount)
  sampleProbs = sampleRates * (1/sum(sampleRates))
  sampleProbsMount = sampleProbs[isMount]
  sampleProbsMount = sampleProbsMount/sum(sampleProbsMount)
  sampleProbsNonMount = sampleProbs[!isMount]
  sampleProbsNonMount = sampleProbsNonMount/sum(sampleProbsNonMount)
  
  # 4. Simulate observations ----
  #    in domain based on the sampling distn and the GRF 
  #    for responses (don't include nugget?)
  gridResX = 1/nx
  gridResY = 1/ny
  sampleIMount = sample((1:nrow(locs))[isMount], nMount, prob=sampleProbsMount, replace=TRUE)
  xsMount = locs[sampleIMount,] + cbind(runif(nMount, max=gridResX)-gridResX/2, 
                                runif(nMount, max=gridResY)-gridResY/2)
  ysMount = truth[sampleIMount] + rnorm(nMount, sd=sqrt(sigmaEpsSqMount))
  sampleINonMount = sample((1:nrow(locs))[!isMount], nNonMount, prob=sampleProbsNonMount, replace=TRUE)
  xsNonMount = locs[sampleINonMount,] + cbind(runif(nNonMount, max=gridResX)-gridResX/2, 
                                        runif(nNonMount, max=gridResY)-gridResY/2)
  ysNonMount = truth[sampleINonMount] + rnorm(nNonMount, sd=sqrt(sigmaEpsSqNonMount))
  
  xs = rbind(xsMount, xsNonMount)
  ys = c(ysMount, ysNonMount)
  mounts = c(rep(TRUE, nMount), rep(FALSE, nNonMount))
  sigmaEpsSqWrong1 = c(rep(sigmaEpsSqMountWrong1, nMount), rep(sigmaEpsSqNonMountWrong1, nNonMount))
  sigmaEpsSqWrong2 = c(rep(sigmaEpsSqMountWrong2, nMount), rep(sigmaEpsSqNonMountWrong2, nNonMount))
  sigmaEpsSqTrue = c(rep(sigmaEpsSqMount, nMount), rep(sigmaEpsSqNonMount, nNonMount))
  
  # 5. Get covariance matrices ----
  #    for sample and crossCov to all locs
  distMatXs = rdist(xs, xs)
  SigmaSample = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsTruth$cov.args$range, 
                               smoothness=rGRFargsTruth$cov.args$smoothness, distMat=distMatXs) * rGRFargsTruth$sigma^2
  SigmaSampleWrong1 = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsWrong1$cov.args$range, 
                                    smoothness=rGRFargsWrong1$cov.args$smoothness, distMat=distMatXs) * rGRFargsWrong1$sigma^2
  SigmaSampleWrong2 = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsWrong2$cov.args$range, 
                                     smoothness=rGRFargsWrong2$cov.args$smoothness, distMat=distMatXs) * rGRFargsWrong2$sigma^2
  SigmaSample = SigmaSample + diag(c(rep(sigmaEpsSq1, n1), rep(sigmaEpsSq2, n2)))
  SigmaSampleWrong1 = SigmaSampleWrong1 + diag(sigmaEpsSqWrong1)
  SigmaSampleWrong2 = SigmaSampleWrong2 + diag(sigmaEpsSqWrong2)
  distLocsToXs = rdist(locs, xs)
  distLocsToSample = distLocsToXs
  SigmaGridToSample = stationary.cov(locs, xs, Covariance="Matern", aRange=rGRFargsTruth$cov.args$range, 
                                     smoothness=rGRFargsTruth$cov.args$smoothness, distMat=distLocsToSample) * rGRFargsTruth$sigma^2
  SigmaGridToSampleWrong1 = stationary.cov(locs, xs, Covariance="Matern", aRange=rGRFargsWrong1$cov.args$range, 
                                          smoothness=rGRFargsWrong1$cov.args$smoothness, distMat=distLocsToSample) * rGRFargsWrong1$sigma^2
  SigmaGridToSampleWrong2 = stationary.cov(locs, xs, Covariance="Matern", aRange=rGRFargsWrong2$cov.args$range, 
                                           smoothness=rGRFargsWrong2$cov.args$smoothness, distMat=distLocsToSample) * rGRFargsWrong2$sigma^2
  
  # 6. for i in 1:number block resolutions: ----
  
  # determine grid cells associated with observations
  getCellI = function(thisLoc) {
    xI = match(FALSE, thisLoc[1] > highs)
    yI = match(FALSE, thisLoc[2] > highs)
    
    indMat[yI, xI]
  }
  
  griddedCVs = numeric(length(gridNs))
  griddedCVsWrong1 = numeric(length(gridNs))
  griddedCVsWrong2 = numeric(length(gridNs))
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
    cellIsSample = apply(xs, 1, getCellI) # only leave out observations from first dataset
    uniqueCellIs = sort(unique(cellIsSample))
    blockCVs = numeric(length(uniqueCellIs))
    blockCVsWrong1 = numeric(length(uniqueCellIs))
    blockCVsWrong2 = numeric(length(uniqueCellIs))
    ns = numeric(length(uniqueCellIs))
    
    # 8. Get gridded CV MSE ----
    for(j in 1:length(uniqueCellIs)) {
      if(printProgress) {
        print(paste0("Leaving out cell ", j, "/", length(uniqueCellIs)))
      }
      
      cellI = uniqueCellIs[j]
      
      # separate data in test/train
      isTest = cellIsSample == cellI
      testYs = ys[isTest]
      trainYs = ys[!isTest]
      
      # check number of observations in block
      ns[j] = sum(isTest)
      
      # predict test data for both true and wrong model
      SigmaAB = matrix(SigmaSample[isTest, !isTest], nrow=ns[j])
      SigmaBB = SigmaSample[!isTest, !isTest]
      SigmaAA = SigmaSample[isTest, isTest]
      
      SigmaABWrong1 = matrix(SigmaSampleWrong1[isTest, !isTest], nrow=ns[j])
      SigmaBBWrong1 = SigmaSampleWrong1[!isTest, !isTest]
      SigmaAAWrong1 = SigmaSampleWrong1[isTest, isTest]
      
      SigmaABWrong2 = matrix(SigmaSampleWrong2[isTest, !isTest], nrow=ns[j])
      SigmaBBWrong2 = SigmaSampleWrong2[!isTest, !isTest]
      SigmaAAWrong2 = SigmaSampleWrong2[isTest, isTest]
      
      if((ns[j] > 0) && (sum(!isTest) > 0)) {
        condDistn = condMeanMVN(SigmaAA=SigmaAA, SigmaAB=SigmaAB, SigmaBB=SigmaBB, 
                                ysB=trainYs, getFullCov=FALSE, getCondVar=TRUE)
        muAcondB = condDistn$muAcondB
        varAcondB = condDistn$varAcondB
        
        condDistnWrong1 = condMeanMVN(SigmaAA=SigmaAAWrong1, SigmaAB=SigmaABWrong1, SigmaBB=SigmaBBWrong1, 
                                     ysB=trainYs, getFullCov=FALSE, getCondVar=TRUE)
        muAcondBWrong = condDistnWrong1$muAcondB
        varAcondBWrong1 = condDistnWrong1$varAcondB
        
        condDistnWrong2 = condMeanMVN(SigmaAA=SigmaAAWrong2, SigmaAB=SigmaABWrong2, SigmaBB=SigmaBBWrong2, 
                                      ysB=trainYs, getFullCov=FALSE, getCondVar=TRUE)
        muAcondBWrong = condDistnWrong2$muAcondB
        varAcondBWrong2 = condDistnWrong2$varAcondB
        
        # # calculate MSE for the grid cell
        # blockCVs[j] = mean((testYs - muAcondB)^2)
        # blockCVsWrong1[j] = mean((testYs - muAcondBWrong1)^2)
        # blockCVsWrong2[j] = mean((testYs - muAcondBWrong2)^2)
        
        # calculate CRPS for the grid cell
        blockCVs[j] = crps(truth=testYs, est=muAcondB, my.var=varAcondB)
        blockCVsWrong1[j] = crps(truth=testYs, est=muAcondBWrong1, my.var=varAcondBWrong1)
        blockCVsWrong2[j] = crps(truth=testYs, est=muAcondBWrong2, my.var=varAcondBWrong2)
      } else {
        blockCVs[j] = NA
        blockCVsWrong[j] = NA
      }
    }
    
    # average over MSEs of each grid cell
    griddedCVs[i] = mean(blockCVs, na.rm=TRUE)
    griddedCVsWrong[i] = mean(blockCVsWrong, na.rm=TRUE)
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
    
    SigmaABWrong1 = matrix(SigmaSampleWrong1[isTest, !isTest], nrow=1)
    SigmaBBWrong1 = SigmaSampleWrong1[!isTest, !isTest]
    SigmaAAWrong1 = SigmaSampleWrong1[isTest, isTest]
    
    SigmaABWrong2 = matrix(SigmaSampleWrong2[isTest, !isTest], nrow=1)
    SigmaBBWrong2 = SigmaSampleWrong2[!isTest, !isTest]
    SigmaAAWrong2 = SigmaSampleWrong2[isTest, isTest]
    
    condDistn = condMeanMVN(SigmaAA=SigmaAA, SigmaAB=SigmaAB, SigmaBB=SigmaBB, 
                            ysB=trainYs, getFullCov=FALSE, getCondVar=TRUE)
    muAcondB = condDistn$muAcondB
    varAcondB = condDistn$varAcondB
    
    condDistnWrong1 = condMeanMVN(SigmaAA=SigmaAAWrong1, SigmaAB=SigmaABWrong1, SigmaBB=SigmaBBWrong1, 
                                 ysB=trainYs, getFullCov=FALSE, getCondVar=TRUE)
    muAcondBWrong1 = condDistnWrong1$muAcondB
    varAcondBWrong1 = condDistnWrong1$varAcondB
    
    condDistnWrong2 = condMeanMVN(SigmaAA=SigmaAAWrong2, SigmaAB=SigmaABWrong2, SigmaBB=SigmaBBWrong2, 
                                  ysB=trainYs, getFullCov=FALSE, getCondVar=TRUE)
    muAcondBWrong2 = condDistnWrong2$muAcondB
    varAcondBWrong2 = condDistnWrong2$varAcondB
    
    # calculate MSE for the grid cell
    LOOCVs[i] = crps(truth=testYs, est=muAcondB, my.var=varAcondB)
    LOOCVsWrong1[i] = crps(truth=testYs, est=muAcondBWrong1, my.var=varAcondBWrong1)
    LOOCVsWrong2[i] = crps(truth=testYs, est=muAcondBWrong2, my.var=varAcondBWrong2)
  }
  LOOCV = mean(LOOCVs)
  LOOCVWrong1 = mean(LOOCVsWrong1)
  LOOCVWrong2 = mean(LOOCVsWrong2)
  
  # get control variates
  minDistsGrid = apply(distLocsToXs, 1, min)
  minDistMu = mean(minDistsGrid)
  minDistsGrid2 = apply(distLocsToXs, 1, function(x) {sort(x)[2]})
  minDistMu2 = mean(minDistsGrid2)
  minDists = apply(distXsToXs, 1, function(x) {min(x[x != 0])}) # distance to nearest neighbor
  minDists2 = apply(distXsToXs, 1, function(x) {sort(x[x != 0])[2]}) # distance to second nearest neighbor
  
  # LOORCV = getRobustIWCV(scores=LOOCVs, IWs=rep(1, n1), controlVarMat=minDists1, 
  #                        controlVarMeans=minDist1Mu, wtdReg=FALSE)
  LOORCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=minDists, muf=minDistMu)[1]
  # LOORCVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=rep(1, n1), controlVarMat=minDists1, 
  #                             controlVarMeans=minDist1Mu, wtdReg=FALSE)
  LOORCVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=minDists, muf=minDistMu)[1]
  LOORCVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=minDists, muf=minDistMu)[1]
  
  # LOOR2CV = getRobustIWCV(scores=LOOCVs, IWs=rep(1, n1), controlVarMat=cbind(minDists1, minDists2), 
  #                        controlVarMeans=c(minDist1Mu, minDist2Mu), wtdReg=FALSE)
  LOOR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=cbind(minDists, minDists2), muf=c(minDist1Mu, minDist2Mu))[1]
  # LOOR2CVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=rep(1, n1), controlVarMat=cbind(minDists1, minDists2), 
  #                              controlVarMeans=c(minDist1Mu, minDist2Mu), wtdReg=FALSE)
  LOOR2CVWrong = getWeightedControlVarEst(Ys=LOOCVsWrong, Zf=cbind(minDists, minDists2), muf=c(minDist1Mu, minDistMu2))[1]
  
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
  LOOISPCVWrong1 = sum(LOOCVsWrong1 * (1/LOOISrates^p))/sum(1/LOOISrates^p)
  LOOISPCVWrong2 = sum(LOOCVsWrong2 * (1/LOOISrates^p))/sum(1/LOOISrates^p)
  
  # LOOISRCV = getRobustIWCV(scores=LOOCVs, IWs=1/LOOISrates, printVarRatio=FALSE, printRoot="LOOISR", controlVarMeans=1)
  # LOOISRCV = getRobustIWCV(scores=LOOCVs, IWs=1/LOOISrates, controlVarMat=LOOISrates, printVarRatio=FALSE, printRoot="LOOISR", controlVarMeans=1, wtdReg=TRUE)
  # LOOISRCV = getRobustIWCV(scores=LOOCVs, IWs=1/LOOISrates, controlVarMat=minDists1, printVarRatio=FALSE, printRoot="LOOISR", controlVarMeans=minDist1Mu, wtdReg=FALSE)
  # LOOISRCV = getGREGCV(scores=LOOCVs, IWs=1/LOOISrates, controlVarMat=cbind(1, minDists1), controlVarMeans=c(1, minDist1Mu), normalizeWeights=FALSE, 
  #                      shrinkWeights=FALSE)
  
  LOOISRCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates)[1]
  # pGregInfo = getPGREG(LOOCVs, IWs=1/LOOISrates, controlVarMat=minDists1, controlVarMeans=minDist1Mu)
  
  # LOOISRCVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMeans=1)
  # LOOISRCVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMat=LOOISrates, controlVarMeans=1, wtdReg=TRUE)
  # LOOISRCVWrong = getRobustIWCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMat=minDists1, controlVarMeans=minDist1Mu, wtdReg=FALSE)
  # LOOISRCVWrong = getGREGCV(scores=LOOCVsWrong, IWs=1/LOOISrates, controlVarMat=cbind(1, minDists1), controlVarMeans=c(1, minDist1Mu), normalizeWeights=FALSE, 
  #                           shrinkWeights=FALSE)
  LOOISRCVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates)[1]
  LOOISRCVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates)[1]
  
  LOOISR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=cbind(minDists, minDists2), muf=c(minDisMu, minDisMu2), Zg=cbind(1/LOOISrates), mug=c(1), ws=1/LOOISrates)[1]
  LOOISR2CVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=cbind(minDists, minDists2), muf=c(minDisMu, minDisMu2), Zg=cbind(1/LOOISrates), mug=c(1), ws=1/LOOISrates)[1]
  LOOISR2CVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=cbind(minDists, minDists2), muf=c(minDisMu, minDisMu2), Zg=cbind(1/LOOISrates), mug=c(1), ws=1/LOOISrates)[1]
  
  LOOISPRCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, shrinkWeights=TRUE)[1]
  LOOISPRCVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, shrinkWeights=TRUE)[1]
  LOOISPRCVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, shrinkWeights=TRUE)[1]
  LOOISPR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=c(minDists), muf=c(minDistMu), Zg=c(1/LOOISrates, LOOISrates2/LOOISrates), mug=c(1, 1), ws=1/LOOISrates, shrinkWeights=TRUE)[1]
  LOOISPR2CVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=c(1/LOOISrates), mug=c(1, 1), ws=1/LOOISrates, shrinkWeights=TRUE)[1]
  LOOISPR2CVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=c(1/LOOISrates), mug=c(1, 1), ws=1/LOOISrates, shrinkWeights=TRUE)[1]
  
  cors1IS = cor(1/LOOISrates, LOOCVs/LOOISrates)
  
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
  
  LOOVCRCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=c(minDists1), muf=c(minDist1Mu), ws=1/rateEsts, shrinkWeights=FALSE)[1]
  LOOVCRCVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=FALSE)[1]
  LOOVCRCVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=FALSE)[1]
  
  LOOVCR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), mug=1, ws=1/rateEsts, shrinkWeights=FALSE)[1]
  LOOVCR2CVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), mug=1, ws=1/rateEsts, shrinkWeights=FALSE)[1]
  LOOVCR2CVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), mug=1, ws=1/rateEsts, shrinkWeights=FALSE)[1]
  
  # 12. Calculate model Scores ----
  condDistn = condMeanMVN(SigmaAA=rep(rGRFargsTruth$sigma^2, nrow(locs)), SigmaAB=SigmaGridToSample, SigmaBB=SigmaSample, 
                          ysB=ys, getFullCov=FALSE, getCondVar=TRUE)
  muAcondB = condDistn$muAcondB
  varAcondB = condDistn$varAcondB
  condDistn = condMeanMVN(SigmaAA=rep(rGRFargsWrong1$sigma^2, nrow(locs)), SigmaAB=SigmaGridToSampleWrong1, SigmaBB=SigmaSampleWrong1, 
                          ysB=ys, getFullCov=FALSE, getCondVar=TRUE)
  muAcondBwrong1 = condDistn$muAcondB
  varAcondBwrong1 = condDistn$varAcondB
  condDistn = condMeanMVN(SigmaAA=rep(rGRFargsWrong2$sigma^2, nrow(locs)), SigmaAB=SigmaGridToSampleWrong2, SigmaBB=SigmaSampleWrong2, 
                          ysB=ys, getFullCov=FALSE, getCondVar=TRUE)
  muAcondBwrong2 = condDistn$muAcondB
  varAcondBwrong2 = condDistn$varAcondB
  
  # calculate Scores for the grid cells
  # trueMSE = mean((truth - muAcondB)^2) + sigmaEpsSq1
  # wrongMSE = mean((truth - muAcondBwrong)^2) + sigmaEpsSq1 # add sigmaEpsSq1, the true error variance, not sigmaEpsSq2
  trueMSE = crps(truth=truth, est=muAcondB, my.var=varAcondB + sigmaEpsSqTrue)
  wrongMSE1 = crps(truth=truth, est=muAcondBWrong1, my.var=varAcondBWrong1 + sigmaEpsSqWrong1)
  wrongMSE2 = crps(truth=truth, est=muAcondBWrong2, my.var=varAcondBWrong2 + sigmaEpsSqWrong2)
  
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
  trueVarW = -1 + mean(1/sampleProbs)
  
  # 13. Save result ----
  save(trueMSE, wrongMSE, LOOCVs, LOOCVsWrong1, LOOCVsWrong2, 
       LOOCV, LOORCV, LOOR2CV, 
       LOOCVWrong1, LOORCVWrong1, LOOR2CVWrong1, 
       LOOCVWrong2, LOORCVWrong2, LOOR2CVWrong2, 
       LOOISCV, LOOISRCV, LOOISR2CV, 
       LOOISPCV, LOOISPRCV, LOOISPR2CV, 
       LOOISPCVWrong1, LOOISPRCVWrong1, LOOISPR2CVWrong1, 
       LOOISCVWrong1, LOOISRCVWrong1, LOOISR2CVWrong1, 
       LOOISPCVWrong2, LOOISPRCVWrong2, LOOISPR2CVWrong2, 
       LOOISCVWrong2, LOOISRCVWrong2, LOOISR2CVWrong2, 
       LOOVCCV, LOOVCRCV, LOOVCR2CV, 
       LOOCVCCV, LOOCVCRCV, LOOCVCR2CV, 
       LOOVCCVWrong1, LOOVCRCVWrong1, LOOVCR2CVWrong1, 
       LOOCVCCVWrong1, LOOCVCRCVWrong1, LOOCVCR2CVWrong1, 
       LOOVCCVWrong2, LOOVCRCVWrong2, LOOVCR2CVWrong2, 
       LOOCVCCVWrong2, LOOCVCRCVWrong2, LOOCVCR2CVWrong2, 
       griddedCVs, griddedRCVs, griddedR2CVs, 
       griddedCVsWrong1, griddedRCVsWrong1, griddedR2CVsWrong1, 
       griddedCVsWrong2, griddedRCVsWrong2, griddedR2CVsWrong2, gridNs, 
       iter, rGRFargsTruth, rGRFargsWrong1, rGRFargsWrong2, 
       n1, nx, ny, propMount, propSamplesMount, 
       sigmaEpsSqNonMount, sigmaEpsSqMount, 
       sigmaEpsSqNonMountWrong1, sigmaEpsSqMountWrong1, 
       sigmaEpsSqNonMountWrong2, sigmaEpsSqMountWrong2, 
       allSeeds, trueVarW1, trueVarW2, 
       file=paste0("savedOutput/griddedCVtestNonstatErr/n1", n1, "_pMount", propMount, "_pSMount", propSamplesMount, 
                   "_s2M", sigmaEpsSqMount, "_", sigmaEpsSqMountWrong1, "_", sigmaEpsSqMountWrong2, 
                   "_s2NM", sigmaEpsSqNonMount, "_", sigmaEpsSqNonMountWrong1, "_", sigmaEpsSqNonMountWrong2, 
                   "_iter", iter, ".RData"))
  
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