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
                                          n1=50, gridNs=2^(1:6), Ks=c(9, 25), iter=1, 
                                          rGRFargsWrong1=NULL, rGRFargsWrong2=NULL, 
                                          nx=100, ny=100, sigmaEpsSqNonMount=.1^2, sigmaEpsSqMount=1^2, 
                                          sigmaEpsSqNonMountWrong1=.1^2, sigmaEpsSqMountWrong1=.1^2, 
                                          sigmaEpsSqNonMountWrong2=.15^2, sigmaEpsSqMountWrong2=.9^2, 
                                          allSeeds=123, printProgress=FALSE, printPhat=FALSE) {
  
  if(iter > 1) {
    currT = proc.time()[3]
    tLeft = estTimeLeft(startT, currT, finishedIter=iter-1, totalIter=totalIter)/60 # convert to minutes
    print(paste0("iteration ", iter, "/", length(allSeeds), ". Est minutes left: ", round(tLeft, digits=1)))
  } else {
    print(paste0("iteration ", iter, "/", length(allSeeds)))
  }
  
  
  doLOO = TRUE
  if(n1 > 50) {
    gridNs = c(3, 5)
    
    if(n1 > 500) {
      doLOO = FALSE
      if(n1 > 2000) {
        Ks = Ks[Ks <= 10]
      }
    }
  }
  
  
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
  
  # 2.5. Simulate 1 GRF for mountains ----
  #    on unit square
  mountGRF = do.call("rGRF", rGRFargsMount)
  mountCov = mountGRF$truth
  isMount = mountCov > quantile(mountCov, .7, type=1) # exactly 30% of domain is mountainous
  
  # 3. Simulate sample distribution ----
  #    on mountainous and non-mountainous part of domain
  sampleRates = rep(1, nrow(locs))
  nMount = round(n1*propSamplesMount)
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
  sampleI = c(sampleIMount, sampleINonMount)
  mounts = c(rep(TRUE, nMount), rep(FALSE, nNonMount))
  sigmaEpsSqWrong1 = c(rep(sigmaEpsSqMountWrong1, nMount), rep(sigmaEpsSqNonMountWrong1, nNonMount))
  sigmaEpsSqWrong2 = c(rep(sigmaEpsSqMountWrong2, nMount), rep(sigmaEpsSqNonMountWrong2, nNonMount))
  sigmaEpsSqTrue = c(rep(sigmaEpsSqMount, nMount), rep(sigmaEpsSqNonMount, nNonMount))
  
  sigmaEpsSqTrueLoc = rep(sigmaEpsSqNonMount, nrow(locs))
  sigmaEpsSqTrueLoc[isMount] = sigmaEpsSqMount
  sigmaEpsSqWrong1Loc = rep(sigmaEpsSqNonMountWrong1, nrow(locs))
  sigmaEpsSqWrong1Loc[isMount] = sigmaEpsSqMountWrong1
  sigmaEpsSqWrong2Loc = rep(sigmaEpsSqNonMountWrong2, nrow(locs))
  sigmaEpsSqWrong2Loc[isMount] = sigmaEpsSqMountWrong2
  
  # 4.5. Set folds ----
  foldI = matrix(nrow=n1, ncol=length(Ks))
  for(i in 1:length(Ks)) {
    K = Ks[i]
    numPerFold = diff(floor(seq(0, n1, by=n1/K)))
    indsLeft = 1:n1
    
    for(j in 1:length(numPerFold)) {
      thisFoldN = numPerFold[j]
      thisFoldI = sample(indsLeft, thisFoldN, replace=FALSE)
      foldI[thisFoldI,i] = j
      indsLeft = indsLeft[!(indsLeft %in% thisFoldI)]
    }
  }
  
  # 5. Get covariance matrices ----
  #    for sample and crossCov to all locs
  distMatXs = rdist(xs, xs)
  SigmaSample = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsTruth$cov.args$range, 
                               smoothness=rGRFargsTruth$cov.args$smoothness, distMat=distMatXs) * rGRFargsTruth$sigma^2
  SigmaSampleWrong1 = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsWrong1$cov.args$range, 
                                    smoothness=rGRFargsWrong1$cov.args$smoothness, distMat=distMatXs) * rGRFargsWrong1$sigma^2
  SigmaSampleWrong2 = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsWrong2$cov.args$range, 
                                     smoothness=rGRFargsWrong2$cov.args$smoothness, distMat=distMatXs) * rGRFargsWrong2$sigma^2
  SigmaSample = SigmaSample + diag(sigmaEpsSqTrue)
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
    
    # 8. Get gridded CV CRPS ----
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
        muAcondBWrong1 = condDistnWrong1$muAcondB
        varAcondBWrong1 = condDistnWrong1$varAcondB
        
        condDistnWrong2 = condMeanMVN(SigmaAA=SigmaAAWrong2, SigmaAB=SigmaABWrong2, SigmaBB=SigmaBBWrong2, 
                                      ysB=trainYs, getFullCov=FALSE, getCondVar=TRUE)
        muAcondBWrong2 = condDistnWrong2$muAcondB
        varAcondBWrong2 = condDistnWrong2$varAcondB
        
        # # calculate MSE for the grid cell
        # blockCVs[j] = mean((testYs - muAcondB)^2)
        # blockCVsWrong1[j] = mean((testYs - muAcondBWrong1)^2)
        # blockCVsWrong2[j] = mean((testYs - muAcondBWrong2)^2)
        
        # calculate CRPS for the grid cell
        blockCVs[j] = crps(truth=testYs, est=muAcondB, est.var=varAcondB)
        blockCVsWrong1[j] = crps(truth=testYs, est=muAcondBWrong1, est.var=varAcondBWrong1)
        blockCVsWrong2[j] = crps(truth=testYs, est=muAcondBWrong2, est.var=varAcondBWrong2)
      } else {
        blockCVs[j] = NA
        blockCVsWrong1[j] = NA
        blockCVsWrong2[j] = NA
      }
    }
    
    # average over MSEs of each grid cell
    griddedCVs[i] = mean(blockCVs, na.rm=TRUE)
    griddedCVsWrong1[i] = mean(blockCVsWrong1, na.rm=TRUE)
    griddedCVsWrong2[i] = mean(blockCVsWrong2, na.rm=TRUE)
  }
  
  # 6. for i in 1:number fold resolutions: ----
  
  KFoldCVs = matrix(nrow=n1, ncol=length(Ks))
  KFoldCVsWrong1 = matrix(nrow=n1, ncol=length(Ks))
  KFoldCVsWrong2 = matrix(nrow=n1, ncol=length(Ks))
  
  KFoldCV = numeric(length(Ks))
  KFoldCVWrong1 = numeric(length(Ks))
  KFoldCVWrong2 = numeric(length(Ks))
  
  KFoldCV = numeric(length(Ks))
  KFoldCVWrong1 = numeric(length(Ks))
  KFoldCVWrong2 = numeric(length(Ks))
  
  KFoldRCV = numeric(length(Ks))
  KFoldRCVWrong1 = numeric(length(Ks))
  KFoldRCVWrong2 = numeric(length(Ks))
  
  KFoldR2CV = numeric(length(Ks))
  KFoldR2CVWrong1 = numeric(length(Ks))
  KFoldR2CVWrong2 = numeric(length(Ks))
  
  cellArea = 1/(nx*ny)
  KFoldISrates = numeric(length(Ks))
  KFoldISCV = numeric(length(Ks))
  KFoldISCVWrong1 = numeric(length(Ks))
  KFoldISCVWrong2 = numeric(length(Ks))
  
  KFoldISRCV = numeric(length(Ks))
  KFoldISRCVWrong1 = numeric(length(Ks))
  KFoldISRCVWrong2 = numeric(length(Ks))
  
  KFoldISR2CV = numeric(length(Ks))
  KFoldISR2CVWrong1 = numeric(length(Ks))
  KFoldISR2CVWrong2 = numeric(length(Ks))
  
  KFoldISPCV = numeric(length(Ks))
  KFoldISPCVWrong1 = numeric(length(Ks))
  KFoldISPCVWrong2 = numeric(length(Ks))
  
  KFoldISPRCV = numeric(length(Ks))
  KFoldISPRCVWrong1 = numeric(length(Ks))
  KFoldISPRCVWrong2 = numeric(length(Ks))
  
  KFoldISPR2CV = numeric(length(Ks))
  KFoldISPR2CVWrong1 = numeric(length(Ks))
  KFoldISPR2CVWrong2 = numeric(length(Ks))
  
  KFoldVCCV = numeric(length(Ks))
  KFoldVCCVWrong1 = numeric(length(Ks))
  KFoldVCCVWrong2 = numeric(length(Ks))
  
  KFoldVCRCV = numeric(length(Ks))
  KFoldVCRCVWrong1 = numeric(length(Ks))
  KFoldVCRCVWrong2 = numeric(length(Ks))
  
  KFoldVCR2CV = numeric(length(Ks))
  KFoldVCR2CVWrong1 = numeric(length(Ks))
  KFoldVCR2CVWrong2 = numeric(length(Ks))
  
  KFoldVCPCV = numeric(length(Ks))
  KFoldVCPCVWrong1 = numeric(length(Ks))
  KFoldVCPCVWrong2 = numeric(length(Ks))
  
  KFoldVCPRCV = numeric(length(Ks))
  KFoldVCPRCVWrong1 = numeric(length(Ks))
  KFoldVCPRCVWrong2 = numeric(length(Ks))
  
  KFoldVCPR2CV = numeric(length(Ks))
  KFoldVCPR2CVWrong1 = numeric(length(Ks))
  KFoldVCPR2CVWrong2 = numeric(length(Ks))
  for(i in 1:length(Ks)) {
    K = Ks[i]
    if(printProgress) {
      print(paste0("K: ", K, " (", i, "/", length(Ks), ")"))
    }
    
    # 7. Group data by fold ----
    thisFoldI = foldI[,i]
    
    # 8. Get K fold CV CRPS ----
    for(j in 1:K) {
      if(printProgress) {
        print(paste0("Leaving out fold ", j, "/", K))
      }
      
      # separate data in test/train
      isTest = thisFoldI == j
      testYs = ys[isTest]
      trainYs = ys[!isTest]
      
      # predict test data for both true and wrong model
      SigmaAB = matrix(SigmaSample[isTest, !isTest], nrow=sum(isTest))
      SigmaBB = SigmaSample[!isTest, !isTest]
      SigmaAA = SigmaSample[isTest, isTest]
      
      SigmaABWrong1 = matrix(SigmaSampleWrong1[isTest, !isTest], nrow=sum(isTest))
      SigmaBBWrong1 = SigmaSampleWrong1[!isTest, !isTest]
      SigmaAAWrong1 = SigmaSampleWrong1[isTest, isTest]
      
      SigmaABWrong2 = matrix(SigmaSampleWrong2[isTest, !isTest], nrow=sum(isTest))
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
      
      # # calculate MSE for the fold
      # foldCVs[i] = mean((testYs - muAcondB)^2)
      # foldCVsWrong1[i] = mean((testYs - muAcondBWrong1)^2)
      # foldCVsWrong2[i] = mean((testYs - muAcondBWrong2)^2)
      
      # calculate CRPS for each observation in the fold
      KFoldCVs[isTest,i] = crps(truth=testYs, est=muAcondB, est.var=varAcondB, getAverage = FALSE)
      KFoldCVsWrong1[isTest,i] = crps(truth=testYs, est=muAcondBWrong1, est.var=varAcondBWrong1, getAverage = FALSE)
      KFoldCVsWrong2[isTest,i] = crps(truth=testYs, est=muAcondBWrong2, est.var=varAcondBWrong2, getAverage = FALSE)
    }
    
    # average over CRPSs of each observation over all folds, for each value of K
    KFoldCV[i] = mean(KFoldCVs[,i], na.rm=TRUE)
    KFoldCVWrong1[i] = mean(KFoldCVsWrong1[,i], na.rm=TRUE)
    KFoldCVWrong2[i] = mean(KFoldCVsWrong2[,i], na.rm=TRUE)
    
    # get control variates
    minDistsGrid = apply(distLocsToXs, 1, min)
    minDistMu = mean(minDistsGrid)
    minDistsGrid2 = apply(distLocsToXs, 1, function(x) {sort(x)[2]})
    minDistMu2 = mean(minDistsGrid2)
    minDists = apply(distMatXs, 1, function(x) {min(x[x != 0])}) # distance to nearest neighbor
    minDists2 = apply(distMatXs, 1, function(x) {sort(x[x != 0])[2]}) # distance to second nearest neighbor
    
    KFoldRCV[i] = getWeightedControlVarEst(Ys=KFoldCVs[,i], Zf=minDists, muf=minDistMu, printPhat=printPhat)[1]
    KFoldRCVWrong1[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong1[,i], Zf=minDists, muf=minDistMu, printPhat=printPhat)[1]
    KFoldRCVWrong2[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong2[,i], Zf=minDists, muf=minDistMu, printPhat=printPhat)[1]
    
    KFoldR2CV[i] = getWeightedControlVarEst(Ys=KFoldCVs[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), printPhat=printPhat)[1]
    KFoldR2CVWrong1[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong1[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), printPhat=printPhat)[1]
    KFoldR2CVWrong2[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong2[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), printPhat=printPhat)[1]
    
    # 10. Calculate IS-CV CRPS ----
    cellArea = 1/(nx*ny)
    LOOISrates = (sampleProbs[sampleI]/sum(sampleProbs)) / cellArea # convert from areal rate to rate density
    # KFoldISCV = weighted.mean(KFoldCVs[,i], w=cellArea/KFoldISrates)
    KFoldISCV[i] = sum(KFoldCVs[,i] * (1/LOOISrates))/n1
    KFoldISCVWrong1[i] = sum(KFoldCVsWrong1[,i] * (1/LOOISrates))/n1
    KFoldISCVWrong2[i] = sum(KFoldCVsWrong2[,i] * (1/LOOISrates))/n1
    
    # estimate ratio of CRPSs
    out = getP(KFoldCVs[,i], 1/LOOISrates)
    p = out$p
    KFoldISPCV[i] = sum(KFoldCVs[,i] * (1/LOOISrates^p))/sum(1/LOOISrates^p)
    out = getP(KFoldCVsWrong1[,i], 1/LOOISrates)
    p = out$p
    KFoldISPCVWrong1[i] = sum(KFoldCVsWrong1[,i] * (1/LOOISrates^p))/sum(1/LOOISrates^p)
    out = getP(KFoldCVsWrong2[,i], 1/LOOISrates)
    p = out$p
    KFoldISPCVWrong2[i] = sum(KFoldCVsWrong2[,i] * (1/LOOISrates^p))/sum(1/LOOISrates^p)
    
    KFoldISRCV[i] = getWeightedControlVarEst(Ys=KFoldCVs[,i], Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, printPhat=printPhat)[1]
    KFoldISRCVWrong1[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong1[,i], Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, printPhat=printPhat)[1]
    KFoldISRCVWrong2[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong2[,i], Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, printPhat=printPhat)[1]
    
    KFoldISR2CV[i] = getWeightedControlVarEst(Ys=KFoldCVs[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=cbind(1/LOOISrates), mug=c(1), ws=1/LOOISrates, printPhat=printPhat)[1]
    KFoldISR2CVWrong1[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong1[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=cbind(1/LOOISrates), mug=c(1), ws=1/LOOISrates, printPhat=printPhat)[1]
    KFoldISR2CVWrong2[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong2[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=cbind(1/LOOISrates), mug=c(1), ws=1/LOOISrates, printPhat=printPhat)[1]
    
    KFoldISPRCV[i] = getWeightedControlVarEst(Ys=KFoldCVs[,i], Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, shrinkWeights=TRUE, printPhat=printPhat)[1]
    KFoldISPRCVWrong1[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong1[,i], Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, shrinkWeights=TRUE, printPhat=printPhat)[1]
    KFoldISPRCVWrong2[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong2[,i], Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, shrinkWeights=TRUE, printPhat=printPhat)[1]
    
    KFoldISPR2CV[i] = getWeightedControlVarEst(Ys=KFoldCVs[,i], Zf=c(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=c(1/LOOISrates), mug=c(1), ws=1/LOOISrates, shrinkWeights=TRUE, printPhat=printPhat)[1]
    KFoldISPR2CVWrong1[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong1[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=c(1/LOOISrates), mug=c(1), ws=1/LOOISrates, shrinkWeights=TRUE, printPhat=printPhat)[1]
    KFoldISPR2CVWrong2[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong2[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=c(1/LOOISrates), mug=c(1), ws=1/LOOISrates, shrinkWeights=TRUE, printPhat=printPhat)[1]
    
    # cors1IS = cor(1/LOOISrates, KFoldCVs[,i]/LOOISrates)
    
    # 11. Calculate K fold VC-CV MSE ----
    domainPoly = rbind(c(0, 0), 
                       c(0, 1), 
                       c(1, 1), 
                       c(1, 0), 
                       c(0, 0))
    vcellInfo = getVCellAreas(xs, domainPoly=domainPoly)
    vcellArea = vcellInfo$area
    rateEsts = 1/vcellArea
    rateEsts = rateEsts/n1 # divide by n to get rate for a single pt draw instead of n pts
    
    KFoldVCCV[i] = sum(KFoldCVs[,i] * (1/rateEsts))/n1
    KFoldVCCVWrong1[i] = sum(KFoldCVsWrong1[,i] * (1/rateEsts))/n1
    KFoldVCCVWrong2[i] = sum(KFoldCVsWrong2[,i] * (1/rateEsts))/n1
    
    out = getP(KFoldCVs[,i], 1/rateEsts)
    p = out$p
    KFoldVCPCV[i] = sum(KFoldCVs[,i] * (1/rateEsts^p))/sum(1/rateEsts^p)
    out = getP(KFoldCVsWrong1[,i], 1/rateEsts)
    p = out$p
    KFoldVCPCVWrong1[i] = sum(KFoldCVsWrong1[,i] * (1/rateEsts^p))/sum(1/rateEsts^p)
    out = getP(KFoldCVsWrong2[,i], 1/rateEsts)
    p = out$p
    KFoldVCPCVWrong2[i] = sum(KFoldCVsWrong2[,i] * (1/rateEsts^p))/sum(1/rateEsts^p)
    
    # I would add 1/rateEsts as an argument Zg with mug=1, but sum(1/rateEsts) == 1 for VC rate estimation
    KFoldVCRCV[i] = getWeightedControlVarEst(Ys=KFoldCVs[,i], Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=FALSE, printPhat=printPhat)[1]
    KFoldVCRCVWrong1[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong1[,i], Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=FALSE, printPhat=printPhat)[1]
    KFoldVCRCVWrong2[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong2[,i], Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=FALSE, printPhat=printPhat)[1]
    
    KFoldVCR2CV[i] = getWeightedControlVarEst(Ys=KFoldCVs[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), ws=1/rateEsts, shrinkWeights=FALSE, printPhat=printPhat)[1]
    KFoldVCR2CVWrong1[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong1[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), ws=1/rateEsts, shrinkWeights=FALSE, printPhat=printPhat)[1]
    KFoldVCR2CVWrong2[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong2[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), ws=1/rateEsts, shrinkWeights=FALSE, printPhat=printPhat)[1]
    
    KFoldVCPRCV[i] = getWeightedControlVarEst(Ys=KFoldCVs[,i], Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=TRUE, printPhat=printPhat)[1]
    KFoldVCPRCVWrong1[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong1[,i], Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=TRUE, printPhat=printPhat)[1]
    KFoldVCPRCVWrong2[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong2[,i], Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=TRUE, printPhat=printPhat)[1]
    
    KFoldVCPR2CV[i] = getWeightedControlVarEst(Ys=KFoldCVs[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), ws=1/rateEsts, shrinkWeights=TRUE, printPhat=printPhat)[1]
    KFoldVCPR2CVWrong1[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong1[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), ws=1/rateEsts, shrinkWeights=TRUE, printPhat=printPhat)[1]
    KFoldVCPR2CVWrong2[i] = getWeightedControlVarEst(Ys=KFoldCVsWrong2[,i], Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), ws=1/rateEsts, shrinkWeights=TRUE, printPhat=printPhat)[1]
  }
  
  if(doLOO) {
    # 9. get LOO-CV CRPS ----
    LOOCVs = numeric(n1)
    LOOCVsWrong1 = numeric(n1)
    LOOCVsWrong2 = numeric(n1)
    for(i in 1:n1) {
      if(printProgress) {
        print(paste0("Leaving out obs ", i, "/", n1))
      }
      
      # separate data in test/train
      isTest = (1:(n1)) == i
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
      LOOCVs[i] = crps(truth=testYs, est=muAcondB, est.var=varAcondB)
      LOOCVsWrong1[i] = crps(truth=testYs, est=muAcondBWrong1, est.var=varAcondBWrong1)
      LOOCVsWrong2[i] = crps(truth=testYs, est=muAcondBWrong2, est.var=varAcondBWrong2)
    }
    LOOCV = mean(LOOCVs)
    LOOCVWrong1 = mean(LOOCVsWrong1)
    LOOCVWrong2 = mean(LOOCVsWrong2)
    
    # get control variates
    minDistsGrid = apply(distLocsToXs, 1, min)
    minDistMu = mean(minDistsGrid)
    minDistsGrid2 = apply(distLocsToXs, 1, function(x) {sort(x)[2]})
    minDistMu2 = mean(minDistsGrid2)
    minDists = apply(distMatXs, 1, function(x) {min(x[x != 0])}) # distance to nearest neighbor
    minDists2 = apply(distMatXs, 1, function(x) {sort(x[x != 0])[2]}) # distance to second nearest neighbor
    
    LOORCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=minDists, muf=minDistMu, printPhat=printPhat)[1]
    LOORCVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=minDists, muf=minDistMu, printPhat=printPhat)[1]
    LOORCVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=minDists, muf=minDistMu, printPhat=printPhat)[1]
    
    LOOR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), printPhat=printPhat)[1]
    LOOR2CVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), printPhat=printPhat)[1]
    LOOR2CVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), printPhat=printPhat)[1]
    
    # 10. Calculate LOOIS-CV CRPS ----
    cellArea = 1/(nx*ny)
    LOOISrates = (sampleProbs[sampleI]/sum(sampleProbs)) / cellArea # convert from single obs sample prob to single obs sample rate
    # LOOISCV = weighted.mean(LOOCVs, w=cellArea/LOOISrates)
    LOOISCV = sum(LOOCVs * (1/LOOISrates))/n1
    LOOISCVWrong1 = sum(LOOCVsWrong1 * (1/LOOISrates))/n1
    LOOISCVWrong2 = sum(LOOCVsWrong2 * (1/LOOISrates))/n1
    
    # estimate ratio of MSEs
    out = getP(LOOCVs, 1/LOOISrates)
    p = out$p
    LOOISPCV = sum(LOOCVs * (1/LOOISrates^p))/sum(1/LOOISrates^p)
    out = getP(LOOCVsWrong1, 1/LOOISrates)
    p = out$p
    LOOISPCVWrong1 = sum(LOOCVsWrong1 * (1/LOOISrates^p))/sum(1/LOOISrates^p)
    out = getP(LOOCVsWrong2, 1/LOOISrates)
    p = out$p
    LOOISPCVWrong2 = sum(LOOCVsWrong2 * (1/LOOISrates^p))/sum(1/LOOISrates^p)
    
    LOOISRCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, printPhat=printPhat)[1]
    LOOISRCVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, printPhat=printPhat)[1]
    LOOISRCVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, printPhat=printPhat)[1]
    
    LOOISR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=cbind(1/LOOISrates), mug=c(1), ws=1/LOOISrates, printPhat=printPhat)[1]
    LOOISR2CVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=cbind(1/LOOISrates), mug=c(1), ws=1/LOOISrates, printPhat=printPhat)[1]
    LOOISR2CVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=cbind(1/LOOISrates), mug=c(1), ws=1/LOOISrates, printPhat=printPhat)[1]
    
    LOOISPRCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, shrinkWeights=TRUE, printPhat=printPhat)[1]
    LOOISPRCVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, shrinkWeights=TRUE, printPhat=printPhat)[1]
    LOOISPRCVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=minDists, muf=minDistMu, Zg=1/LOOISrates, mug=1, ws=1/LOOISrates, shrinkWeights=TRUE, printPhat=printPhat)[1]
    
    LOOISPR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=c(minDists), muf=c(minDistMu), Zg=c(1/LOOISrates), mug=c(1), ws=1/LOOISrates, shrinkWeights=TRUE, printPhat=printPhat)[1]
    LOOISPR2CVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=c(1/LOOISrates), mug=c(1), ws=1/LOOISrates, shrinkWeights=TRUE, printPhat=printPhat)[1]
    LOOISPR2CVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), Zg=c(1/LOOISrates), mug=c(1), ws=1/LOOISrates, shrinkWeights=TRUE, printPhat=printPhat)[1]
    
    # cors1IS = cor(1/LOOISrates, LOOCVs/LOOISrates)
    
    # 11. Calculate LOOVC-CV CRPS ----
    domainPoly = rbind(c(0, 0), 
                       c(0, 1), 
                       c(1, 1), 
                       c(1, 0), 
                       c(0, 0))
    vcellInfo = getVCellAreas(xs, domainPoly=domainPoly)
    vcellArea = vcellInfo$area
    rateEsts = 1/vcellArea
    rateEsts = rateEsts/n1 # divide by n to get rate for a single pt draw instead of n pts
    
    LOOVCCV = sum(LOOCVs * (1/rateEsts))/n1
    LOOVCCVWrong1 = sum(LOOCVsWrong1 * (1/rateEsts))/n1
    LOOVCCVWrong2 = sum(LOOCVsWrong2 * (1/rateEsts))/n1
    
    out = getP(LOOCVs, 1/rateEsts)
    p = out$p
    LOOVCPCV = sum(LOOCVs * (1/rateEsts^p))/sum(1/rateEsts^p)
    out = getP(LOOCVsWrong1, 1/LOOISrates)
    p = out$p
    LOOVCPCVWrong1 = sum(LOOCVsWrong1 * (1/rateEsts^p))/sum(1/rateEsts^p)
    out = getP(LOOCVsWrong2, 1/LOOISrates)
    p = out$p
    LOOVCPCVWrong2 = sum(LOOCVsWrong2 * (1/rateEsts^p))/sum(1/rateEsts^p)
    
    LOOVCRCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=FALSE, printPhat=printPhat)[1]
    LOOVCRCVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=FALSE, printPhat=printPhat)[1]
    LOOVCRCVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=FALSE, printPhat=printPhat)[1]
    
    LOOVCR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), ws=1/rateEsts, shrinkWeights=FALSE, printPhat=printPhat)[1]
    LOOVCR2CVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), ws=1/rateEsts, shrinkWeights=FALSE, printPhat=printPhat)[1]
    LOOVCR2CVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), ws=1/rateEsts, shrinkWeights=FALSE, printPhat=printPhat)[1]
    
    LOOVCPRCV = getWeightedControlVarEst(Ys=LOOCVs, Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=TRUE, printPhat=printPhat)[1]
    LOOVCPRCVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=TRUE, printPhat=printPhat)[1]
    LOOVCPRCVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=c(minDists), muf=c(minDistMu), ws=1/rateEsts, shrinkWeights=TRUE, printPhat=printPhat)[1]
    
    LOOVCPR2CV = getWeightedControlVarEst(Ys=LOOCVs, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), ws=1/rateEsts, shrinkWeights=TRUE, printPhat=printPhat)[1]
    LOOVCPR2CVWrong1 = getWeightedControlVarEst(Ys=LOOCVsWrong1, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), ws=1/rateEsts, shrinkWeights=TRUE, printPhat=printPhat)[1]
    LOOVCPR2CVWrong2 = getWeightedControlVarEst(Ys=LOOCVsWrong2, Zf=cbind(minDists, minDists2), muf=c(minDistMu, minDistMu2), ws=1/rateEsts, shrinkWeights=TRUE, printPhat=printPhat)[1]
  } else {
    LOOCVs[i] = numeric(n1)
    LOOCVsWrong1[i] = numeric(n1)
    LOOCVsWrong2[i] = numeric(n1)
    
    LOOCV = NULL
    LOOCVWrong1 = NULL
    LOOCVWrong2 = NULL
    
    LOORCV = NULL
    LOORCVWrong1 = NULL
    LOORCVWrong2 = NULL
    
    LOOR2CV = NULL
    LOOR2CVWrong1 = NULL
    LOOR2CVWrong2 = NULL
    
    cellArea = 1/(nx*ny)
    LOOISrates = NULL
    LOOISCV = NULL
    LOOISCVWrong = NULL
    
    LOOISPCV = NULL
    LOOISPCVWrong1 = NULL
    LOOISPCVWrong2 = NULL
    
    LOOISRCV = NULL
    LOOISRCVWrong1 = NULL
    LOOISRCVWrong2 = NULL
    
    LOOISR2CV = NULL
    LOOISR2CVWrong1 = NULL
    LOOISR2CVWrong2 = NULL
    
    LOOISPRCV = NULL
    LOOISPRCVWrong1 = NULL
    LOOISPRCVWrong2 = NULL
    
    LOOISPR2CV = NULL
    LOOISPR2CVWrong1 = NULL
    LOOISPR2CVWrong2 = NULL
    
    LOOVCCV = NULL
    LOOVCCVWrong1 = NULL
    LOOVCCVWrong2 = NULL
    
    LOOVCPCV = NULL
    LOOVCPCVWrong1 = NULL
    LOOVCPCVWrong2 = NULL
    
    LOOVCRCV = NULL
    LOOVCRCVWrong1 = NULL
    LOOVCRCVWrong2 = NULL
    
    LOOVCR2CV = NULL
    LOOVCR2CVWrong1 = NULL
    LOOVCR2CVWrong2 = NULL
    
    LOOVCPRCV = NULL
    LOOVCPRCVWrong1 = NULL
    LOOVCPRCVWrong2 = NULL
    
    LOOVCPR2CV = NULL
    LOOVCPR2CVWrong1 = NULL
    LOOVCPR2CVWrong2 = NULL
  }
  
  # 12. Calculate true model Scores ----
  condDistn = condMeanMVN(SigmaAA=rep(rGRFargsTruth$sigma^2, nrow(locs)), SigmaAB=SigmaGridToSample, SigmaBB=SigmaSample, 
                          ysB=ys, getFullCov=FALSE, getCondVar=TRUE)
  muAcondB = condDistn$muAcondB
  varAcondB = condDistn$varAcondB
  condDistn = condMeanMVN(SigmaAA=rep(rGRFargsWrong1$sigma^2, nrow(locs)), SigmaAB=SigmaGridToSampleWrong1, SigmaBB=SigmaSampleWrong1, 
                          ysB=ys, getFullCov=FALSE, getCondVar=TRUE)
  muAcondBWrong1 = condDistn$muAcondB
  varAcondBWrong1 = condDistn$varAcondB
  condDistn = condMeanMVN(SigmaAA=rep(rGRFargsWrong2$sigma^2, nrow(locs)), SigmaAB=SigmaGridToSampleWrong2, SigmaBB=SigmaSampleWrong2, 
                          ysB=ys, getFullCov=FALSE, getCondVar=TRUE)
  muAcondBWrong2 = condDistn$muAcondB
  varAcondBWrong2 = condDistn$varAcondB
  
  # calculate Scores for the grid cells (called MSE for historical reasons...)
  # trueMSE = mean((truth - muAcondB)^2) + sigmaEpsSq1
  # wrongMSE = mean((truth - muAcondBwrong)^2) + sigmaEpsSq1 # add sigmaEpsSq1, the true error variance, not sigmaEpsSq2
  
  trueMSE = expectedCRPS(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondB, est.var=varAcondB + sigmaEpsSqTrueLoc)
  wrongMSE1 = expectedCRPS(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong1, est.var=varAcondBWrong1 + sigmaEpsSqWrong1Loc)
  wrongMSE2 = expectedCRPS(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong2, est.var=varAcondBWrong2 + sigmaEpsSqWrong2Loc)
  
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
  estVarWVC = mean((1/rateEsts - 1)^2)
  
  # 13. Save result ----
  save(trueMSE, wrongMSE1, wrongMSE2, 
       LOOCVs, LOOCVsWrong1, LOOCVsWrong2, 
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
       LOOVCPCV, LOOVCPRCV, LOOVCPR2CV, 
       LOOVCCVWrong1, LOOVCRCVWrong1, LOOVCR2CVWrong1, 
       LOOVCPCVWrong1, LOOVCPRCVWrong1, LOOVCPR2CVWrong1, 
       LOOVCCVWrong2, LOOVCRCVWrong2, LOOVCR2CVWrong2, 
       LOOVCPCVWrong2, LOOVCPRCVWrong2, LOOVCPR2CVWrong2, 
       
       KFoldCVs, KFoldCVsWrong1, KFoldCVsWrong2, 
       KFoldCV, KFoldRCV, KFoldR2CV, 
       KFoldCVWrong1, KFoldRCVWrong1, KFoldR2CVWrong1, 
       KFoldCVWrong2, KFoldRCVWrong2, KFoldR2CVWrong2, 
       KFoldISCV, KFoldISRCV, KFoldISR2CV, 
       KFoldISPCV, KFoldISPRCV, KFoldISPR2CV, 
       KFoldISPCVWrong1, KFoldISPRCVWrong1, KFoldISPR2CVWrong1, 
       KFoldISCVWrong1, KFoldISRCVWrong1, KFoldISR2CVWrong1, 
       KFoldISPCVWrong2, KFoldISPRCVWrong2, KFoldISPR2CVWrong2, 
       KFoldISCVWrong2, KFoldISRCVWrong2, KFoldISR2CVWrong2, 
       KFoldVCCV, KFoldVCRCV, KFoldVCR2CV, 
       KFoldVCCVWrong1, KFoldVCRCVWrong1, KFoldVCR2CVWrong1, 
       KFoldVCCVWrong2, KFoldVCRCVWrong2, KFoldVCR2CVWrong2, 
       KFoldVCPCV, KFoldVCPRCV, KFoldVCPR2CV, 
       KFoldVCPCVWrong1, KFoldVCPRCVWrong1, KFoldVCPR2CVWrong1, 
       KFoldVCPCVWrong2, KFoldVCPRCVWrong2, KFoldVCPR2CVWrong2, 
       
       griddedCVs, griddedCVsWrong1, griddedCVsWrong2, gridNs, 
       iter, rGRFargsTruth, rGRFargsWrong1, rGRFargsWrong2, 
       n1, nx, ny, propMount, propSamplesMount, Ks, 
       sigmaEpsSqNonMount, sigmaEpsSqMount, 
       sigmaEpsSqNonMountWrong1, sigmaEpsSqMountWrong1, 
       sigmaEpsSqNonMountWrong2, sigmaEpsSqMountWrong2, 
       allSeeds, trueVarW, estVarWVC, 
       file=paste0("savedOutput/griddedCVtestNonstatErr/n1", n1, "_pMount", propMount, "_pSMount", propSamplesMount, 
                   "_s2M", sigmaEpsSqMount, "_", sigmaEpsSqMountWrong1, "_", sigmaEpsSqMountWrong2, 
                   "_s2NM", sigmaEpsSqNonMount, "_", sigmaEpsSqNonMountWrong1, "_", sigmaEpsSqNonMountWrong2, 
                   "_iter", iter, ".RData"))
  
  list(trueMSE=trueMSE, wrongMSE1=wrongMSE1, wrongMSE2=wrongMSE2, 
       LOOCVs=LOOCVs, LOOCVsWrong1=LOOCVsWrong1, LOOCVsWrong2=LOOCVsWrong2, 
       LOOCV=LOOCV, LOORCV=LOORCV, LOORCV=LOORCV, 
       LOOCVWrong1=LOOCVWrong1, LOOR2CVWrong1=LOOR2CVWrong1, LOOR2CVWrong1=LOOR2CVWrong1, 
       LOOCVWrong2=LOOCVWrong2, LOOR2CVWrong2=LOOR2CVWrong2, LOOR2CVWrong2=LOOR2CVWrong2, 
       LOOISCV=LOOISCV, LOOISRCV=LOOISRCV, LOOISR2CV=LOOISR2CV, 
       LOOISPCV=LOOISPCV, LOOISPRCV=LOOISPRCV, LOOISPR2CV=LOOISPR2CV, 
       LOOISPCVWrong1=LOOISPCVWrong1, LOOISPRCVWrong1=LOOISPRCVWrong1, LOOISPR2CVWrong1=LOOISPR2CVWrong1, 
       LOOISCVWrong1=LOOISCVWrong1, LOOISRCVWrong1=LOOISRCVWrong1, LOOISR2CVWrong1=LOOISR2CVWrong1, 
       LOOISPCVWrong2=LOOISPCVWrong2, LOOISPRCVWrong2=LOOISPRCVWrong2, LOOISPR2CVWrong2=LOOISPR2CVWrong2, 
       LOOISCVWrong2=LOOISCVWrong2, LOOISRCVWrong2=LOOISRCVWrong2, LOOISR2CVWrong2=LOOISR2CVWrong2, 
       LOOVCCV=LOOVCCV, LOOVCRCV=LOOVCRCV, LOOVCR2CV=LOOVCR2CV, 
       LOOVCCVWrong1=LOOVCCVWrong1, LOOVCRCVWrong1=LOOVCRCVWrong1, LOOVCR2CVWrong1=LOOVCR2CVWrong1, 
       LOOVCCVWrong2=LOOVCCVWrong2, LOOVCRCVWrong2=LOOVCRCVWrong2, LOOVCR2CVWrong2=LOOVCR2CVWrong2, 
       LOOVCPCV=LOOVCPCV, LOOVCPRCV=LOOVCPRCV, LOOVCPR2CV=LOOVCPR2CV, 
       LOOVCPCVWrong1=LOOVCPCVWrong1, LOOVCPRCVWrong1=LOOVCPRCVWrong1, LOOVCPR2CVWrong1=LOOVCPR2CVWrong1, 
       LOOVCPCVWrong2=LOOVCPCVWrong2, LOOVCPRCVWrong2=LOOVCPRCVWrong2, LOOVCPR2CVWrong2=LOOVCPR2CVWrong2, 
       
       KFoldCVs=KFoldCVs, KFoldCVsWrong1=KFoldCVsWrong1, KFoldCVsWrong2=KFoldCVsWrong2, 
       KFoldCV=KFoldCV, KFoldRCV=KFoldRCV, KFoldRCV=KFoldRCV, 
       KFoldCVWrong1=KFoldCVWrong1, KFoldR2CVWrong1=KFoldR2CVWrong1, KFoldR2CVWrong1=KFoldR2CVWrong1, 
       KFoldCVWrong2=KFoldCVWrong2, KFoldR2CVWrong2=KFoldR2CVWrong2, KFoldR2CVWrong2=KFoldR2CVWrong2, 
       KFoldISCV=KFoldISCV, KFoldISRCV=KFoldISRCV, KFoldISR2CV=KFoldISR2CV, 
       KFoldISPCV=KFoldISPCV, KFoldISPRCV=KFoldISPRCV, KFoldISPR2CV=KFoldISPR2CV, 
       KFoldISPCVWrong1=KFoldISPCVWrong1, KFoldISPRCVWrong1=KFoldISPRCVWrong1, KFoldISPR2CVWrong1=KFoldISPR2CVWrong1, 
       KFoldISCVWrong1=KFoldISCVWrong1, KFoldISRCVWrong1=KFoldISRCVWrong1, KFoldISR2CVWrong1=KFoldISR2CVWrong1, 
       KFoldISPCVWrong2=KFoldISPCVWrong2, KFoldISPRCVWrong2=KFoldISPRCVWrong2, KFoldISPR2CVWrong2=KFoldISPR2CVWrong2, 
       KFoldISCVWrong2=KFoldISCVWrong2, KFoldISRCVWrong2=KFoldISRCVWrong2, KFoldISR2CVWrong2=KFoldISR2CVWrong2, 
       KFoldVCCV=KFoldVCCV, KFoldVCRCV=KFoldVCRCV, KFoldVCR2CV=KFoldVCR2CV, 
       KFoldVCCVWrong1=KFoldVCCVWrong1, KFoldVCRCVWrong1=KFoldVCRCVWrong1, KFoldVCR2CVWrong1=KFoldVCR2CVWrong1, 
       KFoldVCCVWrong2=KFoldVCCVWrong2, KFoldVCRCVWrong2=KFoldVCRCVWrong2, KFoldVCR2CVWrong2=KFoldVCR2CVWrong2, 
       KFoldVCPCV=KFoldVCPCV, KFoldVCPRCV=KFoldVCPRCV, KFoldVCPR2CV=KFoldVCPR2CV, 
       KFoldVCPCVWrong1=KFoldVCPCVWrong1, KFoldVCPRCVWrong1=KFoldVCPRCVWrong1, KFoldVCPR2CVWrong1=KFoldVCPR2CVWrong1, 
       KFoldVCPCVWrong2=KFoldVCPCVWrong2, KFoldVCPRCVWrong2=KFoldVCPRCVWrong2, KFoldVCPR2CVWrong2=KFoldVCPR2CVWrong2, 
       
       griddedCVs=griddedCVs, griddedCVsWrong1=griddedCVsWrong1, griddedCVsWrong2=griddedCVsWrong2, 
       gridNs=gridNs, iter=iter, rGRFargsTruth=rGRFargsTruth, 
       rGRFargsWrong1=rGRFargsWrong1, rGRFargsWrong2=rGRFargsWrong2, 
       n1=n1, nx=nx, ny=ny, Ks=Ks, 
       propMount=propMount, propSamplesMount=propSamplesMount, 
       sigmaEpsSqNonMount=sigmaEpsSqNonMount, 
       sigmaEpsSqMount=sigmaEpsSqMount, 
       sigmaEpsSqNonMountWrong1=sigmaEpsSqNonMountWrong1, 
       sigmaEpsSqMountWrong1=sigmaEpsSqMountWrong1, 
       sigmaEpsSqNonMountWrong2=sigmaEpsSqNonMountWrong2, 
       sigmaEpsSqMountWrong2=sigmaEpsSqMountWrong2, 
       allSeeds=allSeeds, trueVarW=trueVarW, estVarWVC=estVarWVC)
}

getNonstatErrECRPSpropGrid = function(n=50, nsim=100, 
                                      propMountSeq=seq(0.05, .95, by=.05), propSamplesMountSeq=propMountSeq, 
                                      rGRFargsTruth=NULL, rGRFargsMount=NULL, 
                                      rGRFargsWrong1=NULL, rGRFargsWrong2=NULL, 
                                      nx=100, ny=100, sigmaEpsSqNonMount=.1^2, sigmaEpsSqMount=1^2, 
                                      sigmaEpsSqNonMountWrong1=.1^2, sigmaEpsSqMountWrong1=.1^2, 
                                      sigmaEpsSqNonMountWrong2=.15^2, sigmaEpsSqMountWrong2=.9^2, 
                                      seed=123, printProgress=TRUE) {
  
  set.seed(seed)
  
  trueECRPSMat = matrix(nrow=length(propMountSeq), ncol=length(propSamplesMountSeq))
  wrongECRPS1Mat = matrix(nrow=length(propMountSeq), ncol=length(propSamplesMountSeq))
  wrongECRPS2Mat = matrix(nrow=length(propMountSeq), ncol=length(propSamplesMountSeq))
  trueECRPSSampleWeightedMat = matrix(nrow=length(propMountSeq), ncol=length(propSamplesMountSeq))
  wrongECRPS1SampleWeightedMat = matrix(nrow=length(propMountSeq), ncol=length(propSamplesMountSeq))
  wrongECRPS2SampleWeightedMat = matrix(nrow=length(propMountSeq), ncol=length(propSamplesMountSeq))
  
  startTime = proc.time()[3]
  totalIter = length(propMountSeq) * length(propSamplesMountSeq)
  for(i in 1:length(propMountSeq)) {
    propMount = propMountSeq[i]
    
    for(j in 1:length(propSamplesMountSeq)) {
      propSamplesMount = propSamplesMountSeq[j]
      
      out = getNonstatErrECRPS(n=n, nsim=nsim, 
                         propMount=propMount, propSamplesMount=propSamplesMount, 
                         rGRFargsTruth=rGRFargsTruth, rGRFargsMount=rGRFargsMount, 
                         rGRFargsWrong1=rGRFargsWrong1, rGRFargsWrong2=rGRFargsWrong2, 
                         nx=nx, ny=ny, 
                         sigmaEpsSqNonMount=sigmaEpsSqNonMount, sigmaEpsSqMount=sigmaEpsSqMount, 
                         sigmaEpsSqNonMountWrong1=sigmaEpsSqNonMountWrong1,sigmaEpsSqMountWrong1=sigmaEpsSqMountWrong1, 
                         sigmaEpsSqNonMountWrong2=sigmaEpsSqNonMountWrong2,sigmaEpsSqMountWrong2=sigmaEpsSqMountWrong2, 
                         seed=NULL, printProgress=FALSE)
      
      trueECRPSMat[i,j] = out$trueECRPS
      wrongECRPS1Mat[i,j] = out$wrongECRPS1
      wrongECRPS2Mat[i,j] = out$wrongECRPS2
      trueECRPSSampleWeightedMat[i,j] = out$trueECRPSSampleWeighted
      wrongECRPS1SampleWeightedMat[i,j] = out$wrongECRPS1SampleWeighted
      wrongECRPS2SampleWeightedMat[i,j] = out$wrongECRPS2SampleWeighted
      
      if(printProgress) {
        currTime = proc.time()[3]
        tleft=estTimeLeft(startTime, currTime, finishedIter=(i-1)*length(propSamplesMountSeq) + j, totalIter)
        print(paste0("finished iter ", i, "/", length(propMountSeq), ", ", j, "/", length(propSamplesMountSeq), ". Est. time left: ", round(tleft/60, 1), " minutes..."))
      }
    }
  }
  
  # plot results ----
  browser()
  
  # return results ----
  out = list(trueECRPSMat=trueECRPSMat, 
             wrongECRPS1Mat=wrongECRPS1Mat, 
             wrongECRPS2Mat=wrongECRPS2Mat, 
             trueECRPSSampleWeightedMat=trueECRPSSampleWeightedMat, 
             wrongECRPS1SampleWeightedMat=wrongECRPS1SampleWeightedMat, 
             wrongECRPS2SampleWeightedMat=wrongECRPS2SampleWeightedMat)
}

# calculate expected CRPS for each of the models
getNonstatErrECRPS = function(n=50, nsim=100, 
                              propMount=.3, propSamplesMount=propMount/5, 
                              rGRFargsTruth=NULL, rGRFargsMount=NULL, 
                              rGRFargsWrong1=NULL, rGRFargsWrong2=NULL, 
                              nx=100, ny=100, sigmaEpsSqNonMount=.1^2, sigmaEpsSqMount=1^2, 
                              sigmaEpsSqNonMountWrong1=.1^2, sigmaEpsSqMountWrong1=.1^2, 
                              sigmaEpsSqNonMountWrong2=.15^2, sigmaEpsSqMountWrong2=.9^2, 
                              seed=NULL, printProgress=TRUE) {
  if(!is.null(seed)) {
    set.seed(123)
  }
  
  # simulate samples, calculate true
  trueMSEs = numeric(nsim)
  wrongMSE1s = numeric(nsim)
  wrongMSE2s = numeric(nsim)
  
  trueMSEsSampleWeighted = numeric(nsim)
  wrongMSE1sSampleWeighted = numeric(nsim)
  wrongMSE2sSampleWeighted = numeric(nsim)
  
  startTime = proc.time()[3]
  for(i in 1:nsim) {
    
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
    
    # 2.5. Simulate 1 GRF for mountains ----
    #    on unit square
    mountGRF = do.call("rGRF", rGRFargsMount)
    mountCov = mountGRF$truth
    isMount = mountCov > quantile(mountCov, .7, type=1) # exactly 30% of domain is mountainous
    
    # 3. Simulate sample distribution ----
    #    on mountainous and non-mountainous part of domain
    sampleRates = rep(1, nrow(locs))
    nMount = round(n*propSamplesMount)
    nNonMount = n - nMount
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
    sampleI = c(sampleIMount, sampleINonMount)
    mounts = c(rep(TRUE, nMount), rep(FALSE, nNonMount))
    sigmaEpsSqWrong1 = c(rep(sigmaEpsSqMountWrong1, nMount), rep(sigmaEpsSqNonMountWrong1, nNonMount))
    sigmaEpsSqWrong2 = c(rep(sigmaEpsSqMountWrong2, nMount), rep(sigmaEpsSqNonMountWrong2, nNonMount))
    sigmaEpsSqTrue = c(rep(sigmaEpsSqMount, nMount), rep(sigmaEpsSqNonMount, nNonMount))
    
    sigmaEpsSqTrueLoc = rep(sigmaEpsSqNonMount, nrow(locs))
    sigmaEpsSqTrueLoc[isMount] = sigmaEpsSqMount
    sigmaEpsSqWrong1Loc = rep(sigmaEpsSqNonMountWrong1, nrow(locs))
    sigmaEpsSqWrong1Loc[isMount] = sigmaEpsSqMountWrong1
    sigmaEpsSqWrong2Loc = rep(sigmaEpsSqNonMountWrong2, nrow(locs))
    sigmaEpsSqWrong2Loc[isMount] = sigmaEpsSqMountWrong2
    
    # 5. Get covariance matrices ----
    #    for sample and crossCov to all locs
    distMatXs = rdist(xs, xs)
    SigmaSample = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsTruth$cov.args$range, 
                                 smoothness=rGRFargsTruth$cov.args$smoothness, distMat=distMatXs) * rGRFargsTruth$sigma^2
    SigmaSampleWrong1 = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsWrong1$cov.args$range, 
                                       smoothness=rGRFargsWrong1$cov.args$smoothness, distMat=distMatXs) * rGRFargsWrong1$sigma^2
    SigmaSampleWrong2 = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsWrong2$cov.args$range, 
                                       smoothness=rGRFargsWrong2$cov.args$smoothness, distMat=distMatXs) * rGRFargsWrong2$sigma^2
    SigmaSample = SigmaSample + diag(sigmaEpsSqTrue)
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
    
    # 12. Calculate true model Scores ----
    condDistn = condMeanMVN(SigmaAA=rep(rGRFargsTruth$sigma^2, nrow(locs)), SigmaAB=SigmaGridToSample, SigmaBB=SigmaSample, 
                            ysB=ys, getFullCov=FALSE, getCondVar=TRUE)
    muAcondB = condDistn$muAcondB
    varAcondB = condDistn$varAcondB
    condDistn = condMeanMVN(SigmaAA=rep(rGRFargsWrong1$sigma^2, nrow(locs)), SigmaAB=SigmaGridToSampleWrong1, SigmaBB=SigmaSampleWrong1, 
                            ysB=ys, getFullCov=FALSE, getCondVar=TRUE)
    muAcondBWrong1 = condDistn$muAcondB
    varAcondBWrong1 = condDistn$varAcondB
    condDistn = condMeanMVN(SigmaAA=rep(rGRFargsWrong2$sigma^2, nrow(locs)), SigmaAB=SigmaGridToSampleWrong2, SigmaBB=SigmaSampleWrong2, 
                            ysB=ys, getFullCov=FALSE, getCondVar=TRUE)
    muAcondBWrong2 = condDistn$muAcondB
    varAcondBWrong2 = condDistn$varAcondB
    
    # calculate Scores for the grid cells (called MSE for historical reasons...)
    # trueMSE = mean((truth - muAcondB)^2) + sigmaEpsSq1
    # wrongMSE = mean((truth - muAcondBwrong)^2) + sigmaEpsSq1 # add sigmaEpsSq1, the true error variance, not sigmaEpsSq2
    
    trueMSEs[i] = expectedCRPS(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondB, est.var=varAcondB + sigmaEpsSqTrueLoc)
    wrongMSE1s[i] = expectedCRPS(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong1, est.var=varAcondBWrong1 + sigmaEpsSqWrong1Loc)
    wrongMSE2s[i] = expectedCRPS(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong2, est.var=varAcondBWrong2 + sigmaEpsSqWrong2Loc)
    
    trueMSEsSampleWeighted[i] = expectedCRPS(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondB, est.var=varAcondB + sigmaEpsSqTrueLoc, 
                                             weights=sampleProbs)
    wrongMSE1sSampleWeighted[i] = expectedCRPS(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong1, est.var=varAcondBWrong1 + sigmaEpsSqWrong1Loc, 
                                               weights=sampleProbs)
    wrongMSE2sSampleWeighted[i] = expectedCRPS(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong2, est.var=varAcondBWrong2 + sigmaEpsSqWrong2Loc, 
                                               weights=sampleProbs)
    
    if(printProgress) {
      currTime = proc.time()[3]
      tleft=estTimeLeft(startTime, currTime, i, nsim)
      print(paste0("finished iter ", i, "/", nsim, ". Est. time left: ", round(tleft/60, 1), " minutes..."))
    }
  }
  
  trueECRPS = mean(trueMSEs)
  wrongECRPS1 = mean(wrongMSE1s)
  wrongECRPS2 = mean(wrongMSE2s)
  
  trueECRPSSampleWeighted = mean(trueMSEsSampleWeighted)
  wrongECRPS1SampleWeighted = mean(wrongMSE1sSampleWeighted)
  wrongECRPS2SampleWeighted = mean(wrongMSE2sSampleWeighted)
  
  list(trueMSEs=trueMSEs, wrongMSE1s=wrongMSE1s, wrongMSE2s=wrongMSE2s, 
       trueMSEsSampleWeighted=trueMSEsSampleWeighted, 
       wrongMSE1sSampleWeighted=wrongMSE1sSampleWeighted, 
       wrongMSE2sSampleWeighted=wrongMSE2sSampleWeighted, 
       trueECRPS=trueECRPS, wrongECRPS1=wrongECRPS1, wrongECRPS2=wrongECRPS2, 
       trueECRPSSampleWeighted=trueECRPSSampleWeighted, 
       wrongECRPS1SampleWeighted=wrongECRPS1SampleWeighted, 
       wrongECRPS2SampleWeighted=wrongECRPS2SampleWeighted)
}






