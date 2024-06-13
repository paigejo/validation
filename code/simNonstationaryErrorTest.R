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
                                          propMount=.3, oversampleMountRatio=1/5, 
                                          n1=50, gridNs=2^(1:6), Ks=c(9, 25), iter=1, 
                                          rGRFargsWrong1=NULL, rGRFargsWrong2=NULL, 
                                          nx=500, ny=nx, sigmaEpsSqNonMount=.1^2, sigmaEpsSqMount=1^2, 
                                          sigmaEpsSqNonMountWrong1=.1^2, sigmaEpsSqMountWrong1=.1^2, 
                                          sigmaEpsSqNonMountWrong2=.1, sigmaEpsSqMountWrong2=.1, 
                                          allSeeds=123, printProgress=FALSE, printPhat=FALSE, 
                                          saveResults=TRUE, subsample=nx/100) {
  
  if(iter > 1) {
    currT = proc.time()[3]
    tLeft = estTimeLeft(startT, currT, finishedIter=iter-1, totalIter=totalIter)/60 # convert to minutes
    print(paste0("iteration ", iter, "/", length(allSeeds), ". Est minutes left: ", round(tLeft, digits=1)))
  } else {
    print(paste0("iteration ", iter, "/", length(allSeeds)))
  }
  
  tmp = oversampleMountRatio * propMount/(1-propMount)
  # propSamplesMount = propMount * propSamplesMountRatio
  propSamplesMount = tmp/(1+tmp)
  
  doLOO = TRUE
  if(n1 > 50) {
    gridNs = c(3, 5)
    
    if(n1 > 500) {
      doLOO = FALSE
      if(n1 > 2000) {
        Ks = Ks[Ks <= 10]
        gridNs = gridNs[gridNs <= 4]
      }
    }
  }
  
  if(!is.null(allSeeds)) {
    set.seed(allSeeds[iter])
  }
  
  # set default GRF parameters
  if(is.null(rGRFargsTruth)) {
    rGRFargsTruth = list(mu=0, sigma=1, 
                         cov.args=list(Covariance="Matern", range=0.1, smoothness=1.5), 
                         delta=10, sigmaEpsSq=0, nx=nx, ny=ny)
  }
  if(is.null(rGRFargsWrong1)) {
    rGRFargsWrong1 = rGRFargsTruth
    # rGRFargsWrong1$cov.args$range = 1.25 * rGRFargsTruth$cov.args$range
    # rGRFargsWrong1$cov.args$smoothness = 1.5
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
  
  if(!is.null(subsample)) {
    truthGRF = downsampleGridI(truthGRF, subsample)$newSimGRF
    locsSub = truthGRF$locs
    truthSub = truthGRF$truth
  }
  
  # 2.5. Simulate 1 GRF for mountains ----
  #    on unit square
  mountGRF = do.call("rGRF", rGRFargsMount)
  mountCov = mountGRF$truth
  isMount = mountCov > quantile(mountCov, .7, type=1) # exactly 30% of domain is mountainous
  
  if(!is.null(subsample)) {
    mountGRF = downsampleGridI(mountGRF, subsample)$newSimGRF
    mountCov = mountGRF$truth
    isMountSub = mountCov > quantile(mountCov, .7, type=1)
  }
  
  
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
  
  if(!is.null(subsample)) {
    # now that we've simulated the data, we don't need the super high resolution 
    # GRF samples anymore
    locs = locsSub
    truth = truthSub
    mountCov = mountCovSub
  }
  
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
        # blockCVs[j] = crps(truth=testYs, est=muAcondB, est.var=varAcondB)
        # blockCVsWrong1[j] = crps(truth=testYs, est=muAcondBWrong1, est.var=varAcondBWrong1)
        # blockCVsWrong2[j] = crps(truth=testYs, est=muAcondBWrong2, est.var=varAcondBWrong2)
        blockCVs[j] = intervalScore(truth=testYs, est=muAcondB, var=varAcondB)
        blockCVsWrong1[j] = intervalScore(truth=testYs, est=muAcondBWrong1, var=varAcondBWrong1)
        blockCVsWrong2[j] = intervalScore(truth=testYs, est=muAcondBWrong2, var=varAcondBWrong2)
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
      # KFoldCVs[isTest,i] = crps(truth=testYs, est=muAcondB, est.var=varAcondB, getAverage = FALSE)
      # KFoldCVsWrong1[isTest,i] = crps(truth=testYs, est=muAcondBWrong1, est.var=varAcondBWrong1, getAverage = FALSE)
      # KFoldCVsWrong2[isTest,i] = crps(truth=testYs, est=muAcondBWrong2, est.var=varAcondBWrong2, getAverage = FALSE)
      KFoldCVs[isTest,i] = intervalScore(truth=testYs, est=muAcondB, var=varAcondB, getAverage = FALSE)
      KFoldCVsWrong1[isTest,i] = intervalScore(truth=testYs, est=muAcondBWrong1, var=varAcondBWrong1, getAverage = FALSE)
      KFoldCVsWrong2[isTest,i] = intervalScore(truth=testYs, est=muAcondBWrong2, var=varAcondBWrong2, getAverage = FALSE)
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
      # LOOCVs[i] = crps(truth=testYs, est=muAcondB, est.var=varAcondB)
      # LOOCVsWrong1[i] = crps(truth=testYs, est=muAcondBWrong1, est.var=varAcondBWrong1)
      # LOOCVsWrong2[i] = crps(truth=testYs, est=muAcondBWrong2, est.var=varAcondBWrong2)
      LOOCVs[i] = intervalScore(truth=testYs, est=muAcondB, var=varAcondB)
      LOOCVsWrong1[i] = intervalScore(truth=testYs, est=muAcondBWrong1, var=varAcondBWrong1)
      LOOCVsWrong2[i] = intervalScore(truth=testYs, est=muAcondBWrong2, var=varAcondBWrong2)
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
    LOOCVs = numeric(n1)
    LOOCVsWrong1 = numeric(n1)
    LOOCVsWrong2 = numeric(n1)
    
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
    LOOISCVWrong1 = NULL
    LOOISCVWrong2 = NULL
    
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
  
  trueMSE = expectedIntervalScore(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondB, est.var=varAcondB + sigmaEpsSqTrueLoc)
  wrongMSE1 = expectedIntervalScore(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong1, est.var=varAcondBWrong1 + sigmaEpsSqWrong1Loc)
  wrongMSE2 = expectedIntervalScore(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong2, est.var=varAcondBWrong2 + sigmaEpsSqWrong2Loc)
  
  # trueCRPS = expectedCRPS(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondB, est.var=varAcondB + sigmaEpsSqTrueLoc)
  # wrongCRPS1 = expectedCRPS(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong1, est.var=varAcondBWrong1 + sigmaEpsSqWrong1Loc)
  # wrongCRPS2 = expectedCRPS(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong2, est.var=varAcondBWrong2 + sigmaEpsSqWrong2Loc)
  
  if(FALSE) {
    
    pdf("figures/nonstatErrorTest/test_Truth.pdf", width=5, height=5)
    plotWithColor(locs[,1], locs[,2], truth, cex=.2, pch=19)
    dev.off()
    pdf("figures/nonstatErrorTest/test_Mount.pdf", width=5, height=5)
    plotWithColor(locs[,1], locs[,2], isMount, cex=.2, pch=19, colScale=c("green", "brown"))
    dev.off()
    pdf("figures/nonstatErrorTest/test_dat.pdf", width=5, height=5)
    plotWithColor(xs[,1], xs[,2], ys, cex=.1, pch=19)
    dev.off()
    pdf("figures/nonstatErrorTest/test_TruePreds.pdf", width=5, height=5)
    plotWithColor(locs[,1], locs[,2], muAcondB, cex=.2, pch=19)
    dev.off()
    pdf("figures/nonstatErrorTest/test_TrueResids.pdf", width=5, height=5)
    plotWithColor(locs[,1], locs[,2], truth-muAcondB, cex=.2, pch=19)
    dev.off()
    diffs = c(0,diff(truth))
    cols = makeRedBlueDivergingColors(64, range(diffs), center=0)
    pdf("figures/nonstatErrorTest/test_TrueDiffs.pdf", width=5, height=5)
    plotWithColor(locs[,1], locs[,2], diffs, cex=.2, pch=19, colScale=cols)
    dev.off()
    pdf("figures/nonstatErrorTest/test_Wrong1Preds.pdf", width=5, height=5)
    plotWithColor(locs[,1], locs[,2], muAcondBWrong1, cex=.2, pch=19)
    dev.off()
    pdf("figures/nonstatErrorTest/test_Wrong1Resids.pdf", width=5, height=5)
    plotWithColor(locs[,1], locs[,2], truth-muAcondBWrong1, cex=.2, pch=19)
    dev.off()
    pdf("figures/nonstatErrorTest/test_Wrong2Preds.pdf", width=5, height=5)
    plotWithColor(locs[,1], locs[,2], muAcondBWrong2, cex=.2, pch=19)
    dev.off()
    pdf("figures/nonstatErrorTest/test_Wrong2Resids.pdf", width=5, height=5)
    plotWithColor(locs[,1], locs[,2], truth-muAcondBWrong2, cex=.2, pch=19)
    dev.off()
  }
  
  # calculate some theoretical properties of the sampling weights
  # mean(1/sampleGRF1$truth)
  # mean(1/sampleGRF2$truth)
  trueVarW = -1 + mean(1/sampleProbs)
  estVarWVC = mean((1/rateEsts - 1)^2)
  
  # 13. Save result ----
  if(saveResults) {
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
         file=paste0("savedOutput/griddedCVtestNonstatErr/n1", n1, "_pMount", propMount, "_Rosamp", oversampleMountRatio, 
                     "_s2M", sigmaEpsSqMount, "_", sigmaEpsSqMountWrong1, "_", sigmaEpsSqMountWrong2, 
                     "_s2NM", sigmaEpsSqNonMount, "_", sigmaEpsSqNonMountWrong1, "_", sigmaEpsSqNonMountWrong2, 
                     "_iter", iter, ".RData"))
  }
  
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
                              sigmaEpsSqNonMountWrong1=.2^2, sigmaEpsSqMountWrong1=.2^2, 
                              sigmaEpsSqNonMountWrong2=.3^2, sigmaEpsSqMountWrong2=.8^2, 
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
    
    if(subsample == 1) {
      distLocsToXs = rdist(locs, xs)
      distLocsToSample = distLocsToXs
      SigmaGridToSample = stationary.cov(locs, xs, Covariance="Matern", aRange=rGRFargsTruth$cov.args$range, 
                                         smoothness=rGRFargsTruth$cov.args$smoothness, distMat=distLocsToSample) * rGRFargsTruth$sigma^2
      SigmaGridToSampleWrong1 = stationary.cov(locs, xs, Covariance="Matern", aRange=rGRFargsWrong1$cov.args$range, 
                                               smoothness=rGRFargsWrong1$cov.args$smoothness, distMat=distLocsToSample) * rGRFargsWrong1$sigma^2
      SigmaGridToSampleWrong2 = stationary.cov(locs, xs, Covariance="Matern", aRange=rGRFargsWrong2$cov.args$range, 
                                               smoothness=rGRFargsWrong2$cov.args$smoothness, distMat=distLocsToSample) * rGRFargsWrong2$sigma^2
    } else {
      distLocsToXs = rdist(locsSub, xs)
      distLocsToSample = distLocsToXs
      SigmaGridToSample = stationary.cov(locsSub, xs, Covariance="Matern", aRange=rGRFargsTruth$cov.args$range, 
                                         smoothness=rGRFargsTruth$cov.args$smoothness, distMat=distLocsToSample) * rGRFargsTruth$sigma^2
      SigmaGridToSampleWrong1 = stationary.cov(locsSub, xs, Covariance="Matern", aRange=rGRFargsWrong1$cov.args$range, 
                                               smoothness=rGRFargsWrong1$cov.args$smoothness, distMat=distLocsToSample) * rGRFargsWrong1$sigma^2
      SigmaGridToSampleWrong2 = stationary.cov(locsSub, xs, Covariance="Matern", aRange=rGRFargsWrong2$cov.args$range, 
                                               smoothness=rGRFargsWrong2$cov.args$smoothness, distMat=distLocsToSample) * rGRFargsWrong2$sigma^2
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
    
    trueMSEs[i] = expectedIntervalScore(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondB, est.var=varAcondB + sigmaEpsSqTrueLoc)
    wrongMSE1s[i] = expectedIntervalScore(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong1, est.var=varAcondBWrong1 + sigmaEpsSqWrong1Loc)
    wrongMSE2s[i] = expectedIntervalScore(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong2, est.var=varAcondBWrong2 + sigmaEpsSqWrong2Loc)
    
    trueMSEsSampleWeighted[i] = expectedIntervalScore(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondB, est.var=varAcondB + sigmaEpsSqTrueLoc, 
                                             weights=sampleProbs)
    wrongMSE1sSampleWeighted[i] = expectedIntervalScore(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong1, est.var=varAcondBWrong1 + sigmaEpsSqWrong1Loc, 
                                               weights=sampleProbs)
    wrongMSE2sSampleWeighted[i] = expectedIntervalScore(truth=truth, truth.var=sigmaEpsSqTrueLoc, est=muAcondBWrong2, est.var=varAcondBWrong2 + sigmaEpsSqWrong2Loc, 
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

getNonstatErrECRPSpropGrid = function(n=50, nsim=100, 
                                      propMountSeq=seq(0.05, .95, by=.05), propSamplesMountSeq=propMountSeq, 
                                      rGRFargsTruth=NULL, rGRFargsMount=NULL, 
                                      rGRFargsWrong1=NULL, rGRFargsWrong2=NULL, 
                                      nx=100, ny=100, sigmaEpsSqNonMount=.1^2, sigmaEpsSqMount=1^2, 
                                      sigmaEpsSqNonMountWrong1=.1^2, sigmaEpsSqMountWrong1=.1^2, 
                                      sigmaEpsSqNonMountWrong2=.15^2, sigmaEpsSqMountWrong2=.9^2, 
                                      seed=123, printProgress=TRUE) {
  
  
  
}

getTransitionErrVar = function(propMount=.3, 
                               oversampleMountRatioSeq=c(1/10, 1/5, 1/2, 1, 2, 5, 10), 
                               sigmaEpsSqNonMount=.1^2, sigmaEpsSqMount=1^2, 
                               score=c("is", "crps"), ns=c(50, 500, 10000), 
                               regenResults=TRUE, saveResults=NULL) {
  
  score = match.arg(score)
  
  
  startTime = proc.time()[3]
  for(i in 1:length(oversampleMountRatioSeq)) {
    print(paste0("i=", i, "/", length(oversampleMountRatioSeq)))
    oversampleMountRatio = oversampleMountRatioSeq[i]
    tmp = oversampleMountRatio * propMount/(1-propMount)
    # propSamplesMount = propMount * propSamplesMountRatio
    propSamplesMount = tmp/(1+tmp)
    
    # get the scoring rule requested
    if(score == "crps") {
      scoreName = "CRPS"
      scoreFunDom = function(sigmaEpsSqNonMountWrong, sigmaEpsSqMountWrong=sigmaEpsSqNonMountWrong) {
        expectedCRPS(0, truth.var = sigmaEpsSqMount, est = 0, est.var = sigmaEpsSqMountWrong, getAverage = FALSE)*propMount + expectedCRPS(0, truth.var = sigmaEpsSqNonMount, est = 0, est.var = sigmaEpsSqNonMountWrong, getAverage = FALSE)*(1-propMount)
      }
      scoreFunDat = function(sigmaEpsSqNonMountWrong, sigmaEpsSqMountWrong=sigmaEpsSqNonMountWrong) {
        expectedCRPS(0, truth.var = sigmaEpsSqMount, est = 0, est.var = sigmaEpsSqMountWrong, getAverage = FALSE)*propSamplesMount + expectedCRPS(0, truth.var = sigmaEpsSqNonMount, est = 0, est.var = sigmaEpsSqNonMountWrong, getAverage = FALSE)*(1-propSamplesMount)
      }
      scoreVarDeltaM = function(sigmaEpsSqMountWrong1, sigmaEpsSqMountWrong2=sigmaEpsSqMountWrong1) {
        stop("not yet implemented")
      }
      scoreVarDeltaP = function(sigmaEpsSqNonMountWrong1, sigmaEpsSqNonMountWrong2=sigmaEpsSqNonMountWrong1) {
        stop("not yet implemented")
      }
    } else {
      scoreName = "IS"
      scoreFunDom = function(sigmaEpsSqNonMountWrong, sigmaEpsSqMountWrong=sigmaEpsSqNonMountWrong) {
        expectedIntervalScore(0, truth.var = sigmaEpsSqMount, est = 0, est.var = sigmaEpsSqMountWrong, getAverage = FALSE)*propMount + expectedIntervalScore(0, truth.var = sigmaEpsSqNonMount, est = 0, est.var = sigmaEpsSqNonMountWrong, getAverage = FALSE)*(1-propMount)
      }
      scoreFunDat = function(sigmaEpsSqNonMountWrong, sigmaEpsSqMountWrong=sigmaEpsSqNonMountWrong) {
        expectedIntervalScore(0, truth.var = sigmaEpsSqMount, est = 0, est.var = sigmaEpsSqMountWrong, getAverage = FALSE)*propSamplesMount + expectedIntervalScore(0, truth.var = sigmaEpsSqNonMount, est = 0, est.var = sigmaEpsSqNonMountWrong, getAverage = FALSE)*(1-propSamplesMount)
      }
      scoreVarDeltaM = function(sigmaEpsSqMountWrong1, sigmaEpsSqMountWrong2=sigmaEpsSqMountWrong1) {
        sigma1SqM = varIntervalScore(0, truth.var = sigmaEpsSqMount, est = 0, est.var = sigmaEpsSqMountWrong1, getAverage = FALSE)
        sigma2SqM = varIntervalScore(0, truth.var = sigmaEpsSqMount, est = 0, est.var = sigmaEpsSqMountWrong2, getAverage = FALSE)
        rho12M = covIntervalScore(0, truth.var = sigmaEpsSqMount, est1 = 0, est.var1 = sigmaEpsSqMountWrong1, 
                                  est2 = 0, est.var2 = sigmaEpsSqMountWrong2, getAverage = FALSE)
        (1/propSamplesMount) * (sigma1SqM + sigma2SqM - 2*rho12M)
      }
      scoreVarDeltaP = function(sigmaEpsSqNonMountWrong1, sigmaEpsSqNonMountWrong2=sigmaEpsSqNonMountWrong1) {
        sigma1SqP = varIntervalScore(0, truth.var = sigmaEpsSqNonMount, est = 0, est.var = sigmaEpsSqNonMountWrong1, getAverage = FALSE)
        sigma2SqP = varIntervalScore(0, truth.var = sigmaEpsSqNonMount, est = 0, est.var = sigmaEpsSqNonMountWrong2, getAverage = FALSE)
        rho12P = covIntervalScore(0, truth.var = sigmaEpsSqNonMount, est1 = 0, est.var1 = sigmaEpsSqNonMountWrong1, 
                                  est2 = 0, est.var2 = sigmaEpsSqNonMountWrong2, getAverage = FALSE)
        (1/(1-propSamplesMount)) * (sigma1SqP + sigma2SqP - 2*rho12P)
      }
    }
    
    correctSelProb = function(n, sigmaEpsSqNonMountWrong1, sigmaEpsSqMountWrong1=sigmaEpsSqNonMountWrong1, 
                              sigmaEpsSqNonMountWrong2, sigmaEpsSqMountWrong2=sigmaEpsSqNonMountWrong2) {
      # calculate MVN mean, variance of joint distribution between LOO score diff and IW score diff
      mu1M = expectedIntervalScore(0, truth.var = sigmaEpsSqMount, est = 0, est.var = sigmaEpsSqMountWrong1, getAverage = FALSE)
      mu2M = expectedIntervalScore(0, truth.var = sigmaEpsSqMount, est = 0, est.var = sigmaEpsSqMountWrong2, getAverage = FALSE)
      mu1P = expectedIntervalScore(0, truth.var = sigmaEpsSqNonMount, est = 0, est.var = sigmaEpsSqNonMountWrong1, getAverage = FALSE)
      mu2P = expectedIntervalScore(0, truth.var = sigmaEpsSqNonMount, est = 0, est.var = sigmaEpsSqNonMountWrong2, getAverage = FALSE)
      muM = mu1M - mu2M
      muP = mu1P - mu2P
      sigmaSqM = scoreVarDeltaM(sigmaEpsSqMountWrong1, sigmaEpsSqMountWrong2)
      sigmaSqP = scoreVarDeltaP(sigmaEpsSqNonMountWrong1, sigmaEpsSqNonMountWrong2)
      
      # first find which model is better
      model1Better = (propMount * muM + (1-propMount) * muP) < 0
      
      # now get probability LOOCV determines the correct model
      meanLOO = (propSamplesMount * muM + (1-propSamplesMount) * muP)
      sigmaSqLOO = (propSamplesMount^2 * sigmaSqM + (1-propSamplesMount)^2 * sigmaSqP)*(1/n)
      selProbsLOO = numeric(length(meanLOO))
      selProbsLOO = pnorm(0, mean=meanLOO, sd=sqrt(sigmaSqLOO))
      selProbsLOO[!model1Better] = 1-selProbsLOO[!model1Better]
      
      # get the same for IWCV
      meanIW = (propMount * muM + (1-propMount) * muP)
      sigmaSqIW = (propMount^2 * sigmaSqM + (1-propMount)^2 * sigmaSqP)*(1/n)
      selProbsIW = numeric(length(meanIW))
      selProbsIW = pnorm(0, mean=meanIW, sd=sqrt(sigmaSqIW))
      selProbsIW[!model1Better] = 1-selProbsIW[!model1Better]
      
      # A = rbind(c(propMount, 1-propMount),
      #           c(propSamplesMount, 1-propSamplesMount))
      # 
      # startTime = proc.time()[3]
      # thisProb = function(i) {
      #   thisSigmaSqM = sigmaSqM[i]
      #   thisSigmaSqP = sigmaSqP[i]
      #   thisMuM = muM[i]
      #   thisMuP = muP[i]
      # 
      #   D = (1/n) * diag(c(thisSigmaSqM, thisSigmaSqP))
      #   Sigma = A %*% D %*% t(A)
      #   mu = A %*% c(thisMuM, thisMuP)
      # 
      #   pneg = pmvnorm(lower=-Inf,
      #                  upper=c(0, 0),
      #                  mean=c(mu),
      #                  sigma=Sigma)[1]
      #   ppos = pmvnorm(lower=0,
      #                  upper=Inf,
      #                  mean=c(mu),
      #                  sigma=Sigma)[1]
      # 
      #   tenthN = floor(length(sigmaSqM)/10)
      #   if(((i %% tenthN) == 0) || i == 1000) {
      #     currTime = proc.time()[3]
      #     tLeft = estTimeLeft(startTime, currTime, i, length(sigmaSqM))
      #     print(paste0("iter ", i, "/", length(sigmaSqM), ", est time left: ", tLeft/60, " minutes"))
      #   }
      #   
      #   pneg + ppos
      # }
      # 
      # selProbs = sapply(1:length(mu1M), thisProb)
      # selProbs
      
      cbind(selProbLOO=selProbsLOO, selProbIW=selProbsIW)
    }
    
    # calculate the best stationary error var according to the data
    sigmaSqTest = seq(sigmaEpsSqNonMount, sigmaEpsSqMount, l=100)
    statDatScores = scoreFunDat(sigmaSqTest)
    bestI = which.min(statDatScores)
    sigmaSqBestDat = sigmaSqTest[bestI]
    scoreBestDat = statDatScores[bestI]
    
    out = optimize(scoreFunDat, lower=0, upper=1)
    sigmaSqBestDat = out$minimum
    scoreBestDat = out$objective
    
    # calculate the best stationary error var on the domain
    statDomScores = scoreFunDom(sigmaSqTest)
    bestI = which.min(statDomScores)
    sigmaSqBestDom = sigmaSqTest[bestI]
    scoreBestDom = statDomScores[bestI]
    
    out = optimize(scoreFunDom, lower=0, upper=1)
    sigmaSqBestDom = out$minimum
    scoreBestDom = out$objective
    
    # Calculate expected scores under the data
    sigmaSqsNM = seq(sigmaSqBestDom, sigmaEpsSqNonMount, l=100)
    sigmaSqsM = seq(sigmaSqBestDom, sigmaEpsSqMount, l=100)
    nonStatDatScores = scoreFunDat(sigmaSqsNM, sigmaSqsM)
    
    if(FALSE) {
      # t = seq(0, 1, l=100)
      # plot(t, nonStatDatScores, type="l")
      plot(sigmaSqsNM, nonStatDatScores, type="l", xlab="Plains error variance", ylab="Score")
      abline(h=scoreBestDat, lty=2)
    }
    
    # sigmaSqsNM = 10^(seq(log10(.1^2), log10(1^2), l=1000))
    # sigmaSqsM = 10^(seq(log10(.1^2), log10(1^2), l=1000))
    sigmaSqsNM = 10^(seq(log10(.01^2), log10(2.9^2), l=1000))
    sigmaSqsM = 10^(seq(log10(.01^2), log10(2.9^2), l=1000))
    
    if(is.null(saveResults)) {
      saveResults = length(sigmaSqsNM) > 200
    }
    
    sigmaSqsMat = expand.grid(sigmaSqsNM=sigmaSqsNM, sigmaSqsM=sigmaSqsM)
    if(regenResults) {
      nonStatDatScoresMat = scoreFunDat(sigmaSqsMat[,1], sigmaSqsMat[,2])
      nonStatDomScoresMat = scoreFunDom(sigmaSqsMat[,1], sigmaSqsMat[,2])
      if(saveResults) {
        save(sigmaSqsMat, nonStatDatScoresMat, nonStatDomScoresMat, 
             file=paste0("savedOutput/griddedCVtestNonstatErr/nonStatScoresMat_oversampM", 
                         oversampleMountRatio, ".RData"))
      }
    } else {
      out = load(paste0("savedOutput/griddedCVtestNonstatErr/nonStatScoresMat_oversampM", 
                        oversampleMountRatio, ".RData"))
      # sigmaSqsM = sort(unique(sigmaSqsMat[,2]))
      # sigmaSqsNM = sort(unique(sigmaSqsMat[,1]))
    }
    
    scoreDiffDat = nonStatDatScoresMat - scoreBestDat
    scoreDiffDom = nonStatDomScoresMat - scoreFunDom(sigmaSqBestDat)
    
    # calculate mean effective distance and mean distance from wrong model for LOO
    meanEffDist = sign(scoreDiffDat) * sign(scoreDiffDom) * apply(abs(cbind(scoreDiffDat, scoreDiffDom)), 1, min)
    fullMat = cbind(sigmaSqsMat, meanEffDist)
    cols = makeRedBlueDivergingColors(64, valRange=range(fullMat[,3]), center=0)
    
    pdf(paste0("figures/nonstatErrorIllustration/effDist_Score", score, "_oversampM", oversampleMountRatio, ".pdf"), width=5.1, height=5)
    par(mar=c(2.8, 1.3, 2, 1), oma=c(0, 0, 0, 1.5), mgp=c(1.9,.7,0))
    plot(log10(fullMat[,1]), log10(fullMat[,2]), type="n", axes=FALSE, xlab="", ylab="",
         main=TeX(paste0("LOO eff. pref to right model ($R_{oversamp}=", oversampleMountRatio, "$)")), asp=1)
    myQuiltPlot(log10(fullMat[,1]), log10(fullMat[,2]), fullMat[,3], col=cols,
                add=TRUE, nx=500, ny=500, legend.mar=1.8)
    correctModelIndexMat = matrix(meanEffDist, nrow=length(sigmaSqsNM), ncol=length(sigmaSqsM))
    contour(log10(sigmaSqsNM), log10(sigmaSqsM), correctModelIndexMat, col=rgb(.4,.4,.4), 
            add=TRUE)
    axis(1, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-.75)
    axis(2, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-1.5)
    mtext(TeX("$\\sigma_{epsilon,Mount}^2"), side=2, line=-.1)
    mtext(TeX("$\\sigma_{epsilon,Plains}^2"), side=1, line=1.3)
    points(log10(c(sigmaSqBestDat, .01)), log10(c(sigmaSqBestDat, 1)), pch=19, col=c("purple", "green"))
    # myQuiltPlot(fullMat[,1], fullMat[,2], fullMat[,3], col=cols, 
    #            nx=500, ny=500, log="xy", asp=1, addColorBar=TRUE, 
    #            xlab=TeX("$\\sigma_{epsilon,Plains}^2"), ylab=TeX("$\\sigma_{epsilon,Mount}^2"), 
    #            main="Correct model selected", legend.mar=0)
    dev.off()
    
    meanDist = abs(scoreDiffDat) * sign(scoreDiffDat) * sign(scoreDiffDom)
    fullMat = cbind(sigmaSqsMat, meanDist)
    
    cols = makeRedBlueDivergingColors(64, valRange=range(fullMat[,3]), center=0)
    pdf(paste0("figures/nonstatErrorIllustration/dist_Score", score, "_oversampM", oversampleMountRatio, ".pdf"), width=5.1, height=5)
    par(mar=c(2.8, 1.3, 2, 1), oma=c(0, 0, 0, 1.5), mgp=c(1.9,.7,0))
    plot(log10(fullMat[,1]), log10(fullMat[,2]), type="n", axes=FALSE, xlab="", ylab="",
         main=TeX(paste0("LOO pref to right model ($R_{oversamp}=", oversampleMountRatio, "$)")), asp=1)
    myQuiltPlot(log10(fullMat[,1]), log10(fullMat[,2]), fullMat[,3], col=cols,
                add=TRUE, nx=500, ny=500, legend.mar=1.8)
    correctBiasIndexMat = matrix(meanDist, nrow=length(sigmaSqsNM), ncol=length(sigmaSqsM))
    contour(log10(sigmaSqsNM), log10(sigmaSqsM), correctBiasIndexMat, col=rgb(.4,.4,.4), 
            add=TRUE)
    axis(1, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-.75)
    axis(2, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-1.5)
    mtext(TeX("$\\sigma_{epsilon,Mount}^2"), side=2, line=-.1)
    mtext(TeX("$\\sigma_{epsilon,Plains}^2"), side=1, line=1.3)
    points(log10(c(sigmaSqBestDat, .01)), log10(c(sigmaSqBestDat, 1)), pch=19, col=c("purple", "green"))
    # myQuiltPlot(fullMat[,1], fullMat[,2], fullMat[,3], col=cols, 
    #            nx=500, ny=500, log="xy", asp=1, addColorBar=TRUE, 
    #            xlab=TeX("$\\sigma_{epsilon,Plains}^2"), ylab=TeX("$\\sigma_{epsilon,Mount}^2"), 
    #            main="Correct model selected", legend.mar=0)
    dev.off()
    
    if(oversampleMountRatio == "0.2") {
      # select some comparison models. First get transect of preferences
      transectI = which.min((sigmaSqsM - .1)^2)
      print(sigmaSqsM[transectI]) # should be about .1
      colI = which.min(abs(correctBiasIndexMat[,transectI]) - as.numeric(sigmaSqsNM > .1^2))
      print(sigmaSqsNM[colI]) # about 0.035
      print(sqrt(sigmaSqsNM[colI])) # about 0.188
      sigmaSqMSims = rep(.1, 3)
      sigmaSqPSims = c(.1^2, 0.035, .1)
      
      
      cols = makeRedBlueDivergingColors(64, valRange=range(fullMat[,3]), center=0)
      pdf(paste0("figures/nonstatErrorIllustration/dist_simsScore", score, "_oversampM", oversampleMountRatio, ".pdf"), width=5.1, height=5)
      par(mar=c(2.8, 1.3, 2, 1), oma=c(0, 0, 0, 1.5), mgp=c(1.9,.7,0))
      plot(log10(fullMat[,1]), log10(fullMat[,2]), type="n", axes=FALSE, xlab="", ylab="",
           main=TeX(paste0("LOO pref to right model ($R_{oversamp}=", oversampleMountRatio, "$)")), asp=1)
      myQuiltPlot(log10(fullMat[,1]), log10(fullMat[,2]), fullMat[,3], col=cols,
                  add=TRUE, nx=500, ny=500, legend.mar=1.8)
      correctBiasIndexMat = matrix(meanDist, nrow=length(sigmaSqsNM), ncol=length(sigmaSqsM))
      contour(log10(sigmaSqsNM), log10(sigmaSqsM), correctBiasIndexMat, col=rgb(.4,.4,.4), 
              add=TRUE)
      axis(1, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-.75)
      axis(2, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-1.5)
      mtext(TeX("$\\sigma_{epsilon,Mount}^2"), side=2, line=-.1)
      mtext(TeX("$\\sigma_{epsilon,Plains}^2"), side=1, line=1.3)
      points(log10(c(sigmaSqBestDat, .01)), log10(c(sigmaSqBestDat, 1)), pch=19, col=c("purple", "green"))
      points(log10(sigmaSqPSims), log10(sigmaSqMSims), pch="x")
      # myQuiltPlot(fullMat[,1], fullMat[,2], fullMat[,3], col=cols, 
      #            nx=500, ny=500, log="xy", asp=1, addColorBar=TRUE, 
      #            xlab=TeX("$\\sigma_{epsilon,Plains}^2"), ylab=TeX("$\\sigma_{epsilon,Mount}^2"), 
      #            main="Correct model selected", legend.mar=0)
      dev.off()
    }
    
    # calculate mean effective distance and mean distance from wrong model for IW
    # noting that they are the same for IW
    meanEffDist = abs(scoreDiffDom)
    fullMat = cbind(sigmaSqsMat, meanEffDist)
    cols = makeRedBlueDivergingColors(64, valRange=c(0, max(fullMat[,3])), center=0)
    
    pdf(paste0("figures/nonstatErrorIllustration/effDistIW_Score", score, "_oversampM", oversampleMountRatio, ".pdf"), width=5.1, height=5)
    par(mar=c(2.8, 1.3, 2, 1), oma=c(0, 0, 0, 1.5), mgp=c(1.9,.7,0))
    plot(log10(fullMat[,1]), log10(fullMat[,2]), type="n", axes=FALSE, xlab="", ylab="",
         main=TeX(paste0("IW eff. pref to right model ($R_{oversamp}=", oversampleMountRatio, "$)")), asp=1)
    myQuiltPlot(log10(fullMat[,1]), log10(fullMat[,2]), fullMat[,3], col=cols,
                add=TRUE, nx=500, ny=500, legend.mar=1.8, zlim=c(0, max(fullMat[,3])))
    correctModelIndexMat = matrix(meanEffDist, nrow=length(sigmaSqsNM), ncol=length(sigmaSqsM))
    # contour(log10(sigmaSqsNM), log10(sigmaSqsM), correctModelIndexMat, col=rgb(.4,.4,.4), 
    #         add=TRUE, zlim=c(0, max(correctModelIndexMat)), 
    #         levels=c(5e-4, seq(.2, max(correctModelIndexMat), by=.2)), 
    #         labels=c("0", as.character(seq(.2, max(correctModelIndexMat), by=.2))))
    # correctModelIndexMat[correctModelIndexMat < 1e-4] = -1e6
    contour(log10(sigmaSqsNM), log10(sigmaSqsM), correctModelIndexMat, col=rgb(.4,.4,.4), 
            add=TRUE, zlim=c(0, max(correctModelIndexMat)), 
            levels=seq(0, max(correctModelIndexMat), by=.2))
    # # add 0 contour by hand since contour doesn't want to plot it...
    # minI = apply(correctModelIndexMat, 1, which.min)
    # toPlotI = c(which(minI < length(sigmaSqsM)), match(TRUE, minI == length(sigmaSqsM)))
    # minI = minI[toPlotI]
    # sigmaSqsMs0 = sigmaSqsM[minI]
    # lines(log10(sigmaSqsNM[toPlotI]), log10(sigmaSqsMs0), col=rgb(.4,.4,.4))
    axis(1, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-.75)
    axis(2, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-1.5)
    mtext(TeX("$\\sigma_{epsilon,Mount}^2"), side=2, line=-.1)
    mtext(TeX("$\\sigma_{epsilon,Plains}^2"), side=1, line=1.3)
    points(log10(c(sigmaSqBestDat, .01)), log10(c(sigmaSqBestDat, 1)), pch=19, col=c("purple", "green"))
    # myQuiltPlot(fullMat[,1], fullMat[,2], fullMat[,3], col=cols, 
    #            nx=500, ny=500, log="xy", asp=1, addColorBar=TRUE, 
    #            xlab=TeX("$\\sigma_{epsilon,Plains}^2"), ylab=TeX("$\\sigma_{epsilon,Mount}^2"), 
    #            main="Correct model selected", legend.mar=0)
    dev.off()
    
    thisStartTime = proc.time()[3]
    for(j in 1:length(ns)) {
      thisN = ns[j]
      
      # calculate asymptotic selection probabilities for this n
      if(regenResults) {
        selProbs = correctSelProb(thisN, sigmaEpsSqNonMountWrong1=fullMat[,1], sigmaEpsSqMountWrong1=fullMat[,2], 
                                  sigmaEpsSqNonMountWrong2=sigmaSqBestDat, sigmaEpsSqMountWrong2=sigmaSqBestDat)
        selProbsLOO = selProbs[,1]
        selProbsIW = selProbs[,2]
        
        if(saveResults) {
          save(selProbsLOO, selProbsIW, 
               file=paste0("savedOutput/griddedCVtestNonstatErr/nonStatSelProbs_oversampM", 
                           oversampleMountRatio, "_n", thisN, ".RData"))
        }
      } else {
        out = load(paste0("savedOutput/griddedCVtestNonstatErr/nonStatSelProbs_oversampM", 
                          oversampleMountRatio, "_n", thisN, ".RData"))
      }
      
      cols = makeRedBlueDivergingColors(64, valRange=c(0, 1), center=0.5)
      
      # LOO selection probabilities
      pdf(paste0("figures/nonstatErrorIllustration/correctModelProbLOO_Score", score, "_oversampM", oversampleMountRatio, "_n", thisN, ".pdf"), width=5.1, height=5)
      par(mar=c(2.8, 1.3, 2, 1), oma=c(0, 0, 0, 1.5), mgp=c(1.9,.7,0))
      plot(log10(fullMat[,1]), log10(fullMat[,2]), type="n", axes=FALSE, xlab="", ylab="",
           main=TeX(paste0("LOO correct probability ($R_{oversamp}=", oversampleMountRatio, "$, $n=$", thisN, ")")), asp=1)
      myQuiltPlot(log10(fullMat[,1]), log10(fullMat[,2]), selProbsLOO, col=cols,
                  add=TRUE, nx=500, ny=500, legend.mar=1.8, zlim=c(0, 1))
      selProbMat = matrix(selProbsLOO, nrow=length(sigmaSqsNM), ncol=length(sigmaSqsM))
      contour(log10(sigmaSqsNM), log10(sigmaSqsM), selProbMat, col=rgb(.4,.4,.4), 
              levels=seq(.2, .8, by=.2), add=TRUE)
      axis(1, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-.75)
      axis(2, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-1.5)
      mtext(TeX("$\\sigma_{epsilon,Mount}^2"), side=2, line=-.1)
      mtext(TeX("$\\sigma_{epsilon,Plains}^2"), side=1, line=1.3)
      points(log10(c(sigmaSqBestDat, .01)), log10(c(sigmaSqBestDat, 1)), pch=19, col=c("purple", "green"))
      # myQuiltPlot(fullMat[,1], fullMat[,2], fullMat[,3], col=cols, 
      #            nx=500, ny=500, log="xy", asp=1, addColorBar=TRUE, 
      #            xlab=TeX("$\\sigma_{epsilon,Plains}^2"), ylab=TeX("$\\sigma_{epsilon,Mount}^2"), 
      #            main="Correct model selected", legend.mar=0)
      dev.off()
      
      # IW selection probabilities
      pdf(paste0("figures/nonstatErrorIllustration/correctModelProbIW_Score", score, "_oversampM", oversampleMountRatio, "_n", thisN, ".pdf"), width=5.1, height=5)
      par(mar=c(2.8, 1.3, 2, 1), oma=c(0, 0, 0, 1.5), mgp=c(1.9,.7,0))
      plot(log10(fullMat[,1]), log10(fullMat[,2]), type="n", axes=FALSE, xlab="", ylab="",
           main=TeX(paste0("IW correct probability ($R_{oversamp}=", oversampleMountRatio, "$, $n=$", thisN, ")")), asp=1)
      myQuiltPlot(log10(fullMat[,1]), log10(fullMat[,2]), selProbsIW, col=cols,
                  add=TRUE, nx=500, ny=500, legend.mar=1.8, zlim=c(0, 1))
      selProbMat = matrix(selProbsIW, nrow=length(sigmaSqsNM), ncol=length(sigmaSqsM))
      contour(log10(sigmaSqsNM), log10(sigmaSqsM), selProbMat, col=rgb(.4,.4,.4), 
              levels=seq(.2, .8, by=.2), add=TRUE)
      axis(1, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-.75)
      axis(2, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-1.5)
      mtext(TeX("$\\sigma_{epsilon,Mount}^2"), side=2, line=-.1)
      mtext(TeX("$\\sigma_{epsilon,Plains}^2"), side=1, line=1.3)
      points(log10(c(sigmaSqBestDat, .01)), log10(c(sigmaSqBestDat, 1)), pch=19, col=c("purple", "green"))
      # myQuiltPlot(fullMat[,1], fullMat[,2], fullMat[,3], col=cols, 
      #            nx=500, ny=500, log="xy", asp=1, addColorBar=TRUE, 
      #            xlab=TeX("$\\sigma_{epsilon,Plains}^2"), ylab=TeX("$\\sigma_{epsilon,Mount}^2"), 
      #            main="Correct model selected", legend.mar=0)
      dev.off()
      
      # IW - LOO selection probabilities
      if(oversampleMountRatio != 1) {
        cols = makeRedBlueDivergingColors(64, range(selProbsIW-selProbsLOO), center=0)
        pdf(paste0("figures/nonstatErrorIllustration/correctModelProbIWvsLOO_Score", score, "_oversampM", oversampleMountRatio, "_n", thisN, ".pdf"), width=5.1, height=5)
        par(mar=c(2.8, 1.3, 2, 1), oma=c(0, 0, 0, 1.5), mgp=c(1.9,.7,0))
        plot(log10(fullMat[,1]), log10(fullMat[,2]), type="n", axes=FALSE, xlab="", ylab="",
             main=TeX(paste0("IW-LOO correct probability ($R_{oversamp}=", oversampleMountRatio, "$, $n=$", thisN, ")")), asp=1)
        myQuiltPlot(log10(fullMat[,1]), log10(fullMat[,2]), selProbsIW-selProbsLOO, col=cols,
                    add=TRUE, nx=500, ny=500, legend.mar=1.8)
        selProbMat = matrix(selProbsIW-selProbsLOO, nrow=length(sigmaSqsNM), ncol=length(sigmaSqsM))
        contour(log10(sigmaSqsNM), log10(sigmaSqsM), selProbMat, col=rgb(.4,.4,.4), 
                levels=seq(-.8, .8, by=.2), add=TRUE)
        axis(1, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-.75)
        axis(2, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-1.5)
        mtext(TeX("$\\sigma_{epsilon,Mount}^2"), side=2, line=-.1)
        mtext(TeX("$\\sigma_{epsilon,Plains}^2"), side=1, line=1.3)
        points(log10(c(sigmaSqBestDat, .01)), log10(c(sigmaSqBestDat, 1)), pch=19, col=c("purple", "green"))
        # myQuiltPlot(fullMat[,1], fullMat[,2], fullMat[,3], col=cols, 
        #            nx=500, ny=500, log="xy", asp=1, addColorBar=TRUE, 
        #            xlab=TeX("$\\sigma_{epsilon,Plains}^2"), ylab=TeX("$\\sigma_{epsilon,Mount}^2"), 
        #            main="Correct model selected", legend.mar=0)
        dev.off()
      }
      
      # model selection score (positive oriented expectation of true score of selected model - unselected model)
      pdf(paste0("figures/nonstatErrorIllustration/modelSelectScoreLOO_Score", score, "_oversampM", oversampleMountRatio, "_n", thisN, ".pdf"), width=5.1, height=5)
      par(mar=c(2.8, 1.3, 2, 1), oma=c(0, 0, 0, 1.5), mgp=c(1.9,.7,0))
      plot(log10(fullMat[,1]), log10(fullMat[,2]), type="n", axes=FALSE, xlab="", ylab="",
           main=TeX(paste0("LOO selection score ($R_{oversamp}=", oversampleMountRatio, "$, $n=$", thisN, ")")), asp=1)
      selProbMat = matrix(selProbsLOO, nrow=length(sigmaSqsNM), ncol=length(sigmaSqsM))
      absDiffMat = matrix(abs(scoreDiffDom), nrow=length(sigmaSqsNM), ncol=length(sigmaSqsM))
      selectScoreLOO = absDiffMat * (2*selProbMat - 1)
      cols = makeRedBlueDivergingColors(64, valRange=range(selectScoreLOO), center=0)
      myQuiltPlot(log10(fullMat[,1]), log10(fullMat[,2]), selectScoreLOO, col=cols,
                  add=TRUE, nx=500, ny=500, legend.mar=1.8)
      contour(log10(sigmaSqsNM), log10(sigmaSqsM), selectScoreLOO, col=rgb(.4,.4,.4), add=TRUE)
      axis(1, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-.75)
      axis(2, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-1.5)
      mtext(TeX("$\\sigma_{epsilon,Mount}^2"), side=2, line=-.1)
      mtext(TeX("$\\sigma_{epsilon,Plains}^2"), side=1, line=1.3)
      points(log10(c(sigmaSqBestDat, .01)), log10(c(sigmaSqBestDat, 1)), pch=19, col=c("purple", "green"))
      # myQuiltPlot(fullMat[,1], fullMat[,2], fullMat[,3], col=cols, 
      #            nx=500, ny=500, log="xy", asp=1, addColorBar=TRUE, 
      #            xlab=TeX("$\\sigma_{epsilon,Plains}^2"), ylab=TeX("$\\sigma_{epsilon,Mount}^2"), 
      #            main="Correct model selected", legend.mar=0)
      dev.off()
      
      pdf(paste0("figures/nonstatErrorIllustration/modelSelectScoreIW_Score", score, "_oversampM", oversampleMountRatio, "_n", thisN, ".pdf"), width=5.1, height=5)
      par(mar=c(2.8, 1.3, 2, 1), oma=c(0, 0, 0, 1.5), mgp=c(1.9,.7,0))
      plot(log10(fullMat[,1]), log10(fullMat[,2]), type="n", axes=FALSE, xlab="", ylab="",
           main=TeX(paste0("IW selection score ($R_{oversamp}=", oversampleMountRatio, "$, $n=$", thisN, ")")), asp=1)
      selProbMat = matrix(selProbsIW, nrow=length(sigmaSqsNM), ncol=length(sigmaSqsM))
      absDiffMat = matrix(abs(scoreDiffDom), nrow=length(sigmaSqsNM), ncol=length(sigmaSqsM))
      selectScoreIW = absDiffMat * (2*selProbMat - 1)
      cols = makeRedBlueDivergingColors(64, valRange=range(selectScoreIW), center=0)
      myQuiltPlot(log10(fullMat[,1]), log10(fullMat[,2]), selectScoreIW, col=cols,
                  add=TRUE, nx=500, ny=500, legend.mar=1.8)
      contour(log10(sigmaSqsNM), log10(sigmaSqsM), selectScoreIW, col=rgb(.4,.4,.4), add=TRUE)
      axis(1, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-.75)
      axis(2, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-1.5)
      mtext(TeX("$\\sigma_{epsilon,Mount}^2"), side=2, line=-.1)
      mtext(TeX("$\\sigma_{epsilon,Plains}^2"), side=1, line=1.3)
      points(log10(c(sigmaSqBestDat, .01)), log10(c(sigmaSqBestDat, 1)), pch=19, col=c("purple", "green"))
      # myQuiltPlot(fullMat[,1], fullMat[,2], fullMat[,3], col=cols, 
      #            nx=500, ny=500, log="xy", asp=1, addColorBar=TRUE, 
      #            xlab=TeX("$\\sigma_{epsilon,Plains}^2"), ylab=TeX("$\\sigma_{epsilon,Mount}^2"), 
      #            main="Correct model selected", legend.mar=0)
      dev.off()
      
      if(oversampleMountRatio != 1) {
        pdf(paste0("figures/nonstatErrorIllustration/modelSelectScoreIWvsLOO_Score", score, "_oversampM", oversampleMountRatio, "_n", thisN, ".pdf"), width=5.1, height=5)
        par(mar=c(2.8, 1.3, 2, 1), oma=c(0, 0, 0, 1.5), mgp=c(1.9,.7,0))
        plot(log10(fullMat[,1]), log10(fullMat[,2]), type="n", axes=FALSE, xlab="", ylab="",
             main=TeX(paste0("IW-LOO selection score ($R_{oversamp}=", oversampleMountRatio, "$, $n=$", thisN, ")")), asp=1)
        selectScoreDiff = selectScoreIW - selectScoreLOO
        cols = makeRedBlueDivergingColors(64, valRange=range(selectScoreDiff), center=0)
        myQuiltPlot(log10(fullMat[,1]), log10(fullMat[,2]), selectScoreDiff, col=cols,
                    add=TRUE, nx=500, ny=500, legend.mar=1.8)
        contour(log10(sigmaSqsNM), log10(sigmaSqsM), selectScoreDiff, col=rgb(.4,.4,.4), add=TRUE)
        axis(1, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-.75)
        axis(2, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-1.5)
        mtext(TeX("$\\sigma_{epsilon,Mount}^2"), side=2, line=-.1)
        mtext(TeX("$\\sigma_{epsilon,Plains}^2"), side=1, line=1.3)
        points(log10(c(sigmaSqBestDat, .01)), log10(c(sigmaSqBestDat, 1)), pch=19, col=c("purple", "green"))
        # myQuiltPlot(fullMat[,1], fullMat[,2], fullMat[,3], col=cols, 
        #            nx=500, ny=500, log="xy", asp=1, addColorBar=TRUE, 
        #            xlab=TeX("$\\sigma_{epsilon,Plains}^2"), ylab=TeX("$\\sigma_{epsilon,Mount}^2"), 
        #            main="Correct model selected", legend.mar=0)
        dev.off()
      }
      
      if(oversampleMountRatio == "0.2") {
        
        pdf(paste0("figures/nonstatErrorIllustration/modelSelectScoreIWvsLOOsims_Score", score, "_oversampM", oversampleMountRatio, "_n", thisN, ".pdf"), width=5.1, height=5)
        par(mar=c(2.8, 1.3, 2, 1), oma=c(0, 0, 0, 1.5), mgp=c(1.9,.7,0))
        plot(log10(fullMat[,1]), log10(fullMat[,2]), type="n", axes=FALSE, xlab="", ylab="",
             main=TeX(paste0("IW-LOO selection score ($R_{oversamp}=", oversampleMountRatio, "$, $n=$", thisN, ")")), asp=1)
        selectScoreDiff = selectScoreIW - selectScoreLOO
        cols = makeRedBlueDivergingColors(64, valRange=range(selectScoreDiff), center=0)
        myQuiltPlot(log10(fullMat[,1]), log10(fullMat[,2]), selectScoreDiff, col=cols,
                    add=TRUE, nx=500, ny=500, legend.mar=1.8)
        points(log10(sigmaSqPSims), log10(sigmaSqMSims), pch="x")
        contour(log10(sigmaSqsNM), log10(sigmaSqsM), selectScoreDiff, col=rgb(.4,.4,.4), add=TRUE)
        axis(1, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-.75)
        axis(2, at=seq(-4, 0, by=2), labels=c(".01^2", ".1^2", "1"), line=-1.5)
        mtext(TeX("$\\sigma_{epsilon,Mount}^2"), side=2, line=-.1)
        mtext(TeX("$\\sigma_{epsilon,Plains}^2"), side=1, line=1.3)
        points(log10(c(sigmaSqBestDat, .01)), log10(c(sigmaSqBestDat, 1)), pch=19, col=c("purple", "green"))
        # myQuiltPlot(fullMat[,1], fullMat[,2], fullMat[,3], col=cols, 
        #            nx=500, ny=500, log="xy", asp=1, addColorBar=TRUE, 
        #            xlab=TeX("$\\sigma_{epsilon,Plains}^2"), ylab=TeX("$\\sigma_{epsilon,Mount}^2"), 
        #            main="Correct model selected", legend.mar=0)
        dev.off()
        
        # calculate select scores and difference in scores
        transectI = which.min((sigmaSqsM - .1)^2)
        sigmaSqMSims = rep(.1, 3)
        sigmaSqPSims = c(.1^2, 0.035, .1)
        colI1 = which.min(abs(sigmaSqsNM - sigmaSqPSims[1]))
        colI2 = which.min(abs(sigmaSqsNM - sigmaSqPSims[2]))
        colI3 = which.min(abs(sigmaSqsNM - sigmaSqPSims[3]))
        colIs = c(colI1, colI2, colI3)
        print(sigmaSqsNM[colIs]) # about 0.01, .035, .1
        print(sqrt(sigmaSqsNM[colIs])) # about .100, 0.188, .317
        
        
        print(paste("", round(selectScoreLOO[colIs, transectI], digits=3)))
        # n=50:
        # [1] " 0.457"  " 0.008"  " -0.134"
        # n=500:
        # [1] " 0.457"  " 0.026"  " -0.134"
        # n=10000:
        # " 0.457"  " 0.114"  " -0.134"
        print(paste("", round(selectScoreIW[colIs, transectI], digits=3)))
        # n=50:
        # [1] " 0.457" " 0.359" " 0.086"
        # n=500:
        # " 0.457" " 0.364" " 0.134"
        # n=10000:
        # " 0.457" " 0.364" " 0.134"
        print(paste("", round(selectScoreDiff[colIs, transectI], digits=3)))
        # n=50:
        # [1] " 0"     " 0.351" " 0.22" 
        # n=500:
        # [1] " 0"     " 0.337" " 0.268"
        # n=10000:
        # [1] " 0"     " 0.249" " 0.268"
        browser()
      }
      
      thisCurrTime = proc.time()[3]
      estT = estTimeLeft(thisStartTime, thisCurrTime, j, length(ns))
      print(paste0("iter j=", j, "/", length(ns), ", est minutes left: ", round(estT/60, 1)))
    }
    
    currTime = proc.time()[3]
    estT = estTimeLeft(startTime, currTime, i, length(oversampleMountRatioSeq))
    print(paste0("est minutes left: ", round(estT/60, 1)))
  }
  
  invisible(NULL)
}


