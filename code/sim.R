
# simulate GRF on unit square. See ?circulantEmbedding from fields
rGRF = function(nx=100, ny=100, mu=0, sigma=1, 
                cov.args=list(Covariance="Matern", range=0.3, smoothness=1.0), 
                delta=5, sigmaEpsSq=0) {
  require(fields)
  
  xres = 1/nx
  yres = 1/ny
  xs = seq(xres/2, 1-xres/2, l=nx)
  ys = seq(yres/2, 1-yres/2, l=ny)
  locs = make.surface.grid(list(x = xs, y=ys))
  
  obj<- circulantEmbeddingSetup(list(x = xs, y=ys), cov.args=cov.args, delta=delta)
  out = circulantEmbedding(obj)
  truth = mu + out*sigma
  obs = truth + rnorm(length(out), sd=sqrt(sigmaEpsSq))
  
  list(locs=locs, truth=truth, obs=obs, xs=xs, ys=ys, nx=nx, ny=ny, cov.args=cov.args, mu=mu, sigma=sigma)
}

# NOTE: SigmaAA can be a matrix with 1 column of the marginal variances if 
#       getFullCov == FALSE, and is not necessary to include if 
#       getCondVar == FALSE.
condMeanMVN = function(SigmaAA=NULL, SigmaAB, SigmaBB, ysB, getFullCov=TRUE, 
                       muA=rep(0, nrow(SigmaAB)), muB=rep(0, length(ysB)), 
                       U=NULL, returnU=FALSE, getCondVar=TRUE) {
  
  # first get predictive/conditional mean:
  # muA|B = muA + SigmaAB SigmaBB^-1 (yB - muB)
  # muAcondBtrue = muTrue + SigmaABtrue %*% solve(SigmaBBtrue, residsTrue)
  # muAcondBfalse = muFalse + SigmaABfalse %*% solve(SigmaBBfalse, residsFalse)
  if(is.null(U)) {
    U = chol(SigmaBB)
  }
  
  muAcondB = muA + SigmaAB %*% t(U) %*% backsolve(U, ysB - muB)
  
  # now get predictive/conditional variance if requested by user:
  if(getCondVar) {
    # SigmaA|B = SigmaAA - SigmaAB SigmaBB^-1 SigmaBA
    #          = SigmaAA - SigmaAB (U' U)^-1 SigmaBA
    #          = SigmaAA - SigmaAB U^-1 (U')^-1 SigmaBA
    #          = SigmaAA - R' R
    RA = forwardsolve(t(U), t(SigmaAB))
    
    if(getFullCov) {
      varAcondB = SigmaAA - t(RA) %*% RA
    } else {
      # (R' R)_ii = (R')_i: R_:i
      varAcondB = myDiag(SigmaAA) - myDiag(apply(RA, 2, function(x) {sum(x^2)}))
    }
  }
  
  out = list(muAcondB=muAcondB)
  if(getCondVar) {
    out = c(out, list(varAcondB=varAcondB))
  }
  if(returnU) {
    out = c(out, list(U=U))
  }
  out
}

# Show why gridded cross validation fails as function of grid resolution. Fix a 
# sampling distribution that is variable on unit square based on GRF, and for a 
# number of different grid resolutions. Error will vary between LOO-CV error 
# and asymptotic extrapolation error, the true predictive error being somewhere 
# in between
#   1.  for iter in 1:number of simulations:
#   2.    Simulate 1 GRF on unit square for responses
#   3.    Simulate 1 GRF on unit square for sample distribution
#   4.    Simulate observations in domain based on the sampling distn and the GRF 
#         for responses (don't include nugget)
#   5.    Calculate full covariate matrix for sample and crossCov to all locs
#   6.    for i in 1:number of block resolutions:
#   7.      Group observations with blocks
#   8.      Leave out one cell at a time, calculate gridded CV MSE
#   9.    Leave out one observation at a time, calculate LOO-CV MSE
#   10.   Calculate actual predictive MSE
#   11.   Save result
#   12. Combine and save all results
# Note that the below function is a single iteration of the for loop 
# for easy parallelization.
griddedResTestIter = function(rGRFargsTruth=NULL, rGRFargsSample=NULL, 
                              n=50, grindNs=2^(1:6), iter=1, seed=123, 
                              nx=100, ny=100, sigmaEpsSq=0) {
  
  set.seed(seed)
  
  # set default GRF parameters
  if(is.null(rGRFargsTruth)) {
    rGRFargsTruth = list(mu=0, sigma=1, 
                         cov.args=list(Covariance="Matern", range=0.2, smoothness=1.0), 
                         delta=5, sigmaEpsSq=0, nx=nx, ny=ny)
  }
  if(is.null(rGRFargsSample)) {
    # smaller range than the truth for more clumpy observation locations (same as rGRFargsPop)
    rGRFargsSample = list(mu=0, sigma=3, 
                          cov.args=list(Covariance="Matern", range=0.1, smoothness=1.0), 
                          delta=3, sigmaEpsSq=0, nx=nx, ny=ny)
  }
  
  # 2. Simulate 1 GRF for responses ----
  #    on unit square
  truthGRF = do.call("rGRF", rGRFargsTruth)
  truth = truthGRF$truth
  locs = sampleGRF$locs
  
  # 3. Simulate sample distribution ----
  #    as 1 GRF on unit square
  sampleGRF = do.call("rGRF", rGRFargsTruth)
  sampleRates = exp(sampleGRF$truth)
  sampleProbs = sampleRates * (1/sum(sampleRates))
  
  # 4. Simulate observations ----
  #    in domain based on the sampling distn and the GRF 
  #    for responses (don't include nugget?)
  gridResX = 1/nx
  gridResy = 1/ny
  sampleI = sample(1:nrow(locs), n, probs=sampleProbs, replace=TRUE)
  xs = locs[sampleI,] + cbind(runif(n, max=gridResX)-gridResX/2, 
                                      runif(n, max=gridResY)-gridResY/2)
  ys = truth[sampleI] + rnorm(n, sd=sqrt(sigmaEpsSq))
  
  # 5. Get covariance matrices ----
  #    for sample and crossCov to all locs
  
  SigmaSample = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsTruth$range, 
                               smoothness=rGRFargsTruth$smoothness)
  SigmaGridToSample = stationary.cov(locs, xs, Covariance="Matern", aRange=rGRFargsTruth$range, 
                               smoothness=rGRFargsTruth$smoothness)
  
  # 6. for i in 1:number block resolutions: ----
  
  # determine grid cells associated with observations
  getCellI = function(thisLoc) {
    xI = which(thisLoc[1] > highs) - 1
    yI = which(thisLoc[2] > highs) - 1
    
    thisInd = indMat[yI, xI]
    thisInd %in% testCellInds
  }
  
  griddedCVs = numeric(length(gridNs))
  for(i in 1:length(gridNs)) {
    gridN = gridNs[i]
    print(paste0("grid resolution ", gridN, " (", i, "/", length(gridNs), ")"))
    
    # 7. Group data by block ----
    highs = seq(0, 1, l=gridN+1)[-1]
    cellIsSample = apply(sampleLocs, 1, getCellI)
    uniqueCellIs = sort(unique(cellIsSample))
    blockCVs = numeric(length(uniqueCellIs))
    
    # 8. Get gridded CV MSE ----
    for(j in 1:length(uniqueCellIs)) {
      print(paste0("Leaving out cell ", j, "/", length(uniqueCellIs)))
      
      cellI = uniqueCellIs[j]
      
      # separate data in test/train
      isTest = cellIsSample == cellI
      testYs = ys[isTest]
      trainYs = ys[!isTest]
      
      # predict test data
      SigmaAB = SigmaSample[isTest, !isTest]
      SigmaBB = SigmaSample[!isTest, !isTest]
      SigmaAA = SigmaSample[isTest, isTest]
      
      condDistn = condMeanMVN(SigmaAA=SigmaAA, SigmaAB=SigmaAB, SigmaBB=SigmaBB, 
                              ysB=trainYs, getFullCov=FALSE, getAnyCov=FALSE, 
                              getCondVar=FALSE)
      muAcondB = condDistn$muAcondB
      
      # calculate MSE for the grid cell
      blockCVs[j] = mean((testYs - muAcondB)^2)
    }
    
    # average over MSEs of each grid cell
    griddedCVs[i] = mean(blockCVs)
  }
  
  # 9. get LOO-CV MSE ----
  LOOCVs = numeric(n)
  for(i in 1:n) {
    print(paste0("Leaving out obs ", i, "/", n))
    
    # separate data in test/train
    isTest = (1:n) == i
    testYs = ys[isTest]
    trainYs = ys[!isTest]
    
    # predict test data
    SigmaAB = matrix(SigmaSample[isTest, !isTest], nrow=1)
    SigmaBB = SigmaSample[!isTest, !isTest]
    SigmaAA = SigmaSample[isTest, isTest]
    
    condDistn = condMeanMVN(SigmaAA=SigmaAA, SigmaAB=SigmaAB, SigmaBB=SigmaBB, 
                            ysB=trainYs, getFullCov=FALSE, getAnyCov=FALSE, 
                            getCondVar=FALSE)
    muAcondB = condDistn$muAcondB
    
    # calculate MSE for the grid cell
    LOOCVs[i] = mean((testYs - muAcondB)^2)
  }
  
  # 10. Calculate true MSE ----
  condDistn = condMeanMVN(SigmaAA=rep(1, nrow(locs)), SigmaAB=SigmaGridToSample, SigmaBB=SigmaSample, 
                          ysB=ys, getFullCov=FALSE, getAnyCov=FALSE, 
                          getCondVar=FALSE)
  muAcondB = condDistn$muAcondB
  
  # calculate MSE for the grid cell
  trueMSE = mean((truth - muAcondB)^2)
  
  # 11. Save result ----
  # rGRFargsTruth=NULL, rGRFargsSample=NULL, 
  # n=50, grindNs=2^(1:6), iter=1, seed=123, 
  # nx=100, ny=100, sigmaEpsSq=0
  save(trueMSE, LOOCVs, griddedCVs, grindNs, 
       iter, seed, rGRFargsTruth, rGRFargsSample, n, nx, ny, sigmaEpsSq, seed)
}

# fit 2 GRFs on field unit square, one with true parameters, other with false parameters. 
# Goal is twofold:
#  1) to be able to tell the true model from the false model
#  2) to be able to estimate the expected value of the scoring rule on the unit square
# Test a variety of:
#  1) Sampling schemes: uniform and non-uniform ("unif" and "lgcp" are 
#     currently supported)
#  2) Estimators of the sampling distribution: currently only the true 
#     distribution is supported ("truth"). These are not always used, however, 
#     depending on the validation type
#  3) Population distributions: What we wish to calculate predictive error with 
#     respect to. Can be either "unif" (uniform), "sampleDistn" (the sampling 
#     distribution, i.e. we assume sampling is proportional to the population), 
#     or a log Gaussian process that is not identical to the sampling 
#     distribution, but with identical distribution parameters.
#  4) Validation types: a method for partitioning sample into train and test 
#     data, and for weighting and combining prediction scores for the test data. 
#     ("basic", "thinning", "gridded", "wtGrid", "IS", and "PSIS" are currently 
#     supported)
# Given a sampling scheme and a true model, simulates a GRF and samples from it. 
# For that sample, tests all applicable combinations of sampling distn 
# estimators and validation types. Note that not all validation types can be 
# combined with all sampling distribution estimators. Combinations not supported 
# are skipped with scores replaced by NAs or by those of an identical method 
# when applicable:
#   - (sample scheme / est. of sample scheme / pop distn / valid. type)
#   - "unif" + "sampleDistn" replaced by "unif" + "unif"
#   - "vcell" + "basic" and "kernel" + "basic" replaced by "truth" + "basic"
#   - "vcell" + "gridded" and "kernel" + "gridded" replaced by "truth" + "gridded"
#   - "unif" + "unif" + "thinning" replaced by "unif" + "unif" + "basic"
#   - "sampleDistn" + "thinning" replaced by "sampleDistn" + "basic"
# Note that in some cases "sampleDistn" + "thinning" may result in NA scores 
# due to large differences between sample and pop distributions
unitSquareTestIter = function(rGRFargsTruth=NULL, rGRFargsFalse=NULL, 
                              rGRFargsSample=NULL, ntot=500, ntest=50, 
                              rGRFargsPop=NULL, 
                              # sampleRateEstMethods=c("truth", "vcell", "kernel"), 
                              # validationTypes=c("basic", "thinning", "gridded", "IS", "PSIS"), 
                              sampleScheme=c("lgcp", "unif", "horizStrip"), 
                              sampleRateEstMethods=c("truth"), 
                              popDistns=c("unif", "sampleDistn", "lgp"), 
                              validationTypes=c("basic", "thinning", "gridded", "wtGrid", "IS", "PSIS"), 
                              gridRes=.05, nx=100, ny=100, seed=NULL, fileNameRoot="uSqTest", 
                              iter=1) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # make sure input args make sense
  if(nx != ny) {
    stop("nx must be the same as ny")
  }
  if((nx %% (1/gridRes)) != 0) {
    stop("nx must be divisible by 1/gridRes so each grid cell has same number of points")
  }
  
  # Set input args ----
  sampleScheme = match.arg(sampleScheme)
  validationType = match.arg(validationType)
  ntrain = ntot - ntest
  
  fracTest = ntest/ntot # only used for gridded validation, since ntest is random for that method
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  if(is.null(rGRFargsFalse)) {
    rGRFargsFalse = list(mu=.1, sigma=1, 
                         cov.args=list(Covariance="Matern", range=0.2, smoothness=1.0), 
                         delta=5, sigmaEpsSq=0, nx=nx, ny=ny)
  }
  if(is.null(rGRFargsTruth)) {
    rGRFargsTruth = list(mu=0, sigma=1, 
                         cov.args=list(Covariance="Matern", range=0.2, smoothness=1.0), 
                         delta=5, sigmaEpsSq=0, nx=nx, ny=ny)
  }
  if(is.null(rGRFargsPop)) {
    # smaller range than the truth for more clumpy observation locations (same as rGRFargsSample)
    rGRFargsPop = list(mu=0, sigma=3, 
                       cov.args=list(Covariance="Matern", range=0.1, smoothness=1.0), 
                       delta=3, sigmaEpsSq=0, nx=nx, ny=ny)
  }
  
  # Generate true field ----
  truth = do.call("rGRF", rGRFargsTruth)
  
  # Do spatial sampling ----
  if(sampleScheme == "lgcp") {
    if(is.null(rGRFargsSample)) {
      # smaller range than the truth for more clumpy observation locations (same as rGRFargsPop)
      rGRFargsSample = list(mu=0, sigma=3, 
                            cov.args=list(Covariance="Matern", range=0.1, smoothness=1.0), 
                            delta=3, sigmaEpsSq=0, nx=nx, ny=ny)
    }
    
    # generate log Gaussian Cox process
    sampleGRF = do.call("rGRF", rGRFargsSample)
    lambdas = exp(sampleGRF$truth)
    lambdas = lambdas * (ntot/sum(lambdas))
    # sampleI = sample(1:length(lambdas), ntot, prob=lambdas/sum(lambdas))
    sampleI = which(UPmidzuno(lambdas))
    sampleLambdas = lambdas[sampleI]
  } else if(sampleScheme == "unif") {
    lambdas = rep(ntot/(truth$nx*truth$ny), truth$nx*truth$ny)
    sampleLambdas = rep(1, ntot)
    sampleI = sample(1:length(lambdas), ntot)
  } else if(sampleScheme == "horizStrip") {
    # a uniform sampled horizontal strip across the middle third of the unit square
    stop("sampleScheme 'horizStrip' not currently supported")
    # sampleLambdas = rep(1, ntot)
    # sampleI = which((truth$locs[,2] > 1/3) & (truth$locs[,2] < 1/3))
    # sampleI = sample(sampleI, ntot)
  }
  
  # get samples of the true field
  sampleLocs = truth$locs[sampleI,]
  sampleYs = truth$truth[sampleI]
  
  # Loop: pop distns ----
  for(k in 1:length(popDistns)) {
    popDistn = popDistns[k]
    print(paste0("pop distn: ", popDistn, " (", k, "/", length(popDistns), ")"))
    
    if((sampleScheme == "unif") && (popDistn == "sampleDistn")) {
      # "unif" + "sampleDistn" replaced by "unif" + "unif"
      next
    }
    
    # Get population distn ----
    # popDistns=c("unif", "sampleDistn", "lgp")
    if(popDistn == "unif") {
      lambdasPop = rep(1, truth$nx*truth$ny) / (truth$nx*truth$ny)
      sampleLambdasPop = rep(1, ntot)
    } else if(popDistn == "sampleDistn") {
      lambdasPop = lambdas
      sampleLambdasPop = sampleLambdas
    } else if(popDistn == "lgp") {
      # simulate the population (log Gaussian process)
      popGRF = do.call("rGRF", rGRFargsPop)
      lambdasPop = exp(popGRF$truth)
      sampleLambdasPop = lambdasPop[sampleI]
    }
    
    # Loop: sample distn ests ----
    for(i in 1:length(sampleRateEstMethods)) {
      sampleRateEstMethod = sampleRateEstMethods[i]
      print(paste0("sample distn est: ", sampleRateEstMethod, " (", i, "/", length(sampleRateEstMethods), ")"))
      
      # Estimate sample probability ----
      # (or at least something roughly proportional to it)
      if(sampleRateEstMethod == "truth") {
        ps = sampleLambdas
      } else if(sampleRateEstMethod == "vcell") {
        stop("Voronoi cell sample rate estimates not yet supported")
      } else if(sampleRateEstMethod == "kernel") {
        stop("kernel sample rate estimates not yet supported")
      }
      # ps = ps/sum(ps)
      
      # Loop: validation types ----
      for(j in 1:length(validationTypes)) {
        validationType = validationTypes[j]
        print(paste0("validation type: ", validationType, " (", j, "/", length(validationTypes), ")"))
        
        if((sampleRateEstMethod == "kernel") && (validationType == "basic")) {
          # "vcell" + "basic" and "kernel" + "basic" replaced by "truth" + "basic"
          next
        } else if((sampleRateEstMethod == "vcell") && (validationType == "basic")) {
          # "vcell" + "basic" and "kernel" + "basic" replaced by "truth" + "basic"
          next
        } else if((sampleRateEstMethod == "kernel") && (validationType == "gridded")) {
          # "vcell" + "gridded" and "kernel" + "gridded" replaced by "truth" + "gridded"
          next
        } else if((sampleRateEstMethod == "vcell") && (validationType == "gridded")) {
          # "vcell" + "gridded" and "kernel" + "gridded" replaced by "truth" + "gridded"
          next
        } else if((sampleScheme == "unif") && (popDistn == "unif") && (validationType == "thinning")) {
          # "unif" + "unif" + "thinning" replaced by "unif" + "unif" + "basic"
          next
        } else if((popDistn == "sampleDistn") && (validationType == "thinning")) {
          # "sampleDistn" + "thinning" replaced by "sampleDistn" + "basic"
          next
        }
        
        # Partition sample into train/test ----
        # also calculate validation weights
        if(validationType == "basic") {
          # SRS sampling of test indices from sample. Don't account for sample rates or pop distn
          testI = sample(1:ntot, ntest, replace=FALSE)
          testPs = rep(ntest/ntot, ntest)
          ws = rep(1/ntest, ntest)
        } else if(validationType == "thinning") {
          # sample test indices with probability proportional to inverse of sample 
          # probability, ps. It should sum to ntest
          testPs = 1/ps
          testPs = testPs/sum(testPs) * ntest
          testI = UPmidzuno(testPs)
          ws = rep(1/ntest, ntest)
        } else if((validationType == "gridded") || validationType == "wtGrid") {
          # divide domain into grid of resolution gridRes. Use fracTest of grid 
          # cells for test data
          
          gridN = ceiling(1/gridRes)
          pointRes = 1/nx
          highs = seq(0, 1, l=gridN+1)[-1]
          nCellsTot = length(highs)^2
          nCellsTest = fracTest * nCellsTot
          
          nCellsTestLow = floor(nCellsTest)
          if(nCellsTestLow != nCellsTest) {
            # set random number of cells if fracTest corresponds to non-integer number 
            # of cells to ensure E[ntest] is correct (so this is comparable to other 
            # methods)
            prob = nCellsTest - nCellsTestLow
            nCellsTest = nCellsTestLow + as.numeric(runif(1) < prob)
          }
          
          # index cells and select test/train cells:
          # cell i is in:
          # row floor(i) + 1 (row counts upward along unit square)
          # col floor(i/row) + i
          indMat = matrix(1:nCellsTot, nrow=gridN)
          testCellInds = sample(1:nCellsTot, nCellsTest)
          
          # determine grid cells associated with observations
          getCellI = function(thisLoc) {
            xI = which(thisLoc[1] > highs) - 1
            yI = which(thisLoc[2] > highs) - 1
            
            thisInd = indMat[yI, xI]
            thisInd %in% testCellInds
          }
          cellIsSample = apply(sampleLocs, 1, getCellI)
          testI = which(cellIs %in% testCellInds)
          cellIsTest = cellIsSample[testI]
          
          # calculate sample prob of each test cell
          if(sampleScheme == "unif") {
            # uniform sample case is simple
            testCellProb = rep(nCellsTest / nCellsTot, nCellsTot)
          } else if(sampleScheme == "lgcp") {
            # numerically integrate sampling rates within each grid cell
            
            # first get which rates are associated with which cell
            allCellInds = apply(truth$locs, 1, getCellI)
            
            # now integrate/sum sampling rates in each cell
            out = aggregate(lambdas, by=list(cellInd=allCellInds), FUN=sum)
            testCellProb = out$x / sum(out$x)
          }
          
          # get the sampling probably of the cell associated with each test observation
          testCellProbExtended = testCellProb[cellIsTest]
          
          # get inclusion probabilities and weights
          nPerTestCell = sapply(testCellInds, function(x) {sum(cellIsTest == x)})
          testPs = testCellProbExtended
          nPerTestCellExtended = nPerTestCell[match(cellIsTest, testCellInds)]
          testInclusionProb = nCellsTest / nCellsTot
          ws = testInclusionProb * (1/nPerTestCellExtended)
        } else if(validationType == "IS") {
          # SRS sampling of test indices from sample
          testI = sample(1:ntot, ntest, replace=FALSE)
          testPs = rep(ntest/ntot, ntest)
          # test inclusion prob = (sample inclusion prob) * (test inclusion prob)
          testInclusionProb = ps[testI] * testPs
          ws = 1/testInclusionProb
        } else if(validationType == "PSIS") {
          # SRS sampling of test indices from sample
          testI = sample(1:ntot, ntest, replace=FALSE)
          testPs = rep(ntest/ntot, ntest)
          # test inclusion prob = (sample inclusion prob) * (test inclusion prob)
          testInclusionProb = ps[testI] * testPs
          wsOrig = 1/testInclusionProb
          
          # adjust weights using Pareto smoothing
          ws = PSISsimple(wsOrig)
        }
        trainI = setdiff(1:ntot, testI)
        
        # Get train/test xs/ys ----
        trainLocs = sampleLocs[trainI,]
        testLocs = sampleLocs[testI,]
        trainYs = sampleYs[trainI,]
        testYs = sampleYs[testI,]
        
        # Get predictions ----
        # a = test
        # b = train
        # c = full domain
        muTrue = rGRFargsTrue$mu
        sigmaTrue = rGRFargsTrue$sigma
        rangeTrue = rGRFargsTrue$range
        smoothnessTrue = rGRFargsTrue$smoothness
        sigmaEpsTrue = sqrt(rGRFargsTrue$sigmaEpsSq)
        muFalse = rGRFargsFalse$mu
        sigmaFalse = rGRFargsFalse$sigma
        rangeFalse = rGRFargsFalse$range
        smoothnessFalse = rGRFargsFalse$smoothness
        sigmaEpsFalse = sqrt(rGRFargsFalse$sigmaEpsSq)
        
        # get distance matrices and true/false covariances
        distAB = rdist(testLocs, trainLocs)
        distCB = rdist(truth$locs, trainLocs)
        distAA = rdist(testLocs, testLocs)
        distBB = rdist(trainLocs, trainLocs)
        
        # rGRFargsFalse = list(mu=.1, sigma=1, 
        #                      cov.args=list(Covariance="Matern", range=0.2, smoothness=0.5), 
        #                      delta=5, sigmaEpsSq=0, nx=nx, ny=ny)
        
        SigmaABtrue = sigmaTrue^2 * stationary.cov(Covariance="Matern", aRange=rangeTrue, 
                                                   smoothness=smoothnessTrue, 
                                                   distMat=distAB)
        SigmaCBtrue = sigmaTrue^2 * stationary.cov(Covariance="Matern", aRange=rangeTrue, 
                                                   smoothness=smoothnessTrue, 
                                                   distMat=distCB)
        SigmaAAtrue = sigmaTrue^2 * stationary.cov(Covariance="Matern", aRange=rangeTrue, 
                                                   smoothness=smoothnessTrue, 
                                                   distMat=distAA) + 
          sigmaEpsTrue^2 * diag(nrow=nrow(testLocs))
        SigmaBBtrue = sigmaTrue^2 * stationary.cov(Covariance="Matern", aRange=rangeTrue, 
                                                   smoothness=smoothnessTrue, 
                                                   distMat=distBB) + 
          sigmaEpsTrue^2 * diag(nrow=nrow(trainLocs))
        
        SigmaABfalse = sigmaFalse^2 * stationary.cov(Covariance="Matern", aRange=rangeFalse, 
                                                     smoothness=smoothnessFalse, 
                                                     distMat=distAB)
        SigmaCBfalse = sigmaFalse^2 * stationary.cov(Covariance="Matern", aRange=rangeFalse, 
                                                     smoothness=smoothnessFalse, 
                                                     distMat=distCB)
        SigmaAAfalse = sigmaFalse^2 * stationary.cov(Covariance="Matern", aRange=rangeFalse, 
                                                     smoothness=smoothnessFalse, 
                                                     distMat=distAA) + 
          sigmaEpsFalse^2 * diag(nrow=nrow(testLocs))
        SigmaBBfalse = sigmaFalse^2 * stationary.cov(Covariance="Matern", aRange=rangeFalse, 
                                                     smoothness=smoothnessFalse, 
                                                     distMat=distBB) + 
          sigmaEpsFalse^2 * diag(nrow=nrow(trainLocs))
        
        # generate predictive means/variances under true/false covariances
        residsTrue = trainYs - muTrue
        residsFalse = trainYs - muFalse
        
        # muA|B = muA + SigmaAB SigmaBB^-1 (yB - muB)
        # muAcondBtrue = muTrue + SigmaABtrue %*% solve(SigmaBBtrue, residsTrue)
        # muAcondBfalse = muFalse + SigmaABfalse %*% solve(SigmaBBfalse, residsFalse)
        Utrue = chol(SigmaBBtrue)
        Ufalse = chol(SigmaBBfalse)
        predsTrue = muTrue + SigmaABtrue %*% t(Utrue) %*% backsolve(Utrue, residsTrue)
        predsFalse = muFalse + SigmaABfalse %*% t(Ufalse) %*% backsolve(Ufalse, residsFalse)
        predsFulltrue = muTrue + SigmaCBtrue %*% t(Utrue) %*% backsolve(Utrue, residsTrue)
        predsFullfalse = muFalse + SigmaCBfalse %*% t(Ufalse) %*% backsolve(Ufalse, residsFalse)
        
        # SigmaA|B = SigmaAA - SigmaAB SigmaBB^-1 SigmaBA
        #          = SigmaAA - SigmaAB (U' U)^-1 SigmaBA
        #          = SigmaAA - SigmaAB U^-1 (U')^-1 SigmaBA
        #          = SigmaAA - R' R
        RAtrue = forwardsolve(t(Utrue), t(SigmaABtrue))
        RAfalse = forwardsolve(t(Ufalse), t(SigmaABfalse))
        RCtrue = forwardsolve(t(Utrue), t(SigmaCBtrue))
        RCfalse = forwardsolve(t(Ufalse), t(SigmaCBfalse))
        
        # (R' R)_ii = (R')_i: R_:i
        varsTrue = myDiag(SigmaAAtrue) - myDiag(apply(RAtrue, 2, function(x) {sum(x^2)}))
        varsFalse = myDiag(SigmaAAfalse) - myDiag(apply(RAfalse, 2, function(x) {sum(x^2)}))
        varsFullTrue = myDiag(SigmaCCtrue) - myDiag(apply(RCtrue, 2, function(x) {sum(x^2)}))
        varsFullFalse = myDiag(SigmaCCfalse) - myDiag(apply(RCfalse, 2, function(x) {sum(x^2)}))
        
        # get scores (MSE, CRPS, coverage)
        scoresTrue = getScores(truth=testYs, est=predsTrue, var=varsTrue, weights=ws, 
                               significance=.95, getAverage=TRUE)
        scoresFalse = getScores(truth=testYs, est=predsFalse, var=varsFalse, weights=ws, 
                                significance=.95, getAverage=TRUE)
        scoresFullTrue = getScores(truth=truth$ys, est=predsFullTrue, var=varsFullTrue, 
                                   significance=.95, getAverage=TRUE)
        scoresFullFalse = getScores(truth=truth$ys, est=predsFullFalse, var=varsFullFalse, 
                                    significance=.95, getAverage=TRUE)
        
        # Save progress ----
        thisFileName = paste0(fileNameRoot, "_i", iter, 
                              "_popDistn", popDistn, 
                              "_distnEst", sampleRateEstMethod, 
                              "_vType", validationType, 
                              ".RData")
        
        save(scoresTrue, scoresFalse, scoresFullTrue, scoresFullFalse, file=thisFileName)
        
        # save scores under identical schemes if need be
        if((sampleScheme == "unif") && (popDistn == "unif")) {
          # "unif" + "sampleDistn" replaced by "unif" + "unif"
          thisFileName = paste0(fileNameRoot, "_i", iter, 
                                "_popDistn", "sampleDistn", 
                                "_distnEst", sampleRateEstMethod, 
                                "_vType", validationType, 
                                ".RData")
          
          save(scoresTrue, scoresFalse, scoresFullTrue, scoresFullFalse, file=thisFileName)
        } else if((sampleRateEstMethod == "truth") && (validationType == "basic")) {
          # "vcell" + "basic" and "kernel" + "basic" replaced by "truth" + "basic"
          thisFileName = paste0(fileNameRoot, "_i", iter, 
                                "_popDistn", popDistn, 
                                "_distnEst", "kernel", 
                                "_vType", validationType, 
                                ".RData")
          
          # "vcell" and "kernel" are not yet supported
          # save(scoresTrue, scoresFalse, scoresFullTrue, scoresFullFalse, file=thisFileName)
          
          thisFileName = paste0(fileNameRoot, "_i", iter, 
                                "_popDistn", popDistn, 
                                "_distnEst", "vcell", 
                                "_vType", validationType, 
                                ".RData")
          
          # save(scoresTrue, scoresFalse, scoresFullTrue, scoresFullFalse, file=thisFileName)
        } else if((sampleRateEstMethod == "truth") && (validationType == "gridded")) {
          # "vcell" + "gridded" and "kernel" + "gridded" replaced by "truth" + "gridded"
          thisFileName = paste0(fileNameRoot, "_i", iter, 
                                "_popDistn", popDistn, 
                                "_distnEst", "kernel", 
                                "_vType", validationType, 
                                ".RData")
          
          # "vcell" and "kernel" are not yet supported
          # save(scoresTrue, scoresFalse, scoresFullTrue, scoresFullFalse, file=thisFileName)
          
          thisFileName = paste0(fileNameRoot, "_i", iter, 
                                "_popDistn", popDistn, 
                                "_distnEst", "vcell", 
                                "_vType", validationType, 
                                ".RData")
          
          # save(scoresTrue, scoresFalse, scoresFullTrue, scoresFullFalse, file=thisFileName)
        } else if((sampleScheme == "unif") && (popDistn == "unif") && (validationType == "basic")) {
          # "unif" + "unif" + "thinning" replaced by "unif" + "unif" + "basic"
          thisFileName = paste0(fileNameRoot, "_i", iter, 
                                "_popDistn", popDistn, 
                                "_distnEst", sampleRateEstMethod, 
                                "_vType", "thinning", 
                                ".RData")
          
          save(scoresTrue, scoresFalse, scoresFullTrue, scoresFullFalse, file=thisFileName)
        } else if((popDistn == "sampleDistn") && (validationType == "basic")) {
          # "sampleDistn" + "thinning" replaced by "sampleDistn" + "basic"
          thisFileName = paste0(fileNameRoot, "_i", iter, 
                                "_popDistn", popDistn, 
                                "_distnEst", sampleRateEstMethod, 
                                "_vType", "thinning", 
                                ".RData")
          
          save(scoresTrue, scoresFalse, scoresFullTrue, scoresFullFalse, file=thisFileName)
        }
      }
    }
  }
  
  # Combine scores ----
  for(k in 1:length(popDistns)) {
    popDistn = popDistns[k]
    
    for(i in 1:length(sampleRateEstMethods)) {
      sampleRateEstMethod = sampleRateEstMethods[i]
      
      for(j in 1:length(validationTypes)) {
        validationType = validationTypes[j]
        
        # load file
        thisFileName = paste0(fileNameRoot, "_i", iter, 
                              "_distnEst", sampleRateEstMethod, 
                              "_vType", validationType, 
                              ".RData")
        out = load(thisFileName)
        # scoresTrue, scoresFalse, scoresFullTrue, scoresFullFalse
        
        if(j == 1) {
          tempScoreTabTrue = scoresTrue
          tempScoreTabFalse = scoresFalse
          tempScoreTabFullTrue = scoresFullTrue
          tempScoreTabFullFalse = scoresFullFalse
        } else {
          tempScoreTabTrue = rbind(tempScoreTabTrue, scoresTrue)
          tempScoreTabFalse = rbind(tempScoreTabFalse, scoresFalse)
          tempScoreTabFullTrue = rbind(tempScoreTabFullTrue, scoresFullTrue)
          tempScoreTabFullFalse = rbind(tempScoreTabFullFalse, scoresFullFalse)
        }
      }
      
      # label scores
      valType = validationTypes
      tempScoreTabTrue = cbind(sampleScheme=sampleScheme, 
                               popDistn=popDistn, 
                               sampleSchemeEst=sampleRateEstMethod, 
                               valType, tempScoreTabTrue)
      tempScoreTabFalse = cbind(sampleScheme=sampleScheme, 
                                popDistn=popDistn, 
                                sampleSchemeEst=sampleRateEstMethod, 
                                valType, tempScoreTabFalse)
      tempScoreTabFullTrue = cbind(sampleScheme=sampleScheme, 
                                   popDistn=popDistn, 
                                   sampleSchemeEst=sampleRateEstMethod, 
                                   valType, tempScoreTabFullTrue)
      tempScoreTabFullFalse = cbind(sampleScheme=sampleScheme, 
                                    popDistn=popDistn, 
                                    sampleSchemeEst=sampleRateEstMethod, 
                                    valType, tempScoreTabFullFalse)
      
      # concatenate scores to full score table
      if(i == 1) {
        scoreTabTrue = tempScoreTabTrue
        scoreTabFalse = tempScoreTabFalse
        scoreTabFullTrue = tempScoreTabFullTrue
        scoreTabFullFalse = tempScoreTabFullFalse
      } else {
        scoreTabTrue = rbind(scoreTabTrue, tempScoreTabTrue)
        scoreTabFalse = rbind(scoreTabFalse, tempScoreTabFalse)
        scoreTabFullTrue = rbind(scoreTabFullTrue, tempScoreTabFullTrue)
        scoreTabFullFalse = rbind(scoreTabFullFalse, tempScoreTabFullFalse)
      }
    }
  }
  
  # Save final results ----
  fileName = paste0(fileNameRoot, "_i", iter, ".RData")
  
  save(scoreTabTrue, scoreTabFalse, scoreTabFullTrue, scoreTabFullFalse, file=fileName)
  
  list(scoreTabTrue, scoreTabFalse, scoreTabFullTrue, scoreTabFullFalse)
}


