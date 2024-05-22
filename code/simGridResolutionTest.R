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
#   10.   Calculate importance weighted LOO-CV MSE, or LOOIS-CV MSE
#   11.   Calculate actual predictive MSE
#   12.   Save result
#   13. Combine and save all results
# Note that the below function is a single iteration of the for loop 
# for easy parallelization.
griddedResTestIter = function(rGRFargsTruth=NULL, rGRFargsSample=NULL, 
                              n=50, gridNs=2^(1:6), iter=1, rGRFargsWrong=NULL, 
                              nx=100, ny=100, sigmaEpsSq=0, allSeeds=123, 
                              unif=FALSE, preferential=FALSE, alpha=0, beta=1, 
                              printProgress=FALSE) {
  
  print(paste0("iteration ", iter, "/", length(allSeeds)))
  
  if(!is.null(allSeeds)) {
    set.seed(allSeeds[iter])
  }
  
  # set default GRF parameters
  if(is.null(rGRFargsTruth)) {
    rGRFargsTruth = list(mu=0, sigma=1, 
                         cov.args=list(Covariance="Matern", range=0.2, smoothness=1.0), 
                         delta=3, sigmaEpsSq=0, nx=nx, ny=ny)
  }
  if(is.null(rGRFargsWrong)) {
    rGRFargsWrong = rGRFargsTruth
    
    if(n==50) {
      rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt(1.30)
    } else if(n==500){
      rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt(1.10)
    }
  }
  if(is.null(rGRFargsSample)) {
    if(!unif) {
      rGRFargsSample = list(mu=0, sigma=1, 
                            cov.args=list(Covariance="Matern", range=0.2, smoothness=1.0), 
                            delta=3, sigmaEpsSq=0, nx=nx, ny=ny)
    } else {
      rGRFargsSample = list(mu=0, sigma=0, 
                            cov.args=list(Covariance="Matern", range=0.2, smoothness=1.0), 
                            delta=3, sigmaEpsSq=0, nx=nx, ny=ny)
    }
    
  }
  
  # 2. Simulate 1 GRF for responses ----
  #    on unit square
  truthGRF = do.call("rGRF", rGRFargsTruth)
  truth = truthGRF$truth
  locs = truthGRF$locs
  
  # 3. Simulate sample distribution ----
  #    as 1 GRF on unit square
  if(!preferential) {
    sampleGRF = do.call("rGRF", rGRFargsSample)
    sampleRates = exp(sampleGRF$truth)
    sampleProbs = sampleRates * (1/sum(sampleRates))
  } else {
    sampleGRF = truthGRF
    sampleRates = exp(alpha + beta * sampleGRF$truth)
    sampleProbs = sampleRates * (1/sum(sampleRates))
  }
  
  # 4. Simulate observations ----
  #    in domain based on the sampling distn and the GRF 
  #    for responses (don't include nugget?)
  gridResX = 1/nx
  gridResY = 1/ny
  sampleI = sample(1:nrow(locs), n, prob=sampleProbs, replace=TRUE)
  xs = locs[sampleI,] + cbind(runif(n, max=gridResX)-gridResX/2, 
                              runif(n, max=gridResY)-gridResY/2)
  ys = truth[sampleI] + rnorm(n, sd=sqrt(sigmaEpsSq))
  
  # 5. Get covariance matrices ----
  #    for sample and crossCov to all locs
  
  SigmaSample = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsTruth$cov.args$range, 
                               smoothness=rGRFargsTruth$cov.args$smoothness) * rGRFargsTruth$sigma^2
  SigmaGridToSample = stationary.cov(locs, xs, Covariance="Matern", aRange=rGRFargsTruth$cov.args$range, 
                                     smoothness=rGRFargsTruth$cov.args$smoothness) * rGRFargsTruth$sigma^2
  
  # 6. for i in 1:number block resolutions: ----
  
  # determine grid cells associated with observations
  getCellI = function(thisLoc) {
    xI = match(FALSE, thisLoc[1] > highs)
    yI = match(FALSE, thisLoc[2] > highs)
    
    indMat[yI, xI]
  }
  
  griddedCVs = numeric(length(gridNs))
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
    cellIsSample = apply(xs, 1, getCellI)
    uniqueCellIs = sort(unique(cellIsSample))
    blockCVs = numeric(length(uniqueCellIs))
    
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
      
      # predict test data
      SigmaAB = matrix(SigmaSample[isTest, !isTest], nrow=sum(isTest))
      SigmaBB = SigmaSample[!isTest, !isTest]
      SigmaAA = SigmaSample[isTest, isTest]
      
      if((sum(isTest) > 0) && (sum(!isTest) > 0)) {
        condDistn = condMeanMVN(SigmaAA=SigmaAA, SigmaAB=SigmaAB, SigmaBB=SigmaBB, 
                                ysB=trainYs, getFullCov=FALSE, getCondVar=FALSE)
        muAcondB = condDistn$muAcondB
        
        # calculate MSE for the grid cell
        blockCVs[j] = mean((testYs - muAcondB)^2)
      } else {
        blockCVs[j] = NA
      }
    }
    
    # average over MSEs of each grid cell
    griddedCVs[i] = mean(blockCVs, na.rm=TRUE)
  }
  
  # 9. get LOO-CV MSE ----
  LOOCVs = numeric(n)
  for(i in 1:n) {
    if(printProgress) {
      print(paste0("Leaving out obs ", i, "/", n))
    }
    
    # separate data in test/train
    isTest = (1:n) == i
    testYs = ys[isTest]
    trainYs = ys[!isTest]
    
    # predict test data
    SigmaAB = matrix(SigmaSample[isTest, !isTest], nrow=1)
    SigmaBB = SigmaSample[!isTest, !isTest]
    SigmaAA = SigmaSample[isTest, isTest]
    
    condDistn = condMeanMVN(SigmaAA=SigmaAA, SigmaAB=SigmaAB, SigmaBB=SigmaBB, 
                            ysB=trainYs, getFullCov=FALSE, getCondVar=FALSE)
    muAcondB = condDistn$muAcondB
    
    # calculate MSE for the grid cell
    LOOCVs[i] = mean((testYs - muAcondB)^2)
  }
  LOOCV = mean(LOOCVs)
  
  # 10. Calculate LOOIS-CV MSE ----
  cellArea = 1/(nx*ny)
  LOOISrates = (sampleProbs[sampleI]/sum(sampleProbs)) / cellArea
  # LOOISCV = weighted.mean(LOOCVs, w=cellArea/LOOISrates)
  LOOISCV = sum(LOOCVs * (1/LOOISrates))/n
  
  # 11. Calculate LOOVC-CV MSE ----
  vcellInfo = getVCellAreas(xs, domainPoly = rbind(c(0, 0), 
                                                   c(0, 1), 
                                                   c(1, 1), 
                                                   c(1, 0), 
                                                   c(0, 0)))
  vcellArea = vcellInfo$area
  rateEsts = 1/vcellArea
  rateEsts = rateEsts/n # divide by n to get rate for a single pt draw instead of n pts
  LOOVCCV = sum(LOOCVs * (1/rateEsts))/n
  
  # 12. Calculate true MSE ----
  condDistn = condMeanMVN(SigmaAA=rep(1, nrow(locs)), SigmaAB=SigmaGridToSample, SigmaBB=SigmaSample, 
                          ysB=ys, getFullCov=FALSE, getCondVar=FALSE)
  muAcondB = condDistn$muAcondB
  
  # calculate MSE for the grid cell
  trueMSE = mean((truth - muAcondB)^2)
  
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
  
  # 13. Save result ----
  # browser()
  unifText = ifelse(unif, "_unif", "")
  prefText = ifelse(preferential, paste0("_prefA", alpha, "B", beta), "")
  save(trueMSE, LOOCVs, LOOCV, griddedCVs, gridNs, 
       iter, rGRFargsTruth, rGRFargsSample, 
       n, nx, ny, sigmaEpsSq, allSeeds, 
       alpha, beta, 
       file=paste0("savedOutput/griddedCVtest/n1", n, "_n2", n2, "_iter", iter, unifText, prefText, ".RData"))
  
  list(trueMSE=trueMSE, LOOCVs=LOOCVs, LOOCV=LOOCV, 
       LOOISCV=LOOISCV, LOOVCCV=LOOVCCV, griddedCVs=griddedCVs, 
       gridNs=gridNs, iter=iter, rGRFargsTruth=rGRFargsTruth, 
       rGRFargsSample=rGRFargsSample, 
       n=n, nx=nx, ny=ny, sigmaEpsSq=sigmaEpsSq, 
       allSeeds=allSeeds)
}