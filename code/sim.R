
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
  out = c(circulantEmbedding(obj))
  truth = mu + out*sigma
  obs = truth + rnorm(length(out), sd=sqrt(sigmaEpsSq))
  
  list(locs=locs, truth=truth, obs=obs, xs=xs, ys=ys, nx=nx, ny=ny, 
       cov.args=cov.args, mu=mu, sigma=sigma, pixelArea=1/(nx*ny))
}

# simulate 2 correlated GRF on unit square. See ?circulantEmbedding from fields
r2GRFs = function(nx=100, ny=100, mu1=0, mu2=mu1, sigma1=1, sigma2=sigma1, rho=-.8, 
                  cov.args=list(Covariance="Matern", range=0.3, smoothness=1.0), 
                  delta=5, sigmaEpsSq1=0, sigmaEpsSq2=sigmaEpsSq1) {
  
  Z1 = rGRF(nx=nx, ny=ny, mu=0, sigma=1, cov.args=cov.args, delta=delta, sigmaEpsSq=0)
  Z2 = rGRF(nx=nx, ny=ny, mu=0, sigma=1, cov.args=cov.args, delta=delta, sigmaEpsSq=0)
  
  Y1 = Z1
  Y2 = Z1
  Y1$truth = sigma1 * (rho * Z1$truth + sqrt(1-rho^2) * Z2$truth) + mu1
  Y2$truth = sigma2 * Z1$truth + mu2
  Y1$obs = Y1$truth
  Y2$obs = Y2$truth
  if(sigmaEpsSq1 > 0) {
    Y1$obs = Y1$obs + rnorm(length(Y1$obs), 0, sqrt(sigmaEpsSq1))
  }
  if(sigmaEpsSq2 > 0) {
    Y2$obs = Y2$obs + rnorm(length(Y2$obs), 0, sqrt(sigmaEpsSq2))
  }
  
  list(Y1=Y1, Y2=Y2)
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
  RA = forwardsolve(t(U), t(SigmaAB))
  Rpt = forwardsolve(t(U), ysB - muB)
  # Lpt = backsolve(U, Rpt)
  # muAcondB = muA + SigmaAB %*% t(U) %*% backsolve(U, ysB - muB)
  # muAcondB = muA + SigmaAB %*% Lpt
  muAcondB = muA + t(RA) %*% Rpt
  
  # now get predictive/conditional variance if requested by user:
  if(getCondVar) {
    # SigmaA|B = SigmaAA - SigmaAB SigmaBB^-1 SigmaBA
    #          = SigmaAA - SigmaAB (U' U)^-1 SigmaBA
    #          = SigmaAA - SigmaAB U^-1 (U')^-1 SigmaBA
    #          = SigmaAA - R' R
    
    
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
griddedResTestAll = function(rGRFargsTruth=NULL, rGRFargsSample=NULL, rGRFargsWrong=NULL, 
                             n=50, gridNs=2^(1:6), niter=100, seed=123, 
                             twoDatasets=FALSE, 
                             nx=100, ny=100, sigmaEpsSq=ifelse(twoDatasets, .1^2, 0), sigmaEpsSq2=9*sigmaEpsSq, 
                             rho=-.8, regenResults=TRUE, 
                             printProgress=FALSE, relTicks1=NULL, relTickLabs1=NULL, 
                             relTicks2=NULL, relTickLabs2=NULL, unif=FALSE, 
                             preferential=FALSE, alpha=0, beta=1) {
  set.seed(seed)
  seeds = sample(1:100000, niter)
  
  if(twoDatasets) {
    unif=FALSE
    preferential=FALSE
  }
  unifText = ifelse(unif, "_unif", "")
  prefText = ifelse(preferential, paste0("_prefA", alpha, "B", beta), "")
  twoDatText = ifelse(twoDatasets, paste0("_2dat_rho", rho), "")
  unifTitleText = ifelse(unif, ", unif", "")
  prefTitleText = ifelse(preferential, paste0("pref, alpha=", alpha, ", beta=", beta), "")
  twoDatTitleText = ifelse(twoDatasets, paste0(", 2 datasets, rho=", rho), "")
  
  # generate results
  if(regenResults) {
    if(!twoDatasets) {
      totTime = system.time(results <- lapply(1:niter, griddedResTestIter, 
                                              rGRFargsTruth=rGRFargsTruth, 
                                              rGRFargsSample=rGRFargsSample, 
                                              n=n, gridNs=gridNs, allSeeds=seeds, 
                                              nx=nx, ny=ny, sigmaEpsSq=sigmaEpsSq, 
                                              printProgress=printProgress, unif=unif, 
                                              preferential=preferential, alpha=alpha, 
                                              beta=beta))
    } else {
      if(is.null(rGRFargsSample)) {
        rGRFargsSample = list(mu1=0, mu2=0, sigma1=1, sigma2=1, rho=rho, 
                              cov.args=list(Covariance="Matern", range=0.2, smoothness=1.0), 
                              delta=3, sigmaEpsSq1=sigmaEpsSq, sigmaEpsSq2=sigmaEpsSq2, nx=nx, ny=ny)
      }
      
      totTime = system.time(results <- lapply(1:niter, griddedResTestIter2Datasets, 
                                              rGRFargsTruth=rGRFargsTruth, 
                                              rGRFargsSample=rGRFargsSample, 
                                              n1=n, n2=n, gridNs=gridNs, allSeeds=seeds, 
                                              nx=nx, ny=ny, rho=rho, 
                                              sigmaEpsSq1=sigmaEpsSq, sigmaEpsSq2=sigmaEpsSq2, 
                                              alpha=alpha, beta=beta, 
                                              printProgress=printProgress))
    }
  } else {
    # load results
    results = list()
    for(i in 1:niter) {
      if(!twoDatasets) {
        out = load(paste0("savedOutput/griddedCVtest/n", n, "_iter", i, unifText, prefText, twoDatText, ".RData"))
        thisList = list(trueMSE=trueMSE, LOOCVs=LOOCVs, LOOCV=LOOCV, 
                        griddedCVs=griddedCVs, gridNs=gridNs, iter=iter, 
                        rGRFargsTruth=rGRFargsTruth, rGRFargsSample=rGRFargsSample, rGRFargsWrong=rGRFargsWrong, 
                        n=n, nx=nx, ny=ny, sigmaEpsSq=sigmaEpsSq, allSeeds=allSeeds, 
                        alpha=alpha, beta=beta)
      } else {
        out = load(paste0("savedOutput/griddedCVtest2Datasets/n1", n1, "_n2", n2, "_rho", rho, "_iter", iter, ".RData"))
        thisList = list(trueMSE=trueMSE, LOOCVs=LOOCVs, LOOCV=LOOCV, 
                        griddedCVs=griddedCVs, gridNs=gridNs, iter=iter, 
                        rGRFargsTruth=rGRFargsTruth, rGRFargsSample=rGRFargsSample, rGRFargsWrong=rGRFargsWrong, 
                        rho=rho, n1=n1, n2=n2, sigmaEpsSq1=sigmaEpsSq1, sigmaEpsSq2=sigmaEpsSq2, 
                        alpha=alpha, beta=beta, nx=nx, ny=ny, allSeeds=allSeeds)
      }
      
      results = c(results, list(thisList))
    }
  }
  
  # concatenate results
  getName = function(x, thisName) {
    x[[thisName]]
  }
  # list(trueMSE=trueMSE, LOOCVs=LOOCVs, LOOCV=LOOCV, 
  #      griddedCVs=griddedCVs, gridNs=gridNs, 
  #      iter=iter, seed=seed, rGRFargsTruth=rGRFargsTruth, 
  #      rGRFargsSample=rGRFargsSample, 
  #      n=n, nx=nx, ny=ny, sigmaEpsSq=sigmaEpsSq, seed=seed)
  trueMSEs = sapply(results, getName, thisName="trueMSE")
  LOOCVs = sapply(results, getName, thisName="LOOCV")
  LOOISCVs = sapply(results, getName, thisName="LOOISCV")
  LOOVCCVs = sapply(results, getName, thisName="LOOVCCV")
  griddedCVs = sapply(results, getName, thisName="griddedCVs")
  
  LOOCVsWrong = sapply(results, getName, thisName="LOOCVWrong")
  LOOISCVsWrong = sapply(results, getName, thisName="LOOISCVWrong")
  LOOVCCVsWrong = sapply(results, getName, thisName="LOOVCCVWrong")
  griddedCVsWrong = sapply(results, getName, thisName="griddedCVsWrong")
  
  # calculate error
  LOOCVerrs = LOOCVs - trueMSEs
  LOOISCVerrs = LOOISCVs - trueMSEs
  LOOVCCVerrs = LOOVCCVs - trueMSEs
  griddedCVerrs = sweep(griddedCVs, 2, trueMSEs)
  griddedCVmeanErr = rowMeans(griddedCVerrs)
  
  # calculate percent error
  LOOCVpctErrs = 100 * LOOCVerrs/trueMSEs
  LOOISCVpctErrs = 100 * LOOISCVerrs/trueMSEs
  LOOVCCVpctErrs = 100 * LOOVCCVerrs/trueMSEs
  griddedCVpctErrs = sweep(griddedCVerrs, 2, 100/trueMSEs, "*")
  griddedCVmeanPctErr = rowMeans(griddedCVpctErrs)
  
  # calculate relative error
  LOOCVrelErrs = LOOCVs/trueMSEs
  LOOISCVrelErrs = LOOISCVs/trueMSEs
  LOOVCCVrelErrs = LOOVCCVs/trueMSEs
  griddedCVrelErrs = sweep(griddedCVs, 2, trueMSEs, "/")
  griddedCVmeanRelErr = rowMeans(griddedCVrelErrs)
  
  # calculate proportion of time correct model is selected
  LOOCVprop = mean(LOOCVs < LOOCVsWrong)
  LOOISCVprop = mean(LOOISCVs < LOOISCVsWrong)
  LOOVCCVprop = mean(LOOVCCVs < LOOVCCVsWrong)
  griddedCVprop = rowMeans(griddedCVs < griddedCVsWrong)
  
  LOOCVpropMOE = qnorm(.975) * sqrt(LOOCVprop*(1-LOOCVprop) / niter)
  LOOISCVpropMOE = qnorm(.975) * sqrt(LOOISCVprop*(1-LOOISCVprop) / niter)
  LOOVCCVpropMOE = qnorm(.975) * sqrt(LOOVCCVprop*(1-LOOVCCVprop) / niter)
  griddedCVpropMOE = qnorm(.975) * sqrt(griddedCVprop*(1-griddedCVprop) / niter)
  
  LOOCVpropHigh = LOOCVprop + LOOCVpropMOE
  LOOISCVpropHigh = LOOISCVprop + LOOISCVpropMOE
  LOOVCCVpropHigh = LOOVCCVprop + LOOVCCVpropMOE
  griddedCVpropHigh = griddedCVprop + griddedCVpropMOE
  
  LOOCVpropLow = LOOCVprop - LOOCVpropMOE
  LOOISCVpropLow = LOOISCVprop - LOOISCVpropMOE
  LOOVCCVpropLow = LOOVCCVprop - LOOVCCVpropMOE
  griddedCVpropLow = griddedCVprop - griddedCVpropMOE
  
  figureFolder = ifelse(twoDatasets, "figures/twoDatasetsTest/", "figures/gridTest/")
  
  # plot results ----
  # proportion of time selecting right model
  pdf(paste0(figureFolder, "selectProb_n", n, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
  ylim = range(c(LOOCVpropHigh, LOOISCVpropHigh, LOOVCCVpropHigh, griddedCVpropHigh, 
                 LOOCVpropLow, LOOISCVpropLow, LOOVCCVpropLow, griddedCVpropLow))
  plot(gridNs, griddedCVprops, type="n", log="x", axes=FALSE, 
       ylim=ylim, xlab="Blocks per side", 
       ylab="Probability", main=paste0("Gridded CV vs resolution (n=", n, unifTitleText, prefTitleText, twoDatTitleText, ")"))
  axis(side=1, at=gridNs)
  axis(side=2)
  box()
  lines(gridNs, griddedCVprop, type="o", pch=19, col="blue")
  arrows(gridNs, griddedCVprop, gridNs, griddedCVpropHigh, angle=90)
  arrows(gridNs, griddedCVprop, gridNs, griddedCVpropLow, angle=90)
  abline(h=LOOCVprop, lty=1, col="purple")
  abline(h=LOOCVpropHigh, lty=2, col="purple")
  abline(h=LOOCVpropLow, lty=2, col="purple")
  abline(h=LOOISCVprop, lty=1, col="orange")
  abline(h=LOOISCVpropHigh, lty=2, col="orange")
  abline(h=LOOISCVpropLow, lty=2, col="orange")
  abline(h=LOOVCCVprop, lty=1, col="brown")
  abline(h=LOOVCCVpropHigh, lty=2, col="brown")
  abline(h=LOOVCCVpropLow, lty=2, col="brown")
  abline(h=0, lty=2, col="green")
  legend("topright", c("Gridded", "LOO", "LOOVC", "LOOIS", "Truth"), col=c("blue", "purple", "brown", "orange", "green"), 
         pch=c(19, NA, NA, NA, NA), lty=1)
  dev.off()
  
  # absolute bias
  pdf(paste0(figureFolder, "griddedBias_n", n, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
  plot(gridNs, griddedCVmeanErr, type="n", log="x", axes=FALSE, 
       ylim=c(min(c(griddedCVmeanErr, mean(LOOCVerrs))), max(griddedCVmeanErr)), 
       xlab="Blocks per side", 
       ylab="MSE Bias", main=paste0("Gridded CV vs resolution (n=", n, unifTitleText, prefTitleText, twoDatTitleText, ")"))
  axis(side=1, at=gridNs)
  axis(side=2)
  box()
  lines(gridNs, griddedCVmeanErr, type="o", pch=19, col="blue")
  abline(h=mean(LOOCVerrs), lty=2, col="purple")
  abline(h=mean(LOOISCVerrs), lty=2, col="orange")
  abline(h=mean(LOOVCCVerrs), lty=2, col="brown")
  abline(h=0, lty=2, col="green")
  legend("topright", c("Gridded", "LOO", "LOOVC", "LOOIS", "Truth"), col=c("blue", "purple", "brown", "orange", "green"), 
         pch=c(19, NA, NA, NA, NA), lty=c(1, 2, 2, 2, 2))
  dev.off()
  
  pdf(paste0(figureFolder, "griddedBiasExtrap_n", n, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
  plot(gridNs, griddedCVmeanErr, type="n", log="x", axes=FALSE, 
       ylim=c(min(c(griddedCVmeanErr, mean(LOOCVerrs))), 1+sigmaEpsSq-mean(trueMSEs)), 
       xlab="Blocks per side", 
       ylab="MSE Bias", main=paste0("Gridded CV vs resolution (n=", n, unifTitleText, prefTitleText, twoDatTitleText, ")"))
  axis(side=1, at=gridNs)
  axis(side=2)
  box()
  lines(gridNs, griddedCVmeanErr, type="o", pch=19, col="blue")
  abline(h=mean(LOOCVerrs), lty=2, col="purple")
  abline(h=mean(LOOISCVerrs), lty=2, col="orange")
  abline(h=mean(LOOVCCVerrs), lty=2, col="brown")
  abline(h=0, lty=2, col="green")
  abline(h=1+sigmaEpsSq-mean(trueMSEs), lty=2, col="red")
  legend("right", c("Gridded", "LOO", "LOOIS", "LOOVC", "Extrapolation", "Interpolation"), 
         col=c("blue", "purple", "orange", "brown", "red", "green"), 
         pch=c(19, NA, NA, NA, NA, NA), lty=c(1, 2, 2, 2, 2, 2))
  dev.off()
  
  # then relative/percent bias
  pdf(paste0(figureFolder, "griddedPctBias_n", n, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
  plot(gridNs, griddedCVmeanPctErr, type="n", log="x", axes=FALSE, 
       ylim=c(min(c(griddedCVmeanPctErr, mean(LOOCVpctErrs))), max(griddedCVmeanPctErr)), 
       xlab="Blocks per side", 
       ylab="MSE Relative Bias (%)", main=paste0("Gridded CV vs resolution (n=", n, unifTitleText, prefTitleText, twoDatTitleText, ")"))
  axis(side=1, at=gridNs)
  axis(side=2)
  box()
  lines(gridNs, griddedCVmeanPctErr, type="o", pch=19, col="blue")
  abline(h=mean(LOOCVpctErrs), lty=2, col="purple")
  abline(h=mean(LOOISCVpctErrs), lty=2, col="orange")
  abline(h=mean(LOOVCCVpctErrs), lty=2, col="brown")
  abline(h=0, lty=2, col="green")
  legend("topright", c("Gridded", "LOO", "LOOIS", "LOOVC", "Truth"), col=c("blue", "purple", "orange", "brown", "green"), 
         pch=c(19, NA, NA, NA, NA), lty=c(1, 2, 2, 2, 2))
  dev.off()
  
  pdf(paste0(figureFolder, "griddedPctBiasExtrap_n", n, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
  plot(gridNs, griddedCVmeanPctErr, type="n", log="x", axes=FALSE, 
       ylim=c(min(c(griddedCVmeanPctErr, mean(LOOCVpctErrs))), 100*mean((1+sigmaEpsSq-trueMSEs)/trueMSEs)), 
       xlab="Blocks per side", 
       ylab="MSE Relative Bias (%)", main=paste0("Gridded CV vs resolution (n=", n, unifTitleText, prefTitleText, twoDatTitleText, ")"))
  axis(side=1, at=gridNs)
  axis(side=2)
  box()
  lines(gridNs, griddedCVmeanPctErr, type="o", pch=19, col="blue")
  abline(h=mean(LOOCVpctErrs), lty=2, col="purple")
  abline(h=mean(LOOISCVpctErrs), lty=2, col="orange")
  abline(h=mean(LOOVCCVpctErrs), lty=2, col="brown")
  abline(h=0, lty=2, col="green")
  abline(h=100*mean((1+sigmaEpsSq-trueMSEs)/trueMSEs), lty=2, col="red")
  legend("right", c("Gridded", "LOO", "LOOIS", "LOOVC", "Extrapolation", "Interpolation"), 
         col=c("blue", "purple", "orange", "brown", "red", "green"), 
         pch=c(19, NA, NA, NA, NA, NA), lty=c(1, 2, 2, 2, 2, 2))
  dev.off()
  
  # then relative bias
  if(is.null(relTicks1)) {
    if(n == 50) {
      if(!unif && !twoDatasets) {
        relTicks1 = c(.67, 1, 1.5, 2, 2.5, 3)
        relTickLabs1 = c("0.67", "1", "1.5", "2", "2.5", "3")
        relTicks2 = c(.67, 1, 2, 3, 5)
        relTickLabs2 = c("0.67", "1", "2", "3", "5")
      } else if(!twoDatasets) {
        relTicks1 = c(.67, 1, 2, 3, 4)
        relTickLabs1 = c("0.67", "1", "2", "3", "4")
        relTicks2 = c(.67, 1, 2, 4, 6, 8)
        relTickLabs2 = c("0.67", "1", "2", "4", "6", "8")
      } else {
        if(rho == -.8) {
          relTicks1 = seq(.7, 1.5, by=.1)
          relTickLabs1 = as.character(relTicks1)
          relTicks2 = c(.67, 1, 2, 3, 5)
          relTickLabs2 = c("0.67", "1", "2", "3", "5")
        } else if(rho == 0) {
          relTicks1 = seq(.7, 3, by=.1)
          relTickLabs1 = as.character(relTicks1)
          relTicks2 = c(.67, 1, 2, 3, 5)
          relTickLabs2 = c("0.67", "1", "2", "3", "5")
        }
        
      }
    } else if(n == 500) {
      if(!unif && !twoDatasets) {
        relTicks = c(.6, 1, 2, 5, 10, 20, 50, 100)
        relTickLabs = c("0.6", "1", "2", "5", "10", "20", "50", "100")
        relTicks1 = relTicks2 = relTicks
        relTicksLab1 = relTickLabs2 = relTickLabs
      } else if(!twoDatasets) {
        browser()
      } else {
        relTicks1 = seq(.7, 3, by=.1)
        relTickLabs1 = as.character(relTicks1)
        relTicks2 = c(.67, 1, 2, 3, 5, 10, 20, 30, 40, 50, 60)
        relTickLabs2 = c("0.67", "1", "2", "3", "5", "10", "20", "30", "40", "50", "60")
      }
    } else {
      relTicks = c(.67, 1, 2, 5, 10, 20, 50, 100)
      relTickLabs = c("0.67", "1", "2", "5", "10", "20", "50", "100")
      relTicks1 = relTicks2 = relTicks
      relTicksLab1 = relTickLabs2 = relTickLabs
    }
  }
  
  browser()
  
  pdf(paste0(figureFolder, "griddedRelBias_n", n, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
  ylim = c(min(c(griddedCVmeanRelErr, mean(LOOCVrelErrs))), max(griddedCVmeanRelErr))
  plot(gridNs, griddedCVmeanRelErr, type="n", log="xy", axes=FALSE, 
       ylim=ylim, 
       xlab="Blocks per side", 
       ylab="MSE Relative Bias", main=paste0("Gridded CV vs resolution (n=", n, unifTitleText, prefTitleText, twoDatTitleText, ")"))
  axis(side=1, at=gridNs)
  axis(side=2, at=relTicks1, labels=relTickLabs1)
  # axis(side=2)
  box()
  lines(gridNs, griddedCVmeanRelErr, type="o", pch=19, col="blue")
  abline(h=mean(LOOCVrelErrs), lty=2, col="purple")
  abline(h=mean(LOOISCVrelErrs), lty=2, col="orange")
  abline(h=mean(LOOVCCVrelErrs), lty=2, col="brown")
  abline(h=1, lty=2, col="green")
  legend("topright", c("Gridded", "LOO", "LOOIS", "LOOVC", "Truth"), 
         col=c("blue", "purple", "orange", "brown", "green"), 
         pch=c(19, NA, NA, NA, NA), lty=c(1, 2, 2, 2, 2))
  dev.off()
  
  pdf(paste0(figureFolder, "griddedRelBiasExtrap_n", n, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
  plot(gridNs, griddedCVmeanRelErr, type="n", log="xy", axes=FALSE, 
       ylim=c(min(c(griddedCVmeanRelErr, mean(LOOCVrelErrs))), mean((1+sigmaEpsSq)/trueMSEs)), 
       xlab="Blocks per side", 
       ylab="MSE Relative Bias", main=paste0("Gridded CV vs resolution (n=", n, unifTitleText, prefTitleText, twoDatTitleText, ")"))
  axis(side=1, at=gridNs)
  axis(side=2, at=relTicks2, labels=relTickLabs2)
  box()
  lines(gridNs, griddedCVmeanRelErr, type="o", pch=19, col="blue")
  abline(h=mean(LOOCVrelErrs), lty=2, col="purple")
  abline(h=mean(LOOISCVrelErrs), lty=2, col="orange")
  abline(h=mean(LOOVCCVrelErrs), lty=2, col="brown")
  abline(h=1, lty=2, col="green")
  abline(h=mean((1+sigmaEpsSq)/trueMSEs), lty=2, col="red")
  legend("right", c("Gridded", "LOO", "LOOIS", "LOOVC", "Extrapolation", "Interpolation"), 
         col=c("blue", "purple", "orange", "brown", "red", "green"), 
         pch=c(19, NA, NA, NA, NA, NA), lty=c(1, 2, 2, 2, 2, 2))
  dev.off()
  
  # boxplots ----
  pdf(paste0(figureFolder, "CVMSE_n", n, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
  dat = data.frame(
    Method=rep(c("LOO", "LOO-IS", "LOO-VC"), each=n), 
    sqResids=c((LOOCVs - trueMSEs)^2, 
               (LOOISCVs - trueMSEs)^2, 
               (LOOVCCVs - trueMSEs)^2))
  boxplot(sqResids~Method, data=dat, log="y", col="skyblue", ylab="Estimator Sq. Err.")
  dev.off()
  
  pdf(paste0(figureFolder, "CVMRE_n", n, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
  dat = data.frame(
    Method=rep(c("LOO", "LOO-IS", "LOO-VC"), each=n), 
    sqResids=c(LOOCVs/trueMSEs, 
               LOOISCVs/trueMSEs, 
               LOOVCCVs/trueMSEs))
  boxplot(sqResids~Method, data=dat, log="y", col="skyblue", ylab="Estimator Rel. Err.")
  dev.off()
  
  print(paste0("MSE of LOO: ", mean((LOOCVs - trueMSEs)^2)))
  print(paste0("MSE of LOOIS: ", mean((LOOISCVs - trueMSEs)^2)))
  print(paste0("MSE of LOOVC: ", mean((LOOVCCVs - trueMSEs)^2)))
  
  print(paste0("MAE of LOO: ", mean(abs(LOOCVs - trueMSEs))))
  print(paste0("MAE of LOOIS: ", mean(abs(LOOISCVs - trueMSEs))))
  print(paste0("MAE of LOOVC: ", mean(abs(LOOVCCVs - trueMSEs))))
  
  print(paste0("MRE of LOO: ", mean(LOOCVs/trueMSEs)))
  print(paste0("MRE of LOOIS: ", mean(LOOISCVs/trueMSEs)))
  print(paste0("MRE of LOOVC: ", mean(LOOVCCVs/trueMSEs)))
  
  print(paste0("MAPE of LOO: ", 100*mean(abs(LOOCVs/trueMSEs-1))))
  print(paste0("MAPE of LOOIS: ", 100*mean(abs(LOOISCVs/trueMSEs-1))))
  print(paste0("MAPE of LOOVC: ", 100*mean(abs(LOOVCCVs/trueMSEs-1))))
  
  print(paste0("MPE of LOO: ", 100*mean(LOOCVs/trueMSEs-1)))
  print(paste0("MPE of LOOIS: ", 100*mean(LOOISCVs/trueMSEs-1)))
  print(paste0("MPE of LOOVC: ", 100*mean(LOOVCCVs/trueMSEs-1)))
  
  # n = 50:
  # [1] "MSE of LOO: 0.0186379406422529"
  # [1] "MSE of LOOIS: 0.0170571022906853"
  # [1] "MSE of LOOVC: 0.0139447206201104"
  
  # [1] "MAE of LOO: 0.100422072993674"
  # [1] "MAE of LOOIS: 0.0960316894688946"
  # [1] "MAE of LOOVC: 0.0845123591858342"
  
  # [1] "MRE of LOO: 0.671353940243741"
  # [1] "MRE of LOOIS: 0.985694537527453"
  # [1] "MRE of LOOVC: 1.12418200931857"
  
  # [1] "MAPE of LOO: 44.4496728924345"
  # [1] "MAPE of LOOIS: 47.7788979101375"
  # [1] "MAPE of LOOVC: 43.8697950279161"
  
  # [1] "MPE of LOO: -32.8646059756259"
  # [1] "MPE of LOOIS: -1.43054624725475"
  # [1] "MPE of LOOVC: 12.4182009318574"
  
  # n=500:
  # [1] "MSE of LOO: 0.000287803463307309"
  # [1] "MSE of LOOIS: 9.96587931614364e-05"
  # [1] "MSE of LOOVC: 0.000124932796730299"
  
  # [1] "MAE of LOO: 0.0129249857160369"
  # [1] "MAE of LOOIS: 0.00770617706114003"
  # [1] "MAE of LOOVC: 0.00860390665295181"
  
  # [1] "MRE of LOO: 0.594394496022479"
  # [1] "MRE of LOOIS: 1.1115313438249"
  # [1] "MRE of LOOVC: 1.25960186604909"
  
  # [1] "MAPE of LOO: 40.8130559620157"
  # [1] "MAPE of LOOIS: 27.8996173768934"
  # [1] "MAPE of LOOVC: 34.7072374564929"
  
  # [1] "MPE of LOO: -40.5605503977521"
  # [1] "MPE of LOOIS: 11.1531343824896"
  # [1] "MPE of LOOVC: 25.9601866049089"
  
  
  browser()
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
    
    if(n1==50) {
      rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt(1.30)
    } else if(n1==500){
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
       file=paste0("savedOutput/griddedCVtest/n", n, "_iter", iter, unifText, prefText, twoDatText, ".RData"))
  
  list(trueMSE=trueMSE, LOOCVs=LOOCVs, LOOCV=LOOCV, 
       LOOISCV=LOOISCV, LOOVCCV=LOOVCCV, griddedCVs=griddedCVs, 
       gridNs=gridNs, iter=iter, rGRFargsTruth=rGRFargsTruth, 
       rGRFargsSample=rGRFargsSample, 
       n=n, nx=nx, ny=ny, sigmaEpsSq=sigmaEpsSq, 
       allSeeds=allSeeds)
}

# Show why gridded cross validation fails as function of grid resolution. Fix a 
# sampling distribution that is variable on unit square based on GRF, and for a 
# number of different grid resolutions. Error will vary between LOO-CV error 
# and asymptotic extrapolation error, the true predictive error being somewhere 
# in between
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
                         cov.args=list(Covariance="Matern", range=0.2, smoothness=1.0), 
                         delta=3, sigmaEpsSq=0, nx=nx, ny=ny)
  }
  if(is.null(rGRFargsWrong)) {
    rGRFargsWrong = rGRFargsTruth
    
    if(n1 == 50) {
      rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt(4)
    } else if(n1 == 500) {
      rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt(2)
    } else if(n1 == 1000) {
      rGRFargsWrong$sigma = rGRFargsTruth$sigma*sqrt(1.75)
    }
  }
  if(is.null(rGRFargsSample)) {
    rGRFargsSample = list(mu1=0, mu2=0, sigma1=1, sigma2=1, rho=rho, 
                          cov.args=list(Covariance="Matern", range=0.2, smoothness=1.0), 
                          delta=3, sigmaEpsSq1=0, sigmaEpsSq2=0, nx=nx, ny=ny)
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
  
  SigmaSample = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsTruth$cov.args$range, 
                               smoothness=rGRFargsTruth$cov.args$smoothness) * rGRFargsTruth$sigma^2
  SigmaSampleWrong = stationary.cov(xs, xs, Covariance="Matern", aRange=rGRFargsWrong$cov.args$range, 
                               smoothness=rGRFargsWrong$cov.args$smoothness) * rGRFargsWrong$sigma^2
  SigmaSample = SigmaSample + diag(c(rep(sigmaEpsSq1, n1), rep(sigmaEpsSq2, n2)))
  SigmaSampleWrong = SigmaSampleWrong + diag(c(rep(sigmaEpsSq1, n1), rep(sigmaEpsSq2, n2)))
  SigmaGridToSample = stationary.cov(locs, xs, Covariance="Matern", aRange=rGRFargsTruth$cov.args$range, 
                                     smoothness=rGRFargsTruth$cov.args$smoothness) * rGRFargsTruth$sigma^2
  SigmaGridToSampleWrong = stationary.cov(locs, xs, Covariance="Matern", aRange=rGRFargsWrong$cov.args$range, 
                                     smoothness=rGRFargsWrong$cov.args$smoothness) * rGRFargsWrong$sigma^2
  
  # 6. for i in 1:number block resolutions: ----
  
  # determine grid cells associated with observations
  getCellI = function(thisLoc) {
    xI = match(FALSE, thisLoc[1] > highs)
    yI = match(FALSE, thisLoc[2] > highs)
    
    indMat[yI, xI]
  }
  
  griddedCVs = numeric(length(gridNs))
  griddedCVsWrong = numeric(length(gridNs))
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
      
      # predict test data for both true and wrong model
      SigmaAB = matrix(SigmaSample[isTest, !isTest], nrow=sum(isTest))
      SigmaBB = SigmaSample[!isTest, !isTest]
      SigmaAA = SigmaSample[isTest, isTest]
      
      SigmaABWrong = matrix(SigmaSampleWrong[isTest, !isTest], nrow=sum(isTest))
      SigmaBBWrong = SigmaSampleWrong[!isTest, !isTest]
      SigmaAAWrong = SigmaSampleWrong[isTest, isTest]
      
      if((sum(isTest) > 0) && (sum(!isTest) > 0)) {
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
  
  # 10. Calculate LOOIS-CV MSE ----
  cellArea = 1/(nx*ny)
  LOOISrates = (sampleProbs[sampleI1]/sum(sampleProbs)) / cellArea
  # LOOISCV = weighted.mean(LOOCVs, w=cellArea/LOOISrates)
  LOOISCV = sum(LOOCVs * (1/LOOISrates))/n1
  LOOISCVWrong = sum(LOOCVsWrong * (1/LOOISrates))/n1
  
  # 11. Calculate LOOVC-CV MSE ----
  vcellInfo = getVCellAreas(xs1, domainPoly = rbind(c(0, 0), 
                                                    c(0, 1), 
                                                    c(1, 1), 
                                                    c(1, 0), 
                                                    c(0, 0)))
  vcellArea = vcellInfo$area
  rateEsts = 1/vcellArea
  rateEsts = rateEsts/n1 # divide by n to get rate for a single pt draw instead of n pts
  LOOVCCV = sum(LOOCVs * (1/rateEsts))/n1
  LOOVCCVWrong = sum(LOOCVsWrong * (1/rateEsts))/n1
  
  # 12. Calculate true MSE ----
  condDistn = condMeanMVN(SigmaAA=rep(1, nrow(locs)), SigmaAB=SigmaGridToSample, SigmaBB=SigmaSample, 
                          ysB=ys, getFullCov=FALSE, getCondVar=FALSE)
  muAcondB = condDistn$muAcondB
  
  # calculate MSE for the grid cell
  trueMSE = mean((truth - muAcondB)^2) + sigmaEpsSq1
  
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
  save(trueMSE, LOOCVs, LOOCVsWrong, LOOCV, LOOCVWrong, 
       LOOISCV, LOOISCVWrong, LOOVCCV, LOOVCCVWrong, 
       griddedCVs, griddedCVsWrong, gridNs, 
       iter, rGRFargsTruth, rGRFargsSample, rGRFargsWrong, rho, 
       n1, n2, nx, ny, sigmaEpsSq1, sigmaEpsSq2, 
       alpha, beta, allSeeds, 
       file=paste0("savedOutput/griddedCVtest2Datasets/n1", n1, "_n2", n2, "_rho", rho, "_iter", iter, ".RData"))
  
  list(trueMSE=trueMSE, LOOCVs=LOOCVs, LOOCVsWrong=LOOCVsWrong, 
       LOOCV=LOOCV, LOOCVWrong=LOOCVWrong, 
       LOOISCV=LOOISCV, LOOISCVWrong=LOOISCVWrong, 
       LOOVCCV=LOOVCCV, LOOVCCVWrong=LOOVCCVWrong, 
       griddedCVs=griddedCVs, griddedCVsWrong=griddedCVsWrong, 
       gridNs=gridNs, iter=iter, rGRFargsTruth=rGRFargsTruth, 
       rGRFargsSample=rGRFargsSample, rGRFargsWrong=rGRFargsWrong, 
       n1=n1, n2=n2, nx=nx, ny=ny, rho=rho, sigmaEpsSq1=sigmaEpsSq1, 
       sigmaEpsSq2=sigmaEpsSq2, allSeeds=allSeeds)
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
            xI = match(FALSE, thisLoc[1] > highs)
            yI = match(FALSE, thisLoc[2] > highs)
            
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


