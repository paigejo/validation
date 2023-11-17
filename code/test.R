testRGRF = function(nx=100, ny=100, mu=0, sigma=1, testN=50, 
                    cov.args=list(Covariance="Matern", range=0.2, smoothness=1.0), 
                    delta=5, sigmaEpsSq=0, alpha=0, beta=1, seed=123) {
  
  set.seed(123)
  
  out = rGRF(nx=nx, ny=ny, mu=mu, sigma=sigma, 
             cov.args=cov.args, 
             delta=delta, sigmaEpsSq=sigmaEpsSq)
  
  rates = exp(alpha + beta * out$truth)
  sampleProbs = rates/sum(rates)
  gridResX = 1/nx
  gridResY = 1/ny
  sampleI = sample(1:nrow(out$locs), testN, prob=sampleProbs, replace=TRUE)
  xs = out$locs[sampleI,] + cbind(runif(testN, max=gridResX)-gridResX/2, 
                              runif(testN, max=gridResY)-gridResY/2)
  
  pdf(paste0("figures/test/testGRFpref_range", cov.args$range, "_smooth", cov.args$smoothness, 
      "_a", alpha, "_b", beta, ".pdf"), width=8, height=5)
  par(mfrow=c(1,2))
  quilt.plot(out$locs[,1], out$locs[,2], out$truth, nx=nx, ny=ny, 
             main=paste0("GRF, range=", cov.args$range, 
                         ", smoothness=", cov.args$smoothness))
  quilt.plot(out$locs[,1], out$locs[,2], sampleProbs, nx=nx, ny=ny, 
             main=paste0("Sample probs, a=", alpha, ", b=", beta))
  points(xs)
  dev.off()
}

testR2GRFs = function(nx=100, ny=100, mu1=0, mu2=1, sigma1=1, sigma2=2, testN=50, 
                    cov.args=list(Covariance="Matern", range=0.2, smoothness=1.0), 
                    delta=5, sigmaEpsSq1=0, sigmaEpsSq2=.1^2, alpha=0, beta=1, rho=-.8, seed=NULL) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  out = r2GRFs(nx=nx, ny=ny, mu1=mu1, mu2=mu2, sigma1=sigma1, sigma2=sigma2, 
             cov.args=cov.args, rho=rho, 
             delta=delta, sigmaEpsSq1=sigmaEpsSq1, sigmaEpsSq2=sigmaEpsSq2)
  
  rates1 = exp(alpha + beta * out$Y1$truth)
  rates2 = exp(alpha + beta * out$Y2$truth)
  sampleProbs1 = rates1/sum(rates1)
  sampleProbs2 = rates2/sum(rates2)
  gridResX = 1/nx
  gridResY = 1/ny
  sampleI1 = sample(1:nrow(out$Y1$locs), testN, prob=sampleProbs1, replace=TRUE)
  sampleI2 = sample(1:nrow(out$Y1$locs), testN, prob=sampleProbs2, replace=TRUE)
  xs1 = out$Y1$locs[sampleI1,] + cbind(runif(testN, max=gridResX)-gridResX/2, 
                                  runif(testN, max=gridResY)-gridResY/2)
  xs2 = out$Y1$locs[sampleI2,] + cbind(runif(testN, max=gridResX)-gridResX/2, 
                                   runif(testN, max=gridResY)-gridResY/2)
  
  pdf(paste0("figures/test/test2GRFspref_range", cov.args$range, "_smooth", cov.args$smoothness, 
             "_a", alpha, "_b", beta, "_rho", rho, ".pdf"), width=8, height=5)
  par(mfrow=c(2,2))
  quilt.plot(out$Y1$locs[,1], out$Y1$locs[,2], out$Y1$truth, nx=nx, ny=ny, 
             main=paste0("GRF, range=", cov.args$range, 
                         ", smoothness=", cov.args$smoothness))
  quilt.plot(out$Y1$locs[,1], out$Y1$locs[,2], sampleProbs1, nx=nx, ny=ny, 
             main=paste0("Sample probs, a=", alpha, ", b=", beta))
  points(xs1)
  
  quilt.plot(out$Y1$locs[,1], out$Y1$locs[,2], out$Y2$truth, nx=nx, ny=ny, 
             main=paste0("GRF, range=", cov.args$range, 
                         ", smoothness=", cov.args$smoothness))
  quilt.plot(out$Y1$locs[,1], out$Y1$locs[,2], sampleProbs2, nx=nx, ny=ny, 
             main=paste0("Sample probs, a=", alpha, ", b=", beta))
  points(xs2)
  dev.off()
}