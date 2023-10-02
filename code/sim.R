
# simulate GRF on unit square. See ?circulantEmbedding from fields
rGRF = function(nx=100, ny=100, mu=0, sigma=1, 
                cov.args=list(Covariance="Matern", range=0.3, smoothness=1.0), 
                delta=5, sigmaEpsSq=0) {
  require(fields)
  xs = seq(0, 1, l=nx)
  ys = seq(0, 1, l=ny)
  locs = make.surface.grid(list(x = xs, y=ys))
  
  obj<- circulantEmbeddingSetup(list(x = xs, y=ys), cov.args=cov.args, delta=delta)
  out = circulantEmbedding(obj)
  truth = mu + out*sigma
  obs = truth + rnorm(length(out), sd=sqrt(sigmaEpsSq))
  
  list(locs=locs, truth=truth, obs=obs, xs=xs, ys=ys, nx=nx, ny=ny, cov.args=cov.args, mu=mu, sigma=sigma)
}

# fit 2 GRFs on field unit square, one with true parameters, other with false parameters. 
# Goal is twofold:
#  1) to be able to tell the two models apart
#  2) to be able to estimate the expected value of the scoring rule on the unit square
simpleTest = function(rGRFargsTruth=NULL, rGRFargsFalse=NULL, 
                      sampleScenario=c("lgcp", "unif", "horizStrip"), 
                      rGRFargsSample=NULL, ntot=500, ntest=50, 
                      sampleRateEstMethod=c("truth", "vcell", "kernel"), 
                      validationType=c("basic", "thinning", "IS", "PSIS")) {
  
  # Set input args ----
  ntrain = ntot - ntest
  sampleScenario = match.arg(sampleScenario)
  validationType = match.arg(validationType)
  
  if(is.null(rGRFargsFalse)) {
    rGRFargsFalse = list(mu=.1, sigma=1, 
                         cov.args=list(Covariance="Matern", range=0.2, smoothness=0.5), 
                         delta=5, sigmaEpsSq=0)
  }
  if(is.null(rGRFargsTruth)) {
    rGRFargsTruth = list(mu=0, sigma=1, 
                         cov.args=list(Covariance="Matern", range=0.2, smoothness=0.5), 
                         delta=5, sigmaEpsSq=0)
  }
  
  # Generate true field ----
  truth = do.call("rGRF", rGRFargsTruth)
  
  # Do spatial sampling ----
  if(sampleScenario == "lgcp") {
    if(is.null(rGRFargsSample)) {
      # smaller range than the truth for more clumpy observation locations
      rGRFargsSample = list(mu=0, sigma=1, 
                           cov.args=list(Covariance="Matern", range=0.1, smoothness=1.0), 
                           delta=3, sigmaEpsSq=0)
    }
    
    # generate log Gaussian Cox process
    sampleGRF = do.call("rGRF", rGRFargsSample)
    lambdas = exp(sampleGRF$truth)
    # sampleI = sample(1:length(lambdas), ntot, prob=lambdas/sum(lambdas))
    sampleI = which(UPmidzuno(lambdas*(ntot/sum(lambdas))))
    lambdas = lambdas[sampleI]
  } else if(sampleScenario == "unif") {
    lambdas = rep(1, ntot)
    sampleI = sample(1:length(lambdas), ntot)
  } else if(sampleScenario == "horizStrip") {
    # a uniform sampled horizontal strip across the middle third of the unit square
    lambdas = rep(1, ntot)
    sampleI = which((truth$locs[,2] > 1/3) & (truth$locs[,2] < 1/3))
    sampleI = sample(sampleI, ntot)
  }
  
  # get samples of the true field
  sampleLocs = truth$locs[sampleI,]
  sampleYs = truth$truth[sampleI]
  
  # Estimate sample probability ----
  # (or at least something roughly proportional to it)
  if(sampleRateEstMethod == "truth") {
    ps = lambdas
  } else if(sampleRateEstMethod == "vcell") {
    stop("Voronoi cell sample rate estimates not yet supported")
  } else if(sampleRateEstMethod == "kernel") {
    stop("kernel sample rate estimates not yet supported")
  }
  # ps = ps/sum(ps)
  
  # Partition sample into train/test ----
  # also calculate validation weights
  if(validationType == "basic") {
    # SRS sampling of test indices from sample
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
  } else if(validationType == "ISvcell") {
    stop(paste0(validationType, " validation type not yet supported"))
  } else if(validationType == "PSISvcell") {
    stop(paste0(validationType, " validation type not yet supported"))
  } else if(validationType == "ISkernel") {
    stop(paste0(validationType, " validation type not yet supported"))
  } else if(validationType == "PSISkernel") {
    stop(paste0(validationType, " validation type not yet supported"))
  }
  trainI = setdiff(1:ntot, testI)
  
  # Get train/test xs/ys ----
  trainLocs = sampleLocs[trainI,]
  testLocs = sampleLocs[testI,]
  trainYs = sampleYs[trainI,]
  testYs = sampleYs[testI,]
  
  # Get predictions ----
  
}


