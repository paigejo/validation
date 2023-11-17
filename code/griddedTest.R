
# illustrating the problem of gridded/blocked CV



estPredMSErange = function(sigma=1, cov.args=list(Covariance="Matern", range=1, smoothness=1.0), 
                           dists=exp(seq(log(.01), log(cov.args$range*6), l=100)), nPts=500, doPlot=TRUE) {
  
  # make points in unit radius ring around origin
  thetaRes = 2 * pi / nPts
  thetas = seq(0, 2*pi-thetaRes, l=nPts)
  unitXs = cos(thetas)
  unitYs = sin(thetas)
  unitPts = cbind(unitXs, unitYs)
  unitDistsAB = rdist(cbind(0, 0), unitPts)
  unitDistsBB = rdist(unitPts)
  
  # for each distance level, calculate conditional variance
  varMins = numeric(length(dists))
  varMaxes = varMins
  for(i in 1:length(dists)) {
    print(paste0("i: ", i, "/", length(dists)))
    r = dists[i]
    thisDistsAB = r * unitDistsAB
    thisDistsBB = r * unitDistsBB
    
    # calculate covariances
    SigmaAA = sigma^2
    SigmaAB = sigma^2 * do.call("stationary.cov", c(list(cbind(0, 0), unitPts, distMat=thisDistsAB), cov.args))
    SigmaBB = sigma^2 * do.call("stationary.cov", c(list(unitPts, distMat=thisDistsBB), cov.args))
    
    # SigmaA|B = SigmaAA - SigmaAB SigmaBB^-1 SigmaBA
    #          = SigmaAA - SigmaAB (U' U)^-1 SigmaBA
    #          = SigmaAA - SigmaAB U^-1 (U')^-1 SigmaBA
    #          = SigmaAA - R' R
    U = chol(SigmaBB)
    R = forwardsolve(t(U), t(SigmaAB))
    
    # (R' R)_ii = (R')_i: R_:i
    varMins[i] = SigmaAA - myDiag(apply(R, 2, function(x) {sum(x^2)}))
    
    # do the same when conditioning on just a single point
    varMaxes[i] = sigma^2 - SigmaAB[1,2] * (1/sigma^2) * SigmaAB[1,2]
  }
  
  if(doPlot) {
    # calculate expected range
    # corDists = seq(0, 5*cov.args$range, l=3000)
    # cors = do.call("stationary.cov", c(list(cbind(0, 0), cbind(corDists, 0)), cov.args))
    # effInd = match(TRUE, cors < .1)
    # effRange = corDists[effInd]
    # plot(corDists, 1-cors, type="l", xlim=range(dists), ylim=c(0, sigma^2), 
    #      main="Variogram", xlab="Distance", ylab="", log="x")
    # lines(dists, varMins, lty=2, col="blue")
    # lines(dists, varMaxes, lty=2, col="blue")
    # abline(v=effRange, col="green")
    # abline(h=.9, col="green", lty=2)
    # 
    # plot(corDists, sqrt(1-cors), type="l", xlim=range(dists), ylim=c(0, sigma), 
    #      main="Sqrt variogram", xlab="Distance", ylab="", log="x")
    # lines(dists, sqrt(varMins), lty=2, col="blue")
    # lines(dists, sqrt(varMaxes), lty=2, col="blue")
    # abline(v=effRange, col="green")
    # abline(h=.9, col="green", lty=2)
    # 
    plot(dists, varMaxes, type="l", lty=2, col="blue", main="MSE range vs distance", 
         xlab="Distance to nearest observation", ylab="MSE", log="x", ylim=c(0,sigma^2))
    lines(dists, varMins, lty=2, col="blue")
    # abline(v=effRange, col="green")
    # abline(h=.9, col="green", lty=2)
    
    rmseRange = c(0, sqrt(max(varMaxes)))
    plot(dists, sqrt(varMaxes), type="l", lty=2, col="blue", main="RMSE range vs distance", 
         xlab="Distance to nearest observation", ylab="RMSE", log="x", ylim=c(0, sigma))
    lines(dists, sqrt(varMins), lty=2, col="blue")
    # abline(v=effRange, col="green")
  }
  
  list(dists=dists, minMSE=varMins, maxMSE=varMaxes)
}