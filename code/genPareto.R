dGenPareto = function(y, u, k, sigma=theta*k, theta=sigma/k, doLog=FALSE) {
  if(k == 0) {
    out = -log(sigma) + (y - u)/sigma
  } else {
    out = -log(sigma) - (1 + 1/k) * log(1 + k*((y - u)/sigma))
  }
  
  if(doLog) {
    out
  } else {
    exp(out)
  }
}

pGenPareto = function(q, u, k, sigma=theta*k, theta=sigma/k) {
  if(k == 0) {
    1 - exp(-(q - u)/sigma)
  } else {
    1 - (1 + k*(q - u) / sigma)^(-1/k)
  }
}

qGenPareto = function(p, u, k, sigma=theta*k, theta=sigma/k) {
  if(k == 0) {
    u - sigma*log((1-p))
  } else {
    u + sigma*((1-p)^(-k) - 1)/k
  }
}

fitGenParetoThetaK = function(rs, u) {
  
  optFun = function(par) {
    theta = par[1]
    k = par[2]
    -sum(dGenPareto(rs, u, k, theta=theta, doLog=TRUE))
  }
  
  opt = optim(c(1, .5), fn=optFun)
  
  opt$par
}

# construct Pareto smoothed importance weights
# rs: weights
PSISsimple = function(rs) {
  S = length(rs) # sample size
  M = min(c(.2*S, 3*sqrt(S))) # smooth M largest weights
  
  # sort weights, get threshold just lower than highest M weights
  ordI = order(rs)
  rsSort = rs[ordI]
  uTemp = rsSort[S - M + 1]
  u = max(rs[rs < uTemp]) # the threshold
  
  # fit gen Pareto distn to M highest weights
  highRs = rsSort[rsSort > u]
  thetaK = fitGenParetoThetaK(highRs, u=u)
  theta = thetaK[1]
  k = thetaK[2]
  sigma = theta*k
  
  print(paste0("estimated k: ", k))
  
  # M order statistics of the generalized Pareto
  Mqs = (1:M - 0.5) * (1/M)
  Mords = qGenPareto(Mqs, u, k, sigma=sigma)
  
  # replace M highest weights with the M order statistics
  rsOut = rs
  rsOut[ordI[(S - M + 1):S]] = Mords
  
  # return results
  list(rs=rsOut, u=u, k=k, sigma=sigma, theta=theta)
}


