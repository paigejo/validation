
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

# Use IWs:
#   1/rateEsts from full dataset
# Use control variates: 
#   rateEsts (with mean 1)
#   rateEsts2 (with mean 1)
#   -------1/rateEstsVC from train data
#   -------(and possibly rateEstsVC2/rateEstsVC from full dataset)
# given n-vector of scores and n-vector of respective importance weights, 
# uses weights as a control variate to reduce variance of IWCV estimator as in:
# https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8918731
# (title: ROBUST IMPORTANCE-WEIGHTED CROSS-VALIDATION UNDER SAMPLE SELECTION BIAS)
# and:
# https://pdf.sciencedirectassets.com/271709/1-s2.0-S0305054800X01220/1-s2.0-0305054887900244/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEF0aCXVzLWVhc3QtMSJHMEUCIQDYIUilf2rZt3kp5BtZJSc0pbR9LHBxj%2Bi4dP6lgOMwNQIgI1Zj8iiXCIPmm3df5aAEcEEDeiZUwwyU9Fi1uoOpqgkquwUIpv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDIMockogYU8iSbcZRiqPBalptBhbgcZk%2FOHN3L673fm1BPguBvCWT%2BiVMuQsJ3%2FJFed2lJZMajVHvy9K6fzDfdqV5wh1Sw2F1gD%2BQAlcAGMWsyL6KpOeHZFFvY%2FzWZnNncFAoh6YAfvNeS2CJUOEDbugd5UfpkmVrjjg2ZLs%2FrcJPAUJfQ5Qv5pJ1QUIc%2FOVTkMuvm%2F1lVYmLiKWWySniWUZmWDz5ncjo2NnrBQxVXn5mpsKLHHntKoOX0LXDlZtV%2BW8Co7aWDsZ6GuKaX9UXEYlGOBQ91Ub0e%2FB7wt1SXAHjUjrE6Fu3rUQ84Y9XKPTpf4gSuqGv757B%2BhL3CTGjksCFZfVfYr62q%2BFaiu5liUXbrgsajB%2B%2Fs%2F1mUjCiVEw1slA%2B4hxuCEJvqVzvAglnBtZY7YO5IYmGQhbl8Om%2FTLWmsoYBDLeQFIkBlwRCIMIT8BK149BWO9MIhwzMrbdDuBA2stJ7lw55TdlgDOaP3LxDmXbZdqROBdMGOumnS3FS2kYS0S9wBeiNeqnGue%2BLl74Vi%2FSczzhVR2t8mgZ2A%2Bwtwtz6sfkJx%2ByB7jDBZq2INHHweVVxJJlgA60wElWQTDSdyNLYkwGgWT3vFTfnWTtgzP2Amzb9ryND1Jbta%2FvgPthOYXJIx%2B5pz96Iaj%2FAN3%2B9A6IjMzLasMz5shw%2FYO96mYsswa%2Bts81XDlFYSm9F3Tu%2FjhNw2XB8eVmz0OM4BHYW1ulb%2FlJ3Err9M09hADniuH4t%2B4vM6oL6dRMBBjxwn%2FpJpah6%2BhV41VA3AUAaH5V1djx2%2FhwJnpG0VLRtrPeQk544om2pKQlvc5hoDglfsdprBRyAUQijyloqbVeEdCOwqAA9fKP%2BZ2sO744QS4lJECZcwBX41IUmJ6k%2FaswprjdqgY6sQEAQ6hqEkSyROTr0v6G1W2zp2elXtsYPQz304FZ6ZhJDmsVO7oO7LWBGbXP7sXfp1JwC9KzSyOTlSKOY%2BFizTs3W1DN4rF0mIgOA01vdfGVssP6E%2B8Aj7F5KK7%2BSo%2FQwVuQ2vg3MfExbofJPuD356MMTeTsDmNReR4idICgON%2FlBTQHDOT8LtsVlNWWlO%2BDIZ7%2FrKKnjA8eYQVPW2GaXivaBXTpc2eMB56nV0nTr6n5CwQ%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20231117T130334Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYSWMLUGWC%2F20231117%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=2c444bcebe82a96504812f09276f43239537a7928ef30f321a7abb453b602456&hash=426e056f6481a649660b90486dc0fe46510d318a3d1f4e91f3d5a8ea5f265173&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=0305054887900244&tid=spdf-fd01066f-9212-46b8-b3b5-cf666bcfc21b&sid=2bed538a4721394ae01b7b10291f68b4853bgxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=1409595c58575c00030700&rr=82782caffbc3712d&cc=no
# (title: ON CONTROL VARIATE ESTIMATORS, by Nelson)
# and, for multivariate control variates:
# https://pdf.sciencedirectassets.com/271697/1-s2.0-S0167637700X01356/1-s2.0-0167637786900982/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEKP%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIFUBxwUyoIChChNo%2BLywmH52yGjc%2FBLuuZvuOyseYg%2BMAiEA78zBPAmEufDLYXNqp6ZoXalbrt7xPugg41jLVbCM%2FRMquwUI7P%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDPMwSIwuAAjJnPGNACqPBfO6yl%2FQ%2BIyUxZPsCO2ryQTzSRSqdH71EJmzpZUtitsgrUPOxgLKNe41cIRpm9AYUiP6JE6LhEmZpselRd%2FvNAcUh%2BeK5g5xQT8sl%2Ft6%2Bp%2FHYd0bTtm9QlYEcP3QDaEUkjXDVIx159cVvg6ezi5sHW7bNGJmf7%2BTuwqgGkODQ5Squsxwguzg%2ByTSr81qT03rGek%2FAZM0nlEaUeXb0vPEhZqbcFItKUkFJZ5wiExXULmT%2BMDEXOY7SNFU4m9r%2BWvqC%2B%2FnRFJlknn1Z9yB6q4hTDUzVPw8vm1R0aaLUgUWTHHO0ofk%2FDi94htY26HujdrWuWG43QU4EpKodVwr1VjmeFN1Yw7WVAX3GcbY8TCu%2FqhKxZ0scHd%2F9BwMdhxdhgLO3CSvM4EbFwG%2FLZSS6mYianZaT%2F8TjqQcO%2Fmx47iKmENqatFmsp2r7ShGXwlvKk1DiYW4ohX4en7k%2BUwRfzJnThrPF2sfScA6Cwkc1SCwdnT52BFNg22YwzAffHUfG%2FYs9zP8Uq0ll1DcHH03ZsP9ePEd8F%2B5hG%2BeyF57cQUBCI%2BsUWJnB0mAt%2FVHon3m6SPKMi%2FKOLtNFvrNMin8IPEgIqChLlYTltF7lxgivFK8UvgoDb52oSD5eTgZFtgfno37zJSFEaDD%2B62yW70WESeWbsonxxUwPacsIhp%2BXzvAvYaIwsXlSZsju7Gy7WnLy2KXhYBImMUEocd4ZosvPo93%2Bfms7%2BZmXmEr%2F8IK7H8GbZuLZW4sZuJjFXnOVFP%2Bo%2BqKGqhtso%2BUgOXQV%2F1S%2FXZSv1KSKGzIelDIl%2Bjj5%2FgvCPpJuhavUqKKjH54FE%2BMo3L6RwDEH4D%2F1rY3xfVakw1PA6YIEZ9WKa%2BFnywsSf6nHlcwovTsqgY6sQFFr6y4g7yPlCq3VtsUZVlidJV2vz%2B7PrHoko%2Ftzhb8uMyXrC%2BL2FukrohULZQSzL19R0xiqcgGp4r%2BgR1KtwfUZgCougvs1u%2FCOkE17xioVJhq%2FYfKaGTot%2FIXzhz%2Fal5aEIvAanii8MGnFPk2KsOS5oag0XIyXjdctj%2BHS7yeoExfDOQnGSAqLqGosMyOCgWCPSYmN06gvbPUDhrSG5Bv4KtYfmvyx1QlsANW80dLdWg%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20231120T114509Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYV4GG6XJE%2F20231120%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=d8c8004072e5a3954f437a6e9d95fb25dd18653149b47a0f110fb1371a905d04&hash=0350a15bda55c1d3fd34e904d776ecea05d941d219460f179fe792477a2eaed0&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=0167637786900982&tid=spdf-9dc4ae03-ce71-4055-b5c0-12b7e0a178d5&sid=56a9d21780efc14153682c316db1e857d7ffgxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=1409595c565f5952070750&rr=829071ef6832712a&cc=no
# (title: THE EFFICIENCY OF CONTROL VARIATES IN MULTIRESPONSE SIMULATION)
getRobustIWCV = function(scores, IWs=1, controlVarMat=IWs, controlVarMeans=1, wtdReg=FALSE, 
                         tol = 1e-08, printVarRatio=FALSE, printRoot="", wtdMeanCentered=TRUE, 
                         returnInfo=FALSE) {
  k = length(scores)
  
  if(!wtdReg && !wtdMeanCentered) {
    Ys = scores * IWs
    controlVarMat = sweep(cbind(controlVarMat), 1, IWs, "*")
    regWs = rep(1, k)
  } else if(!wtdMeanCentered) {
    Ys = scores
    regWs = IWs
  } else if(wtdMeanCentered) {
    Ys=scores
    regWs = rep(1, k)
  }
  
  if(!is.matrix(controlVarMat)) {
    q = 1 # number of control variates
    if(length(scores) -2 - 1 <= 0) {
      warning("fewer observations than control vars. Using 0 control vars...")
      return(mean(IWs*scores))
    }
    controlVarMat = matrix(controlVarMat, ncol=1)
    # NOTE: it's not exactly a classic SLR problem, since we know E[IWs] = 1
    # mod = lm(weightedScores ~ IWs)
    # beta = mod$coefficients[2]
    # beta = sum((scores*IWs - mean(scores*IWs))*(IWs-1))/sum((IWs-1)^2)
    # Ybar = weighted.mean(Ys)
    # beta = sum((Ys - mean(Ys))*(controlVarMat-controlVarMeans))/sum((controlVarMat-controlVarMeans)^2)
    
    # out = c(scores * IWs - beta * (controlVarMat - controlVarMeans))
  }
  
  q = ncol(controlVarMat)
  if(length(scores) -2 - q <= 0) {
    warning("fewer observations than control vars. Using 0 control vars...")
    return(mean(IWs*scores))
  }
  if(length(controlVarMeans) == 1) {
    controlVarMeans = rep(controlVarMeans, q)
  }
  
  # using notation from "THE EFFICIENCY OF CONTROL VARIATES IN MULTIRESPONSE SIMULATION"
  muC = controlVarMeans
  C = controlVarMat
  
  # empirical covariance matrix of C (here, muC is known)
  SC = matrix(0, nrow=ncol(C), ncol=ncol(C))
  if(!wtdMeanCentered) {
    for(j in 1:k) {
      SC = SC + regWs[j] * outer(C[j,] - muC, C[j,] - muC)
    }
    # meanC = mean(C)
    # for(j in 1:k) {
    #   SC = SC + regWs[j] * outer(C[j,] - meanC, C[j,] - meanC)
    # }
  } else {
    meanC = sum(C*IWs)/k
    for(j in 1:k) {
      SC = SC + regWs[j] * outer(C[j,] - meanC, C[j,] - meanC)
    }
  }
  
  SC = 1/(k-1) * SC
  
  # check if SC singular (e.g. if control variates are constant, as may 
  # occur in blocked CV)
  if(abs(det(SC)) < tol) {
    warning("controlVarMat not full rank. Using 0 control vars...")
    return(mean(IWs*scores))
  }
  
  # empirical covariance matrix of Y and C
  SYC = matrix(0, nrow=1, ncol=ncol(C))
  if(!wtdMeanCentered) {
    Ybar = sum(Ys*regWs)/k
    for(j in 1:length(Ys)) {
      SYC = SYC + regWs[j] * outer(Ys[j] - Ybar, C[j,] - muC)
    }
    # for(j in 1:length(Ys)) {
    #   SYC = SYC + regWs[j] * outer(Ys[j] - Ybar, C[j,] - meanC)
    # }
  } else {
    Ybar = sum(Ys*IWs)/k
    for(j in 1:length(Ys)) {
      SYC = SYC + regWs[j] * outer(Ys[j] - Ybar, C[j,] - meanC)
    }
  }
  
  SYC = matrix((1/(k-1)) * SYC, ncol=ncol(C))
  
  # estimate of betaVec
  B = SYC %*% solve(SC)
  # if(k == 50) {
  #   browser()
  # }
  # calculate controlled IWCV
  out = c(Ys - sweep(controlVarMat, 2, controlVarMeans, FUN="-") %*% c(B))
  
  if(printVarRatio) {
    lambda = (k-2)/(k - 2 - q)
    varControlEst = var(out)/k
    varUncontrolEst = var(Ys)/k
    varRatio = lambda * (varControlEst/varUncontrolEst)
    print(paste0("controlled estimate resulted in variance ratio of ~", round(varRatio, 4), " for ", printRoot))
  }
  
  if(!returnInfo) {
    sum(out * regWs) / k 
  } else {
    # wis = IWs
    # gisFirstPt = matrix(1 + t(c(controlVarMeans - tXhat)) %*% BhatXpt, nrow=1)
    # gis = c(gisFirstPt %*% t(controlVarMat))
    # eis = scores - (tYhat + sweep(controlVarMat, 2, controlVarMeans, "-") %*% Bhat)
    # tYhatGREG = tYhat + sum(c(controlVarMeans - tXhat) * c(Bhat))
    # list(tYhatGREG=tYhatGREG, tXhat=tXhat, tYhat=tYhat, wis=wis, gis=gis, eis=eis)
    
    eis = out
    YhatRobust = sum(out * regWs) / k 
    varHatYhat = var(eis * regWs/k)
    list(YhatRobust=YhatRobust, varHatYhat=varHatYhat)
  }
}

# Could be improved by accounting for known pop means, and using sandwich 
# estimation in MSE estimation.
# ws: importance/survey weights with expected value of 1 over the data distribution. 
#     It is recommended to include ws in Zg
getWeightedControlVarEst = function(Ys, Zf, muf, Zg=NULL, mug=NULL, ws=rep(1, length(Ys)), 
                                    shrinkWeights=FALSE, includeIntercept=TRUE, 
                                    pSeq=c(0, expit(seq(logit(.001), logit(.999), l=1000)), 1), 
                                    estVar=TRUE, printPhat=TRUE) {
  
  n = length(Ys)
  if(!is.matrix(Zf)) {
    Zf = matrix(Zf, nrow=length(Ys))
  }
  if(!is.null(Zg) && !is.matrix(Zg)) {
    Zg = matrix(Zg, nrow=length(Ys))
  }
  
  # calculate original weighted estimate
  thetaHatOrig = sum(Ys*ws)/n
  
  # calculate regularization parameter pHat
  if(shrinkWeights) {
    estMSE = function(thisP) {
      out = getWeightedControlVarEst(Ys, Zf, muf, Zg, mug, ws^thisP, 
                                     shrinkWeights=FALSE, includeIntercept, estVar=TRUE)
      thisThetaHat = out[1]
      thisVarHat = out[2]
      (thisThetaHat - thetaHatOrig)^2 + thisVarHat
    }
    
    mses = sapply(pSeq, estMSE)
    minI = which.min(mses)
    pHat = pSeq[minI]
    
    out = getWeightedControlVarEst(Ys, Zf, muf, Zg, mug, ws^pHat, 
                                   shrinkWeights=FALSE, includeIntercept, estVar=TRUE)
    out[2] = mses[minI]
    
    if(printPhat) {
      print(paste0("pHatC: ", pHat))
    }
    
    return(out)
  }
  
  # concatenate effects, since they will all be estimated together anyway
  Z = cbind(Zf, Zg)
  mu = c(muf, mug)
  
  # fit weighted least squares model
  # get estimated coefficients aside from intercept
  if(includeIntercept) {
    mod = lm(Ys~Z, weights=ws)
    betaHats = coef(mod)[-1]
    # Ztemp = cbind(1, Z)
    # mutemp = c(1, mu)
    # 
    # if(!is.null(mug)) {
    #   gInds = (1:length(mug)) + length(muf) + 1
    # } else {
    #   gInds = NULL
    # }
  } else {
    mod = lm(Ys~Z-1, weights=ws)
    betaHats = coef(mod)
    # Ztemp = Z
    # mutemp = mu
    # if(!is.null(mug)) {
    #   gInds = (1:length(mug)) + length(muf)
    # } else {
    #   gInds = NULL
    # }
  }
  # betaHatsg = betaHats[-(1:length(muf))]
  # betaHatsf = betaHats[1:length(muf)]
  
  # empirical covariance matrix of Z (here, muf might be known, but mug must be estimated)
  # SZ = matrix(0, nrow=ncol(Ztemp), ncol=ncol(Ztemp))
  # if(!wtdMeanCentered) {
  #   if(!is.null(mug)) {
  #     mugHatPop = colSums(sweep(Zg, 1, ws, "*")) * (1/n)
  #     mutemp[gInds] = mugHatPop
  #   }
  #   
  #   for(j in 1:k) {
  #     SZ = SZ + ws[j] * outer(Ztemp[j,] - mutemp, Ztemp[j,] - mutemp)
  #   }
  #   
  #   SZ = 1/(n-length(mug)) * SZ
  # } else {
  #   mutempHat = sum(Ztemp*ws)/n
  #   for(j in 1:k) {
  #     SZ = SZ + ws[j] * outer(Ztemp[j,] - mutempHat, Ztemp[j,] - mutempHat)
  #   }
  #   
  #   if(includeIntercept) {
  #     SZ = 1/(n-length(mutempHat)+1) * SZ # add one since the intercept mean is known
  #   } else {
  #     SZ = 1/(n-length(mutempHat)) * SZ
  #   }
  #   
  # }
  # 
  # # empirical covariance matrix of Y and Ztemp
  # SYZ = matrix(0, nrow=1, ncol=ncol(Ztemp))
  # if(!wtdMeanCentered) {
  #   Ybar = sum(Ys*regWs)/k
  #   for(j in 1:length(Ys)) {
  #     SYZ = SYZ + regWs[j] * outer(Ys[j] - Ybar, Ztemp[j,] - mutemp)
  #   }
  #   # for(j in 1:length(Ys)) {
  #   #   SYZ = SYZ + regWs[j] * outer(Ys[j] - Ybar, Ztemp[j,] - meanC)
  #   # }
  # } else {
  #   Ybar = sum(Ys*IWs)/k
  #   for(j in 1:length(Ys)) {
  #     SYZ = SYZ + regWs[j] * outer(Ys[j] - Ybar, Ztemp[j,] - mutempHat)
  #   }
  # }
  # SYZ = matrix((1/(n-length(c(muf, mug)))) * SYZ, ncol=ncol(Ztemp))
  # 
  # # # estimate of betaVec
  # B = SYZ %*% solve(SZ)
  
  # calculate mugHat and mufHat, expected values of Zg and Zf under population distribution
  if(!is.null(Zf)) {
    mufHat = apply(Zf, 2, weighted.mean, w=ws) * (sum(ws)/n)
  } else {
    mufHat = NULL
  }
  if(!is.null(Zg)) {
    mugHat = colMeans(Zg)
  } else {
    mugHat = NULL
  }
  
  muHat = c(mufHat, mugHat)
  
  # calculate final estimate
  thetaHat = thetaHatOrig - (muHat - mu) %*% betaHats
  
  # calculate the MSE if need be based on V_1 of p. 459 of Lohr
  if(estVar) {
    eis = resid(mod)
    if(!is.null(mug)) {
      yHatgs = Zg %*% betaHats[-(1:length(mufHat))]
      varEst = var(eis*ws + (ws - (1/n))*yHatgs)/n
    } else {
      varEst = var(eis*ws)/n
    }
    
    c(thetaHat=thetaHat, thetaHatVar=varEst)
  } else {
    thetaHat
  }
}

# p. 458 of Lohr Sampling
getGREGCV = function(scores, IWs, controlVarMat=rep(1, length(scores)), controlVarMeans=1, 
                     normalizeWeights=FALSE, returnInfo=FALSE, shrinkWeights=FALSE, 
                     logitPSeq=seq(logit(.001), logit(.999), l=1000), type=1) {
  k = length(scores)
  
  # if(shrinkWeights) {
  #   out = getPGREG(LOOCVs=scores, IWs=IWs, controlVarMat=controlVarMat, controlVarMeans=controlVarMeans, 
  #                  logitPSeq=logitPSeq, normalizeWeights=normalizeWeights, type=type)
  #   IWs = IWs^out$p
  # }
  
  if(normalizeWeights) {
    IWs = IWs / sum(IWs)
  }
  
  if(!is.matrix(controlVarMat)) {
    controlVarMat = matrix(controlVarMat, ncol=1)
  }
  q = ncol(controlVarMat)
  
  # calculate Bhat, one part at a time
  BhatXpt = matrix(0, ncol=q, nrow=q)
  for(i in 1:k) {
    Xi = matrix(controlVarMat[i,], ncol=1)
    wi = IWs[i]
    
    BhatXpt = BhatXpt + wi * Xi %*% t(Xi)
  }
  BhatXpt = solve(BhatXpt)
  
  BhatXYpt = matrix(0, ncol=1, nrow=q)
  for(i in 1:k) {
    Xi = matrix(controlVarMat[i,], ncol=1)
    wi = IWs[i]
    Yi = scores[i]
    
    BhatXYpt = BhatXYpt + wi * Xi * Yi
  }
  
  Bhat = BhatXpt %*% BhatXYpt
  
  # calculate Horvitz-Thompson estimates of X and Y
  if(!normalizeWeights) {
    tYhat = (1/k) * sum(IWs * scores)
    tXhat = (1/k) * apply(controlVarMat, 2, function(x) {sum(x * IWs)})
  } else {
    tYhat = sum(IWs * scores)
    tXhat = apply(controlVarMat, 2, function(x) {sum(x * IWs)})
  }
  
  if(shrinkWeights) {
    wis = IWs
    gisFirstPt = matrix(1 + t(c(controlVarMeans - tXhat)) %*% BhatXpt, nrow=1)
    gis = c(gisFirstPt %*% t(controlVarMat))
    eis = scores - (tYhat + sweep(controlVarMat, 2, controlVarMeans, "-") %*% Bhat)
    tYhatGREG = tYhat + sum(c(controlVarMeans - tXhat) * c(Bhat))
    unbiasedEstimate = tYhat
    
    # browser()
    out = getP(LOOCVs=scores, IWs=wis*gis, normalizeIWs=FALSE, unbiasedEstimate=unbiasedEstimate, 
               logitPSeq=logitPSeq, resids=eis)
  }
  brow
  # calculate GREG estimator for Y based on X
  if(!returnInfo) {
    tYhat + sum(c(controlVarMeans - tXhat) * c(Bhat))
  } else {
    wis = IWs
    gisFirstPt = matrix(1 + t(c(controlVarMeans - tXhat)) %*% BhatXpt, nrow=1)
    gis = c(gisFirstPt %*% t(controlVarMat))
    eis = scores - (tYhat + sweep(controlVarMat, 2, controlVarMeans, "-") %*% Bhat)
    YhatGREG = tYhat + sum(c(controlVarMeans - tXhat) * c(Bhat))
    if(!normalizeWeights) {
      varHatYhat = var(eis*wis/k)
    } else {
      varHatYhat = var(eis*wis)
    }
    list(YhatGREG=YhatGREG, varHatYhat=varHatYhat)
  }
}

# DEPRECATED. Instead, use getP using the adjusted GREG weights (wis * gis) and 
# with normalizeIWs=FALSE. You could also set unbiasedEstimate=YhatW
getPGREG_old = function(LOOCVs, IWs, controlVarMat=IWs, controlVarMeans=1, 
                    logitPSeq=seq(logit(.001), logit(.999), l=1000), 
                    normalizeWeights=TRUE, type=1, test=TRUE) {
  k = length(LOOCVs)
  allPs = c(0, expit(logitPSeq), 1)
  VarY = var(LOOCVs)
  
  YhatW = (1/sum(IWs)) * sum(LOOCVs * IWs)
  
  # commented out because it assumes the weights are fixed
  # estMSE = function(p=NULL, logitP=NULL) {
  #   if(is.null(p)) {
  #     p = expit(logitP)
  #   }
  #   YhatP = (1/sum(IWs^p)) * sum(LOOCVs * IWs^p)
  #   Bp2 = (YhatP - YhatW)^2
  #   
  #   Bp2 + VarY * (sum(IWs^(2*p))/sum(IWs^p)^2)
  # }
  
  # accounts for randomness in the weights and correlation with LOOCVs
  estMSE = function(p=NULL, logitP=NULL) {
    if(is.null(p)) {
      p = expit(logitP)
    }
    ws = IWs^p/sum(IWs^p)
    
    gregInfo = getGREGCV(LOOCVs, IWs, controlVarMat, controlVarMeans, normalizeWeights=TRUE, returnInfo=TRUE)
    gis = gregInfo$gis
    eis = gregInfo$eis
    YhatP = gregInfo$tYhatGREG
    
    Bp2 = (YhatP - YhatW)^2
    
    if(type == 1) {
      VarYhatP = var(eis * ws) * k
    } else {
      VarYhatP = var(eis * gis * ws) * k
    }
    
    Bp2 + VarYhatP
  }
  
  if(test) {
    type=1
    allMSEs1 = sapply(allPs, estMSE)
    bestI1 = which.min(allMSEs1)
    type=2
    allMSEs2 = sapply(allPs, estMSE)
    bestI2 = which.min(allMSEs2)
    if(bestI1 != bestI2) {
      browser()
      allEsts = sapply(allPs, function(p) {
        getGREGCV(LOOCVs, IWs^p, controlVarMat, controlVarMeans, normalizeWeights=TRUE, returnInfo=FALSE)
      })
    }
    allMSEs = allMSEs1
  } else {
    allMSEs = sapply(allPs, estMSE)
  }
  
  bestI = which.min(allMSEs)
  pEst = allPs[bestI]
  MSEp = allMSEs[bestI]
  
  # out = optim(0, estMSE)
  # pEst = expit(out$par)
  # 
  # MSEp = estMSE(logit(pEst))
  
  list(p=pEst, MSEest=MSEp, allPs=allPs, allMSEs=allMSEs)
}

getRatioIWCV = function(scores, IWs=1, controlVarMat=IWs, controlVarMeans=NA, tol = 1e-08, 
                         printVarRatio=FALSE, printRoot="") {
  
  
  Ys = scores * IWs
  Ybar = mean(Ys)
  
  if(!is.matrix(controlVarMat)) {
    if(is.na(controlVarMeans)) {
      controlVarMeans = mean(controlVarMat)
    }
    
    Xbar = controlVarMeans
    R = Ybar/Xbar
    
    
  }
}

# normalizeIWs: if TRUE, calculates YhatW = (1/sum(IWs)) * sum(LOOCVs * IWs). 
#               if FALSE, calculates YhatW = sum(LOOCVs * IWs).
# resids: used for variance estimation for GREG estimators
getP = function(LOOCVs, IWs, logitPSeq=seq(logit(.001), logit(.999), l=1000), normalizeIWs=TRUE, 
                unbiasedEstimate=NULL, resids=LOOCVs) {
  k = length(LOOCVs)
  allPs = c(0, expit(logitPSeq), 1)
  VarY = var(LOOCVs)
  
  if(is.null(unbiasedEstimate)) {
    if(normalizeIWs) {
      unbiasedEstimate = (1/sum(IWs)) * sum(LOOCVs * IWs)
    } else {
      unbiasedEstimate = sum(LOOCVs * IWs)
    }
  }
  
  # commented out because it assumes the weights are fixed
  # YhatW = ifelse(normalizeIWs, (1/sum(IWs)) * sum(LOOCVs * IWs), sum(LOOCVs * IWs))
  # estMSE = function(p=NULL, logitP=NULL) {
  #   if(is.null(p)) {
  #     p = expit(logitP)
  #   }
  #   YhatP = (1/sum(IWs^p)) * sum(LOOCVs * IWs^p)
  #   Bp2 = (YhatP - YhatW)^2
  # 
  #   Bp2 + VarY * (sum(IWs^(2*p))/sum(IWs^p)^2)
  # }
  
  # accounts for randomness in the weights and correlation with LOOCVs
  estMSE = function(p=NULL, logitP=NULL) {
    if(is.null(p)) {
      p = expit(logitP)
    }
    
    if(normalizeIWs) {
      ws = IWs^p/sum(IWs^p)
    } else {
      ws = IWs^p
    }
    
    YhatP = sum(LOOCVs * ws)
    Bp2 = (YhatP - unbiasedEstimate)^2
    
    VarYp = var(resids * ws)
    
    Bp2 + VarYp * k
  }
  
  allMSEs = sapply(allPs, estMSE)
  bestI = which.min(allMSEs)
  pEst = allPs[bestI]
  MSEp = allMSEs[bestI]
  
  # out = optim(0, estMSE)
  # pEst = expit(out$par)
  # 
  # MSEp = estMSE(logit(pEst))
  
  list(p=pEst, MSEest=MSEp, allPs=allPs, allMSEs=allMSEs)
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
                             preferential=FALSE, alpha=0, beta=1, n2=n) {
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
                                              n1=n, n2=n2, gridNs=gridNs, allSeeds=seeds, 
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
        out = load(paste0("savedOutput/griddedCVtest/n1", n, "_n2", n2, "_iter", i, unifText, prefText, twoDatText, ".RData"))
        thisList = list(trueMSE=trueMSE, wrongMSE=wrongMSE, LOOCVs=LOOCVs, LOOCV=LOOCV, 
                        griddedCVs=griddedCVs, gridNs=gridNs, iter=iter, 
                        rGRFargsTruth=rGRFargsTruth, rGRFargsSample=rGRFargsSample, rGRFargsWrong=rGRFargsWrong, 
                        n=n, nx=nx, ny=ny, sigmaEpsSq=sigmaEpsSq, allSeeds=allSeeds, 
                        alpha=alpha, beta=beta)
      } else {
        out = load(paste0("savedOutput/griddedCVtest2Datasets/n1", n, "_n2", n2, "_rho", rho, "_iter", i, ".RData"))
        # thisList = list(trueMSE=trueMSE, wrongMSE=wrongMSE, LOOCVs=LOOCVs, LOOCV=LOOCV, 
        #                 griddedCVs=griddedCVs, gridNs=gridNs, iter=iter, 
        #                 rGRFargsTruth=rGRFargsTruth, rGRFargsSample=rGRFargsSample, rGRFargsWrong=rGRFargsWrong, 
        #                 rho=rho, n1=n1, n2=n2, sigmaEpsSq1=sigmaEpsSq1, sigmaEpsSq2=sigmaEpsSq2, 
        #                 alpha=alpha, beta=beta, nx=nx, ny=ny, allSeeds=allSeeds)\
        thisList = list(trueMSE=trueMSE, wrongMSE=wrongMSE, LOOCVs=LOOCVs, LOOCVsWrong=LOOCVsWrong, 
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
  
  trueVarW1s = sapply(results, getName, thisName="trueVarW1")
  trueVarW2s = sapply(results, getName, thisName="trueVarW2")
  meanIW1s = sapply(results, getName, thisName="meanIW1")
  meanIW2s = sapply(results, getName, thisName="meanIW2")
  meanCVCIW1s = sapply(results, getName, thisName="meanCVCIW1")
  meanCVCIW2s = sapply(results, getName, thisName="meanCVCIW2")
  
  cors1IS = sapply(results, getName, thisName="cors1IS")
  cors2IS = sapply(results, getName, thisName="cors2IS")
  cors1CVC = sapply(results, getName, thisName="cors1CVC")
  cors2CVC = sapply(results, getName, thisName="cors2CVC")
  
  trueMSEs = sapply(results, getName, thisName="trueMSE")
  wrongMSEs = sapply(results, getName, thisName="wrongMSE")
  LOOCVs = sapply(results, getName, thisName="LOOCV")
  LOOISCVs = sapply(results, getName, thisName="LOOISCV")
  LOOVCCVs = sapply(results, getName, thisName="LOOVCCV")
  LOOCVCCVs = sapply(results, getName, thisName="LOOCVCCV")
  griddedCVs = sapply(results, getName, thisName="griddedCVs")
  
  LOOCVsWrong = sapply(results, getName, thisName="LOOCVWrong")
  LOOISCVsWrong = sapply(results, getName, thisName="LOOISCVWrong")
  LOOVCCVsWrong = sapply(results, getName, thisName="LOOVCCVWrong")
  LOOCVCCVsWrong = sapply(results, getName, thisName="LOOCVCCVWrong")
  griddedCVsWrong = sapply(results, getName, thisName="griddedCVsWrong")
  
  LOORCVs = sapply(results, getName, thisName="LOORCV")
  LOOISRCVs = sapply(results, getName, thisName="LOOISRCV")
  LOOVCRCVs = sapply(results, getName, thisName="LOOVCRCV")
  LOOCVCRCVs = sapply(results, getName, thisName="LOOCVCRCV")
  griddedRCVs = sapply(results, getName, thisName="griddedRCVs")
  LOOISRCVsWrong = sapply(results, getName, thisName="LOOISRCVWrong")
  LOOVCRCVsWrong = sapply(results, getName, thisName="LOOVCRCVWrong")
  LOOCVCRCVsWrong = sapply(results, getName, thisName="LOOCVCRCVWrong")
  griddedRCVsWrong = sapply(results, getName, thisName="griddedRCVsWrong")
  
  LOOR2CVs = sapply(results, getName, thisName="LOOR2CV")
  LOOISR2CVs = sapply(results, getName, thisName="LOOISR2CV")
  LOOVCR2CVs = sapply(results, getName, thisName="LOOVCR2CV")
  LOOCVCR2CVs = sapply(results, getName, thisName="LOOCVCR2CV")
  griddedR2CVs = sapply(results, getName, thisName="griddedR2CVs")
  LOOISR2CVsWrong = sapply(results, getName, thisName="LOOISR2CVWrong")
  LOOVCR2CVsWrong = sapply(results, getName, thisName="LOOVCR2CVWrong")
  LOOCVCR2CVsWrong = sapply(results, getName, thisName="LOOCVCR2CVWrong")
  griddedR2CVsWrong = sapply(results, getName, thisName="griddedR2CVsWrong")
  
  LOOISPCVs = sapply(results, getName, thisName="LOOISPCV")
  LOOISPRCVs = sapply(results, getName, thisName="LOOISPRCV")
  LOOISPR2CVs = sapply(results, getName, thisName="LOOISPR2CV")
  LOOISPCVsWrong = sapply(results, getName, thisName="LOOISPCVWrong")
  LOOISPRCVsWrong = sapply(results, getName, thisName="LOOISPRCVWrong")
  LOOISPR2CVsWrong = sapply(results, getName, thisName="LOOISPR2CVWrong")
  
  # calculate error
  LOOCVerrs = LOOCVs - trueMSEs
  LOOISCVerrs = LOOISCVs - trueMSEs
  LOOVCCVerrs = LOOVCCVs - trueMSEs
  LOOCVCCVerrs = LOOCVCCVs - trueMSEs
  griddedCVerrs = sweep(griddedCVs, 2, trueMSEs)
  griddedCVmeanErr = rowMeans(griddedCVerrs)
  LOOISRCVerrs = LOOISRCVs - trueMSEs
  LOOVCRCVerrs = LOOVCRCVs - trueMSEs
  LOOCVCRCVerrs = LOOCVCRCVs - trueMSEs
  griddedRCVerrs = sweep(griddedRCVs, 2, trueMSEs)
  LOOISR2CVerrs = LOOISR2CVs - trueMSEs
  LOOVCR2CVerrs = LOOVCR2CVs - trueMSEs
  LOOCVCR2CVerrs = LOOCVCR2CVs - trueMSEs
  griddedR2CVerrs = sweep(griddedR2CVs, 2, trueMSEs)
  
  LOOISPCVerrs = LOOISPCVs - trueMSEs
  LOOISPRCVerrs = LOOISPRCVs - trueMSEs
  LOOISPR2CVerrs = LOOISPR2CVs - trueMSEs
  
  # calculate percent error
  LOOCVpctErrs = 100 * LOOCVerrs/trueMSEs
  LOOISCVpctErrs = 100 * LOOISCVerrs/trueMSEs
  LOOVCCVpctErrs = 100 * LOOVCCVerrs/trueMSEs
  LOOCVCCVpctErrs = 100 * LOOCVCCVerrs/trueMSEs
  griddedCVpctErrs = sweep(griddedCVerrs, 2, 100/trueMSEs, "*")
  griddedCVmeanPctErr = rowMeans(griddedCVpctErrs)
  LOOISRCVpctErrs = 100 * LOOISRCVerrs/trueMSEs
  LOOCVCRCVpctErrs = 100 * LOOCVCRCVerrs/trueMSEs
  griddedRCVpctErrs = sweep(griddedRCVerrs, 2, 100/trueMSEs, "*")
  LOOISR2CVpctErrs = 100 * LOOISR2CVerrs/trueMSEs
  LOOCVCR2CVpctErrs = 100 * LOOCVCR2CVerrs/trueMSEs
  griddedR2CVpctErrs = sweep(griddedR2CVerrs, 2, 100/trueMSEs, "*")
  
  # calculate relative error
  LOOCVrelErrs = LOOCVs/trueMSEs
  LOOISCVrelErrs = LOOISCVs/trueMSEs
  LOOVCCVrelErrs = LOOVCCVs/trueMSEs
  LOOCVCCVrelErrs = LOOCVCCVs/trueMSEs
  griddedCVrelErrs = sweep(griddedCVs, 2, trueMSEs, "/")
  griddedCVmeanRelErr = rowMeans(griddedCVrelErrs)
  LOOISRCVrelErrs = LOOISRCVs/trueMSEs
  LOOCVCRCVrelErrs = LOOCVCRCVs/trueMSEs
  griddedRCVrelErrs = sweep(griddedRCVs, 2, trueMSEs, "/")
  LOOISR2CVrelErrs = LOOISR2CVs/trueMSEs
  LOOCVCR2CVrelErrs = LOOCVCR2CVs/trueMSEs
  griddedR2CVrelErrs = sweep(griddedR2CVs, 2, trueMSEs, "/")
  
  # calculate proportion of time correct model is selected
  LOOCVprop = mean(LOOCVs < LOOCVsWrong)
  LOOISCVprop = mean(LOOISCVs < LOOISCVsWrong)
  LOOVCCVprop = mean(LOOVCCVs < LOOVCCVsWrong)
  LOOCVCCVprop = mean(LOOCVCCVs < LOOCVCCVsWrong)
  griddedCVprop = rowMeans(griddedCVs < griddedCVsWrong)
  LOOISRCVprop = mean(LOOISRCVs < LOOISRCVsWrong)
  LOOCVCRCVprop = mean(LOOCVCRCVs < LOOCVCRCVsWrong)
  griddedRCVprop = rowMeans(griddedRCVs < griddedRCVsWrong)
  LOOISR2CVprop = mean(LOOISR2CVs < LOOISR2CVsWrong)
  LOOCVCR2CVprop = mean(LOOCVCR2CVs < LOOCVCR2CVsWrong)
  griddedR2CVprop = rowMeans(griddedR2CVs < griddedR2CVsWrong)
  LOOISPCVprop = mean(LOOISPCVs < LOOISPCVsWrong)
  LOOISPRCVprop = mean(LOOISPRCVs < LOOISPRCVsWrong)
  LOOISPR2CVprop = mean(LOOISPR2CVs < LOOISPR2CVsWrong)
  
  LOOCVpropMOE = qnorm(.975) * sqrt(LOOCVprop*(1-LOOCVprop) / niter)
  LOOISCVpropMOE = qnorm(.975) * sqrt(LOOISCVprop*(1-LOOISCVprop) / niter)
  LOOVCCVpropMOE = qnorm(.975) * sqrt(LOOVCCVprop*(1-LOOVCCVprop) / niter)
  LOOCVCCVpropMOE = qnorm(.975) * sqrt(LOOCVCCVprop*(1-LOOCVCCVprop) / niter)
  griddedCVpropMOE = qnorm(.975) * sqrt(griddedCVprop*(1-griddedCVprop) / niter)
  LOOISRCVpropMOE = qnorm(.975) * sqrt(LOOISRCVprop*(1-LOOISRCVprop) / niter)
  LOOCVCRCVpropMOE = qnorm(.975) * sqrt(LOOCVCRCVprop*(1-LOOCVCRCVprop) / niter)
  griddedRCVpropMOE = qnorm(.975) * sqrt(griddedRCVprop*(1-griddedRCVprop) / niter)
  LOOISR2CVpropMOE = qnorm(.975) * sqrt(LOOISR2CVprop*(1-LOOISR2CVprop) / niter)
  LOOCVCR2CVpropMOE = qnorm(.975) * sqrt(LOOCVCR2CVprop*(1-LOOCVCR2CVprop) / niter)
  griddedR2CVpropMOE = qnorm(.975) * sqrt(griddedR2CVprop*(1-griddedR2CVprop) / niter)
  LOOISPCVpropMOE = qnorm(.975) * sqrt(LOOISPCVprop*(1-LOOISPCVprop) / niter)
  LOOISPRCVpropMOE = qnorm(.975) * sqrt(LOOISPRCVprop*(1-LOOISPRCVprop) / niter)
  LOOISPR2CVpropMOE = qnorm(.975) * sqrt(LOOISPR2CVprop*(1-LOOISPR2CVprop) / niter)
  
  LOOCVpropHigh = LOOCVprop + LOOCVpropMOE
  LOOISCVpropHigh = LOOISCVprop + LOOISCVpropMOE
  LOOVCCVpropHigh = LOOVCCVprop + LOOVCCVpropMOE
  LOOCVCCVpropHigh = LOOCVCCVprop + LOOCVCCVpropMOE
  griddedCVpropHigh = griddedCVprop + griddedCVpropMOE
  LOOISRCVpropHigh = LOOISRCVprop + LOOISRCVpropMOE
  LOOCVCRCVpropHigh = LOOCVCRCVprop + LOOCVCRCVpropMOE
  griddedRCVpropHigh = griddedRCVprop + griddedRCVpropMOE
  LOOISR2CVpropHigh = LOOISR2CVprop + LOOISR2CVpropMOE
  LOOCVCR2CVpropHigh = LOOCVCR2CVprop + LOOCVCR2CVpropMOE
  griddedR2CVpropHigh = griddedR2CVprop + griddedR2CVpropMOE
  LOOISPCVpropHigh = LOOISPCVprop + LOOISPCVpropMOE
  LOOISCVpropHigh = LOOISCVprop + LOOISCVpropMOE
  LOOISR2CVpropHigh = LOOISR2CVprop + LOOISR2CVpropMOE
  
  LOOCVpropLow = LOOCVprop - LOOCVpropMOE
  LOOISCVpropLow = LOOISCVprop - LOOISCVpropMOE
  LOOVCCVpropLow = LOOVCCVprop - LOOVCCVpropMOE
  LOOCVCCVpropLow = LOOCVCCVprop - LOOCVCCVpropMOE
  griddedCVpropLow = griddedCVprop - griddedCVpropMOE
  LOOISRCVpropLow = LOOISRCVprop - LOOISRCVpropMOE
  LOOCVCRCVpropLow = LOOCVCRCVprop - LOOCVCRCVpropMOE
  griddedRCVpropLow = griddedRCVprop - griddedRCVpropMOE
  LOOISR2CVpropLow = LOOISR2CVprop - LOOISR2CVpropMOE
  LOOCVCR2CVpropLow = LOOCVCR2CVprop - LOOCVCR2CVpropMOE
  griddedR2CVpropLow = griddedR2CVprop - griddedR2CVpropMOE
  LOOISPCVpropLow = LOOISPCVprop - LOOISPCVpropMOE
  LOOISCVpropLow = LOOISCVprop - LOOISCVpropMOE
  # LOOISR2CVpropLow = LOOISR2CVprop - LOOISRR2CVpropMOE
  
  figureFolder = ifelse(twoDatasets, "figures/twoDatasetsTest/", "figures/gridTest/")
  
  # plot results ----
  # proportion of time selecting right model
  pdf(paste0(figureFolder, "selectProb_n1", n, "_n2", n2, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
  ylim = range(c(LOOCVpropHigh, LOOISCVpropHigh, LOOISRCVpropHigh, LOOVCCVpropHigh, LOOCVCCVpropHigh, LOOCVCRCVpropHigh, griddedCVpropHigh, 
                 LOOCVpropLow, LOOISCVpropLow, LOOISRCVpropLow, LOOVCCVpropLow, LOOCVCCVpropLow, LOOCVCRCVpropLow, griddedCVpropLow))
  plot(gridNs, griddedCVprop, type="n", log="x", axes=FALSE, 
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
  LOORcol = do.call("rgb", as.list(c(col2rgb("purple"))/255 * .8))
  # abline(h=LOORCVprop, lty=1, col="purple")
  # abline(h=LOORCVpropHigh, lty=2, col="purple")
  # abline(h=LOORCVpropLow, lty=2, col="purple")
  abline(h=LOOISCVprop, lty=1, col="orange")
  abline(h=LOOISCVpropHigh, lty=2, col="orange")
  abline(h=LOOISCVpropLow, lty=2, col="orange")
  LOOISRcol = do.call("rgb", as.list(c(col2rgb("orange"))/255 * .8))
  abline(h=LOOISPCVprop, lty=1, col=LOOISRcol)
  abline(h=LOOISPCVpropHigh, lty=2, col=LOOISRcol)
  abline(h=LOOISPCVpropLow, lty=2, col=LOOISRcol)
  abline(h=LOOVCCVprop, lty=1, col="brown")
  abline(h=LOOVCCVpropHigh, lty=2, col="brown")
  abline(h=LOOVCCVpropLow, lty=2, col="brown")
  LOOCVCRcol = do.call("rgb", as.list(c(col2rgb("brown"))/255 * .8))
  abline(h=LOOCVCCVprop, lty=1, col=LOOCVCRcol)
  abline(h=LOOCVCCVpropHigh, lty=2, col=LOOCVCRcol)
  abline(h=LOOCVCCVpropLow, lty=2, col=LOOCVCRcol)
  abline(h=0, lty=2, col="green")
  legend("topright", c("Gridded", "LOO", "LOOVC", "LOOCVC", "LOOIS", "LOOISP", "Truth"), col=c("blue", "purple", "brown", LOOCVCRcol, "orange", LOOISRcol, "green"), 
         pch=c(19, NA, NA, NA, NA, NA, NA), lty=1)
  dev.off()
  
  # bias
  pdf(paste0(figureFolder, "griddedBias_n1", n, "_n2", n2, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
  plot(gridNs, griddedCVmeanErr, type="n", log="x", axes=FALSE, 
       ylim=c(min(c(griddedCVmeanErr, mean(LOOCVerrs))), max(griddedCVmeanErr)), 
       xlab="Blocks per side", 
       ylab="MSE Bias", main=paste0("Gridded CV vs resolution (n=", n, unifTitleText, prefTitleText, twoDatTitleText, ")"))
  axis(side=1, at=gridNs)
  axis(side=2)
  box()
  lines(gridNs, griddedCVmeanErr, type="o", pch=19, col="blue")
  abline(h=mean(LOOCVerrs), lty=2, col="purple")
  abline(h=mean(LOORCVs), lty=2, col=LOORcol)
  abline(h=mean(LOOISCVerrs), lty=2, col="orange")
  abline(h=mean(LOOISPCVerrs), lty=2, col=LOOISRcol)
  abline(h=mean(LOOVCCVerrs), lty=2, col="brown")
  abline(h=mean(LOOCVCCVerrs), lty=2, col=LOOCVCRcol)
  abline(h=0, lty=2, col="green")
  legend("topright", c("Gridded", "LOO", "LOOR", "LOOVC", "LOOCVC", "LOOIS", "LOOISP", "Truth"), 
         col=c("blue", "purple", LOORcol, "brown", LOOCVCRcol, "orange", LOOISRcol, "green"), 
         pch=c(19, NA, NA, NA, NA, NA, NA, NA), lty=c(1, 2, 2, 2, 2, 2, 2, 2))
  dev.off()
  
  pdf(paste0(figureFolder, "griddedBiasExtrap_n1", n, "_n2", n2, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
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
  pdf(paste0(figureFolder, "griddedPctBias_n1", n, "_n2", n2, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
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
  
  pdf(paste0(figureFolder, "griddedPctBiasExtrap_n1", n, "_n2", n2, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
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
        relTicks2 = c(.67, .75, .9, 1, 1.1, 1.25, 1.5, 2, 3, 5)
        relTickLabs2 = c("0.67", "0.75", "0.9", "1", "1.1", "1.25", "1.5", "2", "3", "5")
      } else if(!twoDatasets) {
        relTicks1 = c(.67, 1, 2, 3, 4)
        relTickLabs1 = c("0.67", "1", "2", "3", "4")
        relTicks2 = c(.67, 1, 2, 4, 6, 8)
        relTickLabs2 = c("0.67", "1", "2", "4", "6", "8")
      } else {
        if(rho == -.8) {
          relTicks1 = seq(.7, 1.5, by=.1)
          relTickLabs1 = as.character(relTicks1)
          relTicks2 = c(.67, .75, .9, 1, 1.1, 1.25, 1.5, 2, 3, 5)
          relTickLabs2 = c("0.67", "0.75", "0.9", "1", "1.1", "1.25", "1.5", "2", "3", "5")
        } else if(rho == 0) {
          relTicks1 = seq(.7, 3, by=.1)
          relTickLabs1 = as.character(relTicks1)
          relTicks2 = c(.67, 1, 2, 3, 5)
          relTickLabs2 = c("0.67", "1", "2", "3", "5")
        } else if(rho == .8) {
          relTicks1 = seq(.7, 1.5, by=.1)
          relTickLabs1 = as.character(relTicks1)
          relTicks2 = c(.67, .75, .9, 1, 1.1, 1.25, 1.5, 2, 3, 5)
          relTickLabs2 = c("0.67", "0.75", "0.9", "1", "1.1", "1.25", "1.5", "2", "3", "5")
        }
        
      }
    } else if(n == 500) {
      if(!unif && !twoDatasets) {
        relTicks = c(.6, 1, 2, 5, 10, 20, 50, 100)
        relTickLabs = c("0.6", "1", "2", "5", "10", "20", "50", "100")
        relTicks1 = relTicks2 = relTicks
        relTickLabs1 = relTickLabs2 = relTickLabs
      } else if(!twoDatasets) {
        browser()
      } else {
        relTicks1 = c(.8, .9, 1, 1.1, 1.2, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
        relTickLabs1 = as.character(relTicks1)
        relTicks2 = c(.67, 1, 2, 3, 5, 10, 20, 30, 40, 50, 60)
        relTickLabs2 = c("0.67", "1", "2", "3", "5", "10", "20", "30", "40", "50", "60")
      }
    } else {
      relTicks = c(.67, 1, 2, 5, 10, 20, 50, 100)
      relTickLabs = c("0.67", "1", "2", "5", "10", "20", "50", "100")
      relTicks1 = relTicks2 = relTicks
      relTickLabs1 = relTickLabs2 = relTickLabs
    }
  }
  
  # browser()
  
  pdf(paste0(figureFolder, "griddedRelBias_n1", n, "_n2", n2, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
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
  
  pdf(paste0(figureFolder, "griddedRelBiasExtrap_n1", n, "_n2", n2, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
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
  
  # Overall MSE ----
  methods = c("LOO", 
              "LOO-IW", "LOO-IWR", "LOO-IWR2", 
              "LOO-IWP", "LOO-IWRP", "LOO-IWRP2", 
              "LOO-VC", "LOO-VCR", "LOO-VCR2"
              # "LOO-VCP"
              # "LOO-CVC", "LOO-CVCP", "LOO-CVCR", "LOO-CVCR2"
  )
  # methods = methods[-c(match(c("LOO-VCP"), methods))]
  pdf(paste0(figureFolder, "CVMSE_n1", n, "_n2", n2, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=8, height=6)
  dat = data.frame(
    Method=factor(rep(methods, each=niter), 
                  levels=methods, ordered=FALSE), 
    MSE=c((LOOCVs - trueMSEs)^2, 
          (LOOISCVs - trueMSEs)^2, 
          (LOOISRCVs - trueMSEs)^2, 
          (LOOISR2CVs - trueMSEs)^2, 
          (LOOISPCVs - trueMSEs)^2, 
          (LOOISPRCVs - trueMSEs)^2, 
          (LOOISPR2CVs - trueMSEs)^2, 
          (LOOVCCVs - trueMSEs)^2, 
          (LOOVCRCVs - trueMSEs)^2, 
          (LOOVCR2CVs - trueMSEs)^2
          # (LOOVCPCVs - trueMSEs)^2, 
          # (LOOCVCCVs - trueMSEs)^2, 
          # (LOOCVCPCVs - trueMSEs)^2,
          # (LOOCVCRCVs - trueMSEs)^2, 
          # (LOOCVCR2CVs - trueMSEs)^2)
    ))
  boxplot(MSE~Method, data=dat, col="skyblue", log="y", ylab="Estimator Sq. Err.")
  dev.off()
  
  pdf(paste0(figureFolder, "CVMSEerr_n1", n, "_n2", n2, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=8, height=6)
  datac <- summarySEwithin(dat, measurevar="MSE", withinvars=c("Method"))
  #>    Shape   ColorScheme  N     Time Time_norm       sd        se        ci
  #> 1  Round       Colored 12 43.58333  43.58333 1.212311 0.3499639 0.7702654
  #> 2  Round Monochromatic 12 44.58333  44.58333 1.331438 0.3843531 0.8459554
  #> 3 Square       Colored 12 42.58333  42.58333 1.461630 0.4219364 0.9286757
  #> 4 Square Monochromatic 12 43.58333  43.58333 1.261312 0.3641095 0.8013997
  ggplot(datac, aes(x=Method, y=MSE)) +
    geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=MSE-ci, ymax=MSE+ci)) +
    # coord_cartesian(ylim=c(40,46)) +
    # scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
    # scale_y_continuous(breaks=seq(1:100)) +
    theme_bw()
    # geom_hline(yintercept=38)
  dev.off()
  
  pdf(paste0(figureFolder, "CVMRE_n1", n, "_n2", n2, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=6, height=6)
  dat = data.frame(
    Method=rep(c("LOO", "LOO-IW", "LOO-VC", "LOO-ISR", "LOO-CVCR"), each=n), 
    sqResids=c(LOOCVs/trueMSEs, 
               LOOISCVs/trueMSEs, 
               LOOVCCVs/trueMSEs, 
               LOOISRCVs/trueMSEs, 
               LOOCVCRCVs/trueMSEs))
  boxplot(sqResids~Method, data=dat, log="y", col="skyblue", ylab="Estimator Rel. Err.")
  dev.off()
  
  pdf(paste0(figureFolder, "CVMSE_n1", n, "_n2", n2, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=8, height=6)
  dat = data.frame(
    Method=factor(rep(c("LOO", "LOO-IW", "LOO-IWP", "LOO-IWPR2", "LOO-VC", "LOO-VCR", "LOO-VCR2", "LOO-CVC", "LOO-CVCR", "LOO-CVCR2"), each=niter), 
                  levels=c("LOO", "LOO-IW", "LOO-IWP", "LOO-IWPR2", "LOO-VC", "LOO-VCR", "LOO-VCR2", "LOO-CVC", "LOO-CVCR", "LOO-CVCR2"), ordered=FALSE), 
    MSE=c((LOOCVs - trueMSEs)^2, 
          (LOOISCVs - trueMSEs)^2, 
          (LOOISPCVs - trueMSEs)^2, 
          (LOOISPR2CVs - trueMSEs)^2, 
          (LOOVCCVs - trueMSEs)^2, 
          (LOOVCRCVs - trueMSEs)^2, 
          (LOOVCR2CVs - trueMSEs)^2, 
          (LOOCVCCVs - trueMSEs)^2, 
          (LOOCVCRCVs - trueMSEs)^2, 
          (LOOCVCR2CVs - trueMSEs)^2))
  boxplot(MSE~Method, data=dat, col="skyblue", log="y", ylab="Estimator Sq. Err.")
  dev.off()
  
  pdf(paste0(figureFolder, "selProbErr_n1", n, "_n2", n2, unifText, prefText, twoDatText, "_niter", niter, ".pdf"), width=8, height=6)
  dat = data.frame(
    Method=factor(rep(c("LOO", "LOO-IW", "LOO-IWP", "LOO-IWPR2", "LOO-VC", "LOO-VCR", "LOO-VCR2", "LOO-CVC", "LOO-CVCR", "LOO-CVCR2"), each=niter), 
                  levels=c("LOO", "LOO-IW", "LOO-IWP", "LOO-IWPR2", "LOO-VC", "LOO-VCR", "LOO-VCR2", "LOO-CVC", "LOO-CVCR", "LOO-CVCR2"), ordered=FALSE), 
    Probability=c((LOOCVs - LOOCVsWrong < 0), 
          (LOOISCVs - LOOISCVsWrong < 0), 
          (LOOISPCVs - LOOISPCVsWrong < 0), 
          (LOOISPR2CVs - LOOISPR2CVsWrong < 0), 
          (LOOVCCVs - LOOVCCVsWrong < 0), 
          (LOOVCRCVs - LOOVCRCVsWrong < 0), 
          (LOOVCR2CVs - LOOVCR2CVsWrong < 0), 
          (LOOCVCCVs - LOOCVCCVsWrong < 0), 
          (LOOCVCRCVs - LOOCVCRCVsWrong < 0), 
          (LOOCVCR2CVs - LOOCVCR2CVsWrong < 0)))
  datac <- summarySEwithin(dat, measurevar="Probability", withinvars=c("Method"))
  #>    Shape   ColorScheme  N     Time Time_norm       sd        se        ci
  #> 1  Round       Colored 12 43.58333  43.58333 1.212311 0.3499639 0.7702654
  #> 2  Round Monochromatic 12 44.58333  44.58333 1.331438 0.3843531 0.8459554
  #> 3 Square       Colored 12 42.58333  42.58333 1.461630 0.4219364 0.9286757
  #> 4 Square Monochromatic 12 43.58333  43.58333 1.261312 0.3641095 0.8013997
  ggplot(datac, aes(x=Method, y=Probability)) +
    geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=Probability-ci, ymax=Probability+ci)) +
    # coord_cartesian(ylim=c(40,46)) +
    # scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
    # scale_y_continuous(breaks=seq(1:100)) +
    theme_bw()
  # geom_hline(yintercept=38)
  dev.off()
  
  browser()
  
  print(paste0("MSE of LOO: ", mean((LOOCVs - trueMSEs)^2)))
  print(paste0("MSE of LOOIS: ", mean((LOOISCVs - trueMSEs)^2)))
  print(paste0("MSE of LOOISR: ", mean((LOOISRCVs - trueMSEs)^2)))
  print(paste0("MSE of LOOISR: ", mean((LOOISR2CVs - trueMSEs)^2)))
  print(paste0("MSE of LOOISP: ", mean((LOOISPCVs - trueMSEs)^2)))
  print(paste0("MSE of LOOISPR: ", mean((LOOISPRCVs - trueMSEs)^2)))
  print(paste0("MSE of LOOISPR: ", mean((LOOISPR2CVs - trueMSEs)^2)))
  print(paste0("MSE of LOOVC: ", mean((LOOVCCVs - trueMSEs)^2)))
  print(paste0("MSE of LOOCVC: ", mean((LOOCVCCVs - trueMSEs)^2)))
  print(paste0("MSE of LOOCVCR: ", mean((LOOCVCRCVs - trueMSEs)^2)))
  print(paste0("MSE of LOOCVCR: ", mean((LOOCVCR2CVs - trueMSEs)^2)))
  
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


