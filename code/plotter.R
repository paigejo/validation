# for plotting illustrations

# plot illustrations of blockCV and buffered CV
# n: number of observations to plot
# M: illustrate M x M grid block CV
plotSpatialCVIllustrations = function(seed=1, n=75, M=5, bufferRadius=.2) {
  set.seed(seed)
  
  # simulate observations
  xs = runif(n)
  ys = runif(n)
  
  # Buffered CV ----
  # select a test point for buffered CV
  possibleXs = (xs > bufferRadius) & (xs < (1-bufferRadius))
  possibleYs = (ys > bufferRadius) & (ys < (1-bufferRadius))
  possibleL = possibleXs & possibleYs
  testIBuff = sample(which(possibleL), 1)
  xsTestBuff = xs[testIBuff]
  ysTestBuff = ys[testIBuff]
  
  # find train points for buffered CV
  dists = rdist(cbind(xsTestBuff, ysTestBuff), cbind(xs, ys))
  trainL = dists > bufferRadius
  xsTrainBuff = xs[trainL]
  ysTrainBuff = ys[trainL]
  
  # find left out pounts for buffered CV
  leftOutI = (1:n)[-c(which(trainL), testIBuff)]
  xsLeftOutBuff = xs[leftOutI]
  ysLeftOutBuff = ys[leftOutI]
  
  # Make table
  tabBuff = data.frame(x=c(xsTrainBuff, xsLeftOutBuff, xsTestBuff), 
                       y=c(ysTrainBuff, ysLeftOutBuff, ysTestBuff), 
                       type=c(rep("Train", length(xsTrainBuff)), 
                              rep("Left out", length(xsLeftOutBuff)), 
                              rep("Test", length(xsTestBuff))))
  
  # Block CV ----
  # determine blocks
  lims = seq(0, 1, l=M+1)
  # middleI = ceiling(M/2)
  # lowLim = lims[middleI]
  # highLim = lims[middleI+1]
  lowLimX = rev(lims)[match(TRUE, rev(lims) < xsTestBuff)]
  highLimX = lims[match(TRUE, lims > xsTestBuff)]
  lowLimY = rev(lims)[match(TRUE, rev(lims) < ysTestBuff)]
  highLimY = lims[match(TRUE, lims > ysTestBuff)]
  
  # determine which points are left out of the block with the test point
  goodXs = (xs > lowLimX) & (xs < highLimX)
  goodYs = (ys > lowLimY) & (ys < highLimY)
  testL = goodXs & goodYs
  trainL = !testL
  xsTestBlock = xs[testL]
  ysTestBlock = ys[testL]
  xsTrainBlock = xs[!testL]
  ysTrainBlock = ys[!testL]
  
  tabBlock = data.frame(x=c(xsTrainBlock, xsTestBlock), 
                       y=c(ysTrainBlock, ysTestBlock), 
                       type=c(rep("Train", length(xsTrainBlock)), 
                              rep("Test", length(xsTestBlock))))
  
  # Plotting ----
  # plotting settings
  pchTrain = 19
  pchLeftOut = 1
  pchTest = 17
  colTrain = "blue"
  colLeftOut = "black"
  colTest = "red"
  cex = .5
  
  # plot Block CV illustration
  # pdf("blockCVillustration.pdf", width=5, height=5)
  # plot(1, 1, type="n", xlim=c(0,1), ylim=c(0,1))
  # points(xsTrainBlock, ysTrainBlock, col=colTrain, pch=pchTrain, cex=cex)
  # points(xsTestBlock, ysTestBlock, col=colTest, pch=pchTest, cex=cex)
  # legend()
  # dev.off()
  
  p = ggplot(tabBlock, aes(x=x, y=y, color=type, shape=type)) + 
    geom_point() + 
    scale_color_manual(name="", values=c(colTest, colTrain)) + 
    scale_shape_manual(name="", values=c(pchTest, pchTrain)) + 
    geom_hline(yintercept = lims, linetype="dashed") + 
    geom_vline(xintercept = lims, linetype="dashed") + 
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    scale_x_continuous(name=NULL, limits=c(0, 1), breaks=c(0,1)) +
    scale_y_continuous(name=NULL, limits=c(0, 1), breaks=c(0,1)) +
    theme_classic() + theme(plot.margin = unit(c(0.05, 0.05, 0, 0), "inches"), 
                            legend.position = "none", legend.margin=margin(c(0,0,0,-5)), legend.text=element_text(size=14))
  ggsave("figures/illustrations/blockCV.pdf", plot=p, width=5, height=5)
  
  # plot buffered CV illustration
  require(ggforce)
  p = ggplot(tabBuff, aes(x=x, y=y, color=type, shape=type)) + 
    geom_point() + 
    scale_color_manual(name="", values=c(colLeftOut, colTest, colTrain)) + 
    scale_shape_manual(name="", values=c(pchLeftOut, pchTest, pchTrain)) + 
    geom_circle(aes(x0=xsTestBuff, y0=ysTestBuff, r=bufferRadius), color="black", linetype="dashed", show.legend=FALSE) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    scale_x_continuous(name=NULL, limits=c(0, 1), breaks=c(0,1)) +
    scale_y_continuous(name=NULL, limits=c(0, 1), breaks=c(0,1)) +
    theme_classic() + theme(plot.margin = unit(c(-.05, 0.05, 0, 0), "inches"), 
                            legend.position = "top", legend.margin=margin(c(0,0,0,-5)), legend.text=element_text(size=14))
  ggsave("figures/illustrations/bufferedCV.pdf", plot=p, width=5, height=5.32)
  
  # Voronoi Diagram
  require(deldir)
  
  # p = ggplot(tabBuff, aes(x=x, y=y, group = -1L)) +
  #   geom_voronoi_tile(linetype = "dashed") +
  #   geom_voronoi_segment(linetype = "dashed") + # Type of the lines
  #   geom_point() +
  #   coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
  #   scale_x_continuous(name=NULL, limits=c(0, 1), breaks=lims) +
  #   scale_y_continuous(name=NULL, limits=c(0, 1), breaks=lims) +
  #   theme_classic()
  
  tabVC = tabBuff
  tabVC$type[tabVC$type == "Left out"] = "Train"
  p = ggplot(tabVC, aes(x=x, y=y, group = -1L)) +
    # geom_voronoi_segment(linetype = "dashed", bound = cbind(x=c(0,0,1,1), y=c(0,1,1,0))) + # Type of the lines
    geom_voronoi_tile(linetype = "dotted", bound=cbind(x=c(0,0,1,1), y=c(0,1,1,0)), 
                      alpha=1, colour="black", fill="white") + # Type of the lines
    geom_point(aes(colour=type, shape=type)) +
    scale_color_manual(name="", values=c(colTest, colTrain)) + 
    scale_shape_manual(name="", values=c(pchTest, pchTrain)) + 
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) + 
    scale_x_continuous(name=NULL, limits=c(0, 1), breaks=c(0,1)) +
    scale_y_continuous(name=NULL, limits=c(0, 1), breaks=c(0,1)) +
    theme_classic() + 
    theme(plot.margin = unit(c(0.05, 0.05, 0, 0), "inches"), legend.position = "none")
  ggsave("figures/illustrations/VCCVcolored.pdf", plot=p, width=5, height=5)
  
  p = ggplot(tabVC, aes(x=x, y=y, group = -1L)) +
    # geom_voronoi_segment(linetype = "dashed", bound = cbind(x=c(0,0,1,1), y=c(0,1,1,0))) + # Type of the lines
    geom_voronoi_tile(linetype = "dotted", bound=cbind(x=c(0,0,1,1), y=c(0,1,1,0)), 
                      alpha=1, colour="black", fill="white") + # Type of the lines
    geom_point(color="blue") +
    scale_color_manual(name="", values=c(colTest, colTrain)) + 
    scale_shape_manual(name="", values=c(pchTest, pchTrain)) + 
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) + 
    scale_x_continuous(name=NULL, limits=c(0, 1), breaks=c(0,1)) +
    scale_y_continuous(name=NULL, limits=c(0, 1), breaks=c(0,1)) +
    theme_classic() + theme(plot.margin = unit(c(0.05, 0.05, 0, 0), "inches"))
  ggsave("figures/illustrations/VCCV.pdf", plot=p, width=5, height=5)
  browser()
  invisible(NULL)
}

plotIllustration1D = function(seed=16, subsample=10000/500) {
  set.seed(seed)
  
  # construct and sample grid ----
  sIs = 1:10000
  sgrid = sIs/max(sIs)
  sIsamp = sort(sample(sIs, 10, prob=sgrid/sum(sgrid), replace = FALSE))
  ssamp = sgrid[sIsamp]
  sIrest = sIs[-sIsamp]
  sgridRest = sgrid[sIrest]
  
  # get distances, covs ----
  fullDistMat = rdist(cbind(sgrid))
  fullSigma = exp(-10*fullDistMat)
  Sigma = fullSigma[sIsamp,sIsamp]
  
  # simulate points ----
  L = t(chol(Sigma))
  ys = L %*% rnorm(10)
  
  # get nearest sampled locations to the right and left
  sampleItoRight = sapply(sIrest, function(i) {sIsamp[match(TRUE, sIsamp > i)]})
  sampleItoLeft = sapply(sIrest, function(i) {rev(sIsamp)[match(TRUE, rev(sIsamp) < i)]})
  # sToRight = sgrid[sampleItoRight]
  # sToLeft = sgrid[sampleItoLeft]
  # leftRightDists = cbind(fullDistMat[sIrest, sampleItoLeft], 
  #                        fullDistMat[sIrest, sampleItoRight])
  
  # calculate true MSEs ----
  # (there is an analytical solution, since Markov, but I'm lazy...)
  gridMSEs = rep(0, length(sgrid))
  for(tempI in 1:length(sIrest)) {
    # tempI is index in sIrest, i is index in sI
    i = sIrest[tempI]
    thisS = sgrid[i]
    
    if((i < sIsamp[1]) || (i > sIsamp[10])) {
      # only the nearest sample matters
      if(i < sIsamp[1]) {
        thisCov = fullSigma[i, sIsamp[1]]
      } else if(i > sIsamp[10]) {
        thisCov = fullSigma[i, sIsamp[10]]
      }
      gridMSEs[i] = 1 - thisCov^2
    } else {
      leftI = sampleItoLeft[tempI]
      rightI = sampleItoRight[tempI]
      
      leftCov = fullSigma[i, leftI]
      rightCov = fullSigma[i, rightI]
      leftRightCov = fullSigma[leftI, rightI]
      
      SigmaA = 1
      SigmaAB = rbind(c(leftCov, rightCov))
      SigmaB = rbind(c(1, leftRightCov), 
                     c(leftRightCov, 1))
      
      gridMSEs[i] = SigmaA - SigmaAB %*% solve(SigmaB) %*% t(SigmaAB)
      
      if(is.na(gridMSEs[i])) {
        browser()
      }
    }
  }
  
  # do LOOCV ----
  LOO = numeric(10)
  for(i in 1:10) {
    thisSsamp = ssamp[-i]
    thisYs = ys[-i]
    thisSigma = Sigma[-i,-i]
    pred = condMeanMVN(Sigma[i,i], matrix(Sigma[i,-i], nrow=1), Sigma[-i,-i], thisYs, 
                getFullCov=F, getCondVar=FALSE)$muAcondB
    LOO[i] = (ys[i] - pred)^2
  }
  
  meanMSE = mean(gridMSEs)
  meanLOO = mean(LOO)
  
  # plotting ----
  subSeq = seq(from=subsample, to=length(sgrid), by=subsample)
  df = data.frame(s=ssamp, y=ys)
  dfLines = data.frame(sgrid=c(sgrid[subSeq], rep(c(0,1), 2)), 
                       y=c(gridMSEs[subSeq], rep(meanMSE, 2), rep(meanLOO, 2)), 
                       group=factor(c(rep("MSE at location", length(sgrid[subSeq])), rep("Mean MSE over domain", 2), rep("LOOCV MSE", 2))))
  
  pl <- ggplot(df, aes(s, y)) + 
    # geom_hline(yintercept=meanMSE, linetype="dashed", color="blue") + 
    # geom_hline(yintercept=meanLOO, linetype="dashed", color="red") + 
    geom_line(data=dfLines, aes(x=sgrid, y=y, color=group, linetype=group)) + 
    geom_point() + 
    scale_linetype_manual(values=c("dashed", "dashed", "solid"))+
    scale_color_manual(values=c('red', 'blue','blue'))+
    scale_x_continuous(name="", expand=c(0, 0.01), breaks=c(0, 1), limits=c(0,1)) + 
    scale_y_continuous(name="", expand=c(0, 0.07)) + 
    theme_classic() + 
    theme(legend.position=c(0.01,1.07), legend.justification=c(0,1), legend.title=element_blank())
  pl
  ggsave("figures/illustrations/illustration1D.pdf", pl, height=2.7, width=4.5)
  browser()
}

# plot illustration of covariate shift and heteroscedasticity
plotCovShiftIllustration = function() {
  set.seed(8)
  
  # set parameters ----
  # y = a x + b + eps sqrt(x)
  minX = 0
  maxX = 3
  a = .5
  b = .2
  sigma = .2
  muTr = 1
  sigmaTr=.5
  muTe = 2
  sigmaTe=.5
  nTr = 100
  nTe = 100
  
  # simulate data ----
  require(truncnorm)
  # xTr = rnorm(nTr, muTr, sigmaTr)
  # xTe = rnorm(nTe, muTe, sigmaTe)
  xTr = rtruncnorm(nTr, a=minX, b=maxX, muTr, sigmaTr)
  xTe = rtruncnorm(nTe, a=minX, b=maxX, muTe, sigmaTe)
  yTr = rtruncnorm(nTr, a=0, b=Inf, a * xTr + b, sigma*xTr)
  yTe = rtruncnorm(nTe, a=0, b=Inf, a * xTe + b, sigma*xTe)
  epsTr = yTr - a * xTr + b
  epsTe = yTe - a * xTe + b
  
  # calculate squared error ----
  sqErrTr = epsTr^2
  sqErrTe = epsTe^2
  
  # construct data.frame ----
  df = data.frame(X=c(xTr, xTe), Y=c(yTr, yTe), sqErr=c(sqErrTr, sqErrTe), Group=factor(c(rep("Present", nTr), rep("Future", nTe)), levels=c("Present", "Future")))
  
  # plotting ----
  require(gridExtra)
  require(ggplot2)
  
  # scatter plot of x and y variables
  # color by groups
  scatterPlot <- ggplot(df, aes(X, Y, color=Group, shape=Group)) + 
    geom_point() + 
    scale_color_manual(values = c(rgb(0,0,1,.6), rgb(1,0,0,.6))) + 
    scale_shape_manual(values=c(16, 16)) + 
    geom_abline(slope=a, intercept=b) + 
    scale_x_continuous(name="Precipitation", expand=c(0, .03), breaks=NULL) + 
    scale_y_continuous(name="Mosquito Abundance", breaks=NULL, expand=c(0, .02)) + 
    theme_classic() + 
    guides() + 
    theme(legend.position=c(0.01,1), legend.justification=c(0,1), legend.title=element_blank())
  scatterPlot
  
  # Marginal density plot of x (top panel): 
  require(dplyr)
  labels <- data.frame(xPos=c(.8, 1.95), yPos=.92, Group=df$Group[c(1,nTr+nTe)])
  
  xdensity <- ggplot(df, aes(X, fill=Group)) + 
    geom_density(alpha=.5) + 
    geom_text(data=labels, aes(x=xPos, y=yPos, label=Group), hjust=0) + 
    scale_fill_manual(values = c('blue','red')) + 
    scale_y_continuous(name="Density", breaks = NULL, limits=c(0,1), expand=c(0, .005)) +
    scale_x_continuous(name=NULL, expand=c(0, .03), breaks=NULL) + 
    theme_classic() + 
    theme(legend.position = "none", legend.title=element_blank())
  xdensity
  
  # Marginal violin plot of y (right panel)
  yviolin <- ggplot(df, aes(x=Group, y=sqErr, fill=Group)) + 
    geom_violin(alpha=.5, width = 1) + 
    stat_summary(fun = "mean",
                 geom = "point") + 
    scale_fill_manual(values = c('blue','red')) +
    # scale_y_continuous(name="Square error", breaks = c(seq(0,.5,0.1), seq(.5, 2.5, .5)), trans = "log", expand = c(0,0.005), limits=c(.01, 2.2)) +
    # scale_y_continuous(name="Square error", breaks = c(seq(0,.25,0.05), seq(.5, 2.5, .5)), trans = "log", expand = c(0,0.005), limits=c(.01, 2.2)) +
    scale_y_continuous(name="Square Error", breaks=NULL, expand = c(0,0.005)) +
    scale_x_discrete(name=NULL) + 
    theme_classic() + 
    theme(legend.position = "none", legend.title=element_blank(), 
          # plot.margin = unit(c(0.1, 0, -.15, 0), "inches"), 
          axis.line=element_blank(), axis.text=element_text(size=11)
          )
  yviolin
  
  # put it all together (with a blank plot)
  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank()
    )
  finalFig = arrangeGrob(xdensity, blankPlot, scatterPlot, yviolin, 
               ncol=2, nrow=2, widths=c(3, 1.1), heights=c(.75, 2.4), 
               padding = unit(0.1, "line"))
  ggsave("figures/illustrations/covShiftIllustration.pdf", finalFig, 
         height=3.7, width=5.5)
  
  browser()
}






