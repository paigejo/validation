# for plotting illustrations

# plot illustrations of blockCV and buffered CV
# n: number of observations to plot
# M: illustrate M x M grid block CV
plotSpatialCVIllustrations = function(seed=1, n=75, M=5, bufferRadius=.2) {
  set.seed(seed)
  
  # simulate observations
  xs = runif(n)
  ys = runif(n)
  
  # Block CV ----
  # determine blocks
  lims = seq(0, 1, l=M+1)
  middleI = ceiling(M/2)
  lowLim = lims[middleI]
  highLim = lims[middleI+1]
  
  # determine which points are left out of the center block
  goodXs = (xs > lowLim) & (xs < highLim)
  goodYs = (ys > lowLim) & (ys < highLim)
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
    scale_x_continuous(name="X", limits=c(0, 1), breaks=lims) +
    scale_y_continuous(name="Y", limits=c(0, 1), breaks=lims) +
    theme_classic()
  ggsave("figures/illustrations/blockCV.pdf", plot=p, width=5.85, height=5)
  browser()
  # plot buffered CV illustration
  require(ggforce)
  p = ggplot(tabBuff, aes(x=x, y=y, color=type, shape=type)) + 
    geom_point() + 
    scale_color_manual(name="", values=c(colLeftOut, colTest, colTrain)) + 
    scale_shape_manual(name="", values=c(pchLeftOut, pchTest, pchTrain)) + 
    geom_circle(aes(x0=xsTestBuff, y0=ysTestBuff, r=bufferRadius), color="black", linetype="dashed", show.legend=FALSE) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    scale_x_continuous(name="X", limits=c(0, 1), breaks=lims) +
    scale_y_continuous(name="Y", limits=c(0, 1), breaks=lims) +
    theme_classic()
  ggsave("figures/illustrations/bufferedCV.pdf", plot=p, width=5.85, height=5)
  
  invisible(NULL)
}







