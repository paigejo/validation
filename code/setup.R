# setup script for the jittering project
# install required packages if necessary
if(FALSE) {
  install.packages("Matrix")
  install.packages("spam")
  install.packages("fields")
  install.packages("LatticeKrig")
  install.packages("invgamma")
  install.packages("latex2exp")
  install.packages("xtable")
  install.packages("profvis")
  install.packages("geosphere")
  install.packages("viridis")
  install.packages("sp")
  install.packages("raster")
  install.packages("MCMCpack")
  install.packages("numDeriv")
  install.packages("INLA")
  install.packages("edfun") # for drawing from empirical distributions quickly
  install.packages("data.table")
  install.packages("sampling")
  install.packages("haven")
  install.packages("survey")
  install.packages("abind")
  install.packages("devtools")
  install.packages("ggplot2")
  install.packages("WeightedCluster")
  install.packages("colorspace")
  install.packages("TMB")
  install.packages("plyr")
  install.packages("deldir")
  install.packages("stringr")
  install.packages("truncnorm")
  
  library(devtools)
  if(FALSE) {
    install_github("https://github.com/richardli/SUMMER/")
    install_github("https://github.com/paigejo/SUMMER/")
    install_local("~/git/SUMMER/")
    document("~/git/SUMMER/")
    load_all("~/git/SUMMER/")
  }
}

# load required packages and R scripts
library(Matrix)
library(spam)
library(fields)
library(LatticeKrig)
library(invgamma)
library(latex2exp)
library(xtable)
library(profvis)
library(geosphere)
library(viridis)
library(sp)
library(raster)
library(MCMCpack)
library(numDeriv)
library(INLA)
library(edfun) # for drawing from empirical distributions quickly
library(data.table)
library(sampling)
library(haven)
library(survey)
library(abind)
# install_github("https://github.com/richardli/SUMMER/tree/dev")
library(SUMMER)
# library(Rcpp)
library(ggplot2)
library(WeightedCluster)
library(devtools)
library(colorspace)
library(TMB)
library(plyr)
library(deldir)
library(stringr)
library(truncnorm)

inf = sessionInfo()
if(inf$platform == "x86_64-apple-darwin17.0 (64-bit)") {
  setwd("~/git/jittering/")
  options(error=recover)
} else if(inf$platform != "x86_64-apple-darwin15.6.0 (64-bit)" && inf$platform != "x86_64-w64-mingw32/x64 (64-bit)" && inf$platform != "x86_64-pc-linux-gnu (64-bit)") {
  # INLA:::inla.dynload.workaround() # use inla.binary.install instead
  # avoid setting too many threads and thereby using too much memory
  inla.setOption(num.threads=1)
  options(error=traceback)
  setwd("~/git/jittering/")
} else if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)" && inf$platform != "x86_64-pc-linux-gnu (64-bit)") {
  setwd("~/git/jittering/")
  options(error=recover)
} else {
  setwd("~/git/jittering/")
  inla.setOption(num.threads=1) # consider raising
  options(error=recover)
}
options(warn=1)

setwd("~/git/validation")
source("code/genPareto.R")
source('code/griddedTest.R')
source('code/sim.R')
source('code/scores.R')
source('code/sim2DatasetsTest.R')
source('code/simGridResolutionTest.R')
source('code/simNonstationaryErrorTest.R')
source('code/utilityFuns.R')
source('code/voronoi.R')
source('code/test.R')
source("code/genericSpatialPlottingFunctions.R")



