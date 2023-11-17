# srun --partition=CPUQ --time=02:00:00 --mem-per-cpu=10000 --pty bash

source("setup.R")
options(error=traceback)
index = as.numeric(commandArgs(trailingOnly = TRUE)) # test with index == something

# Rprof("savedOutput/simStudyResults/tempFiles/data.Rprof", interval = 0.01, line.profiling = TRUE,
#       gc.profiling = TRUE, memory.profiling = TRUE)

# p = profvis({

# griddedResTestIter(rGRFargsTruth=NULL, rGRFargsSample=NULL, 
#                    n=50, grindNs=2^(1:6), iter=1, seed=123, 
#                    nx=100, ny=100, sigmaEpsSq=0)
nSimPerRho=1000
rhos = c(-.8, -.4, 0, .4, .8)
ns = c(50, 500, 1000)
nSimPerRho = c(1000, 100, 50)
set.seed(123)
seeds = sample(1:10000000, sum(nSimPerRho))
if(index <= nSimPerRho[1]*length(rhos)) {
  thisN = ns[1]
  rho = rhos[ceiling(index/nSimPerRho[1])]
} else if((index - nSimPerRho[1]*length(rhos)) <= nSimPerRho[2]*length(rhos)) {
  thisN = ns[2]
  rho = rhos[ceiling((index - nSimPerRho[1]*length(rhos))/nSimPerRho[2])]
} else {
  thisN = ns[3]
  rho = rhos[ceiling((index - sum(nSimPerRho[1:2])*length(rhos))/nSimPerRho[3])]
}

system.time(out <- griddedResTestIter2Datasets(iter=index, allSeeds=seeds, 
                                               rho=rho, n1=thisN))
# })
# save(p, file="savedOutput/simStudyResults/tempFiles/profFile.RData")

# Rprof(NULL)
# profvis(prof_input = "data.Rprof")