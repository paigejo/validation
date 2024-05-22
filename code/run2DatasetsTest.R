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
  iter = ((index-1) %% nSimPerRho[1]) + 1
  thisSeeds = seeds[1:nSimPerRho[1]]
} else if((index - nSimPerRho[1]*length(rhos)) <= nSimPerRho[2]*length(rhos)) {
  thisN = ns[2]
  rho = rhos[ceiling((index - nSimPerRho[1]*length(rhos))/nSimPerRho[2])]
  iter = ((index-1 - nSimPerRho[1]*length(rhos)) %% nSimPerRho[2]) + 1
  thisSeeds = seeds[(nSimPerRho[1]+1):(nSimPerRho[1]+nSimPerRho[2])]
} else {
  thisN = ns[3]
  rho = rhos[ceiling((index - sum(nSimPerRho[1:2])*length(rhos))/nSimPerRho[3])]
  iter = ((index-1 - sum(nSimPerRho[1:2])*length(rhos)) %% nSimPerRho[3]) + 1
  thisSeeds = seeds[(nSimPerRho[2]+1):(nSimPerRho[2]+nSimPerRho[3])]
}
# indices 1-5000: n1=50
# indices 5001-5500: n1=500
# indices 5501-5750: n1=1000
system.time(out <- griddedResTestIter2Datasets(iter=iter, allSeeds=thisSeeds, 
                                               rho=rho, n1=thisN))
# for n1=500: 79.679 seconds per iter
# For n1=1000: 685.924 seconds per iter
# })
# save(p, file="savedOutput/simStudyResults/tempFiles/profFile.RData")

# Rprof(NULL)
# profvis(prof_input = "data.Rprof")