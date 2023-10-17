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
nSimulations=10
set.seed(123)
seeds = sample(1:100000, nSimulations)
thisSeed=seeds[index]
system.time(out <- griddedResTestIter(iter=index, seed=thisSeed))
# })
# save(p, file="savedOutput/simStudyResults/tempFiles/profFile.RData")

# Rprof(NULL)
# profvis(prof_input = "data.Rprof")