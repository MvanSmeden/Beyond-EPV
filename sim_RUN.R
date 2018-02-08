# This file is part of combination of R-files that contain code to perform a simulation study
# The code was used for a simulation study for which the results are presented in a currently unpublished manuscript

###################################################
#### Execute simulation runs and store results ####
#### Last modified: Feb 8 2018                 ####
#### Author: M van Smeden                      ####
###################################################

#setwd("path")  # note: set working directory or start a new project in R-studio

# BEGIN CODE ---------------------------------------

source("sim_dependencies.R")
source("sim_gendata.R")
source("sim_LRMs.R")
source("sim_stats.R")
source("sim_simroutine.R")

simlists <- readRDS(file="prelim/mainsimlists.rds")  

iter <- 1 # note: one could build a loop to run multiple iterations

for(i in 1:10){
  condno <- i
  condID <- names(simlists)[condno]
  simlist <-simlists[[condID]]$development
  simlist.validation <- simlists[[condID]]$validation
  seeds <- sample(1:1e6,3)   # note: in our simulation the seeds were pre-defined
  single.run(simlist,simlist.validation,seeds,condno,condID,iter)
}

## --------------------------------------------------------

