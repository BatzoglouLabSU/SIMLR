# required external packages for SIMLR
library(Matrix)
library(parallel)

# load the SIMLR R package
source("./R/SIMLR_Estimate_Number_of_Clusters.R")
source("./R/compute.multiple.kernel.R")
source("./R/network.diffusion.R")

# load the datasets
load(file="./data/Test_2_Kolod.RData")

# in this example we test the estimation of the number of cluster on Test_1_mECS 
# where the number of clusters is supposed to be 3 
set.seed(32411)
cat("Performing analysis for Test_2_Kolod","\n")
NUMC = 2:5
res_example = SIMLR_Estimate_Number_of_Clusters(Test_2_Kolod$in_X,NUMC=NUMC,cores.ratio=1)

# the best number of clusters is the one with lower values for the two heuristics
print(paste0("Best number of clusters, K1 heuristic: ",NUMC[which.min(res_example$K1)],", K2 heuristic: ",NUMC[which.min(res_example$K2)],"."))
print(res_example)
