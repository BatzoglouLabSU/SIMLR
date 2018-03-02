# required external packages for CIMLR
library(Matrix)
library(parallel)

# load the CIMLR R package
source("./R/CIMLR_Estimate_Number_of_Clusters.R")
source("./R/compute.multiple.kernel.R")
source("./R/network.diffusion.R")

# load the datasets
load(file="./data/Test_6_gliomas_multi_omic_data.RData")

# in this example we test the estimation of the number of cluster on Test_6_gliomas_multi_omic_data 
# where the number of clusters is supposed to be 3 
set.seed(86677)
cat("Performing analysis for Test_6_gliomas_multi_omic_data","\n")
NUMC = 2:15
res_example = CIMLR_Estimate_Number_of_Clusters(Test_6_gliomas_multi_omic_data$in_X,NUMC=NUMC,cores.ratio=1)

# the best number of clusters is the one with lower values for the two heuristics
print(paste0("Best number of clusters, K1 heuristic: ",NUMC[which.min(res_example$K1)],", K2 heuristic: ",NUMC[which.min(res_example$K2)]))
print(res_example)
