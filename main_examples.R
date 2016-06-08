# required external packages for SIMLR
library(Matrix)
library(parallel)

# load the igraph package to compute the NMI
library(igraph)

# load the SIMLR R package
source("./R/SIMLR.R")
source("./R/compute.multiple.kernel.R")
source("./R/network.diffusion.R")
source("./R/utils.simlr.R")
source("./R/tsne.R")
source("./R/tsne-internal.R")

# load the C file
dyn.load("./src/projsplx_R.so")

# load the datasets
load(file="./data/Test_1_mECS.RData")
load(file="./data/Test_2_Kolod.RData")
load(file="./data/Test_3_Pollen.RData")
load(file="./data/Test_4_Usoskin.RData")

# test SIMLR.R on example 1
set.seed(11111)
res_example1 = SIMLR(X=Test_1_mECS$in_X,c=Test_1_mECS$n_clust)
nmi_1 = compare(Test_1_mECS$true_labs[,1],res_example1$y$cluster,method="nmi")

# test SIMLR.R on example 2
set.seed(22222)
res_example2 = SIMLR(X=Test_2_Kolod$in_X,c=Test_2_Kolod$n_clust)
nmi_2 = compare(Test_2_Kolod$true_labs[,1],res_example2$y$cluster,method="nmi")

# test SIMLR.R on example 3
set.seed(33333)
res_example3 = SIMLR(X=Test_3_Pollen$in_X,c=Test_3_Pollen$n_clust)
nmi_3 = compare(Test_3_Pollen$true_labs[,1],res_example3$y$cluster,method="nmi")

# test SIMLR.R on example 4
set.seed(44444)
res_example4 = SIMLR(X=Test_4_Usoskin$in_X,c=Test_4_Usoskin$n_clust)
nmi_4 = compare(Test_4_Usoskin$true_labs[,1],res_example4$y$cluster,method="nmi")
