# required external Packages
library(Matrix)
library(parallel)

# load the SIMLR R package
source("SIMLR.R")
source("compute.multiple.kernel.R")
source("network.diffusion.R")
source("utils.simlr.R")

# load the c file
dyn.load("./../src/projsplx_R.so")

# load the R scripts
source("./../src/tsne.R")
source("./../src/tsne-internal.R")

# load the datasets
load(file="./../data/Test_1_mECS.RData")
load(file="./../data/Test_2_Kolod.RData")
load(file="./../data/Test_3_Pollen.RData")
load(file="./../data/Test_4_Usoskin.RData")

# test SIMLR.R on example 1
set.seed(11111)
res_example1 = SIMLR(X=Test_1_mECS$in_X,c=Test_1_mECS$n_clust)

# test SIMLR.R on example 2
set.seed(22222)
res_example2 = SIMLR(X=Test_2_Kolod$in_X,c=Test_2_Kolod$n_clust)

# test SIMLR.R on example 3
set.seed(33333)
res_example3 = SIMLR(X=Test_3_Pollen$in_X,c=Test_3_Pollen$n_clust)

# test SIMLR.R on example 4
set.seed(44444)
res_example4 = SIMLR(X=Test_4_Usoskin$in_X,c=Test_4_Usoskin$n_clust)
