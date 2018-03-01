# required external packages for CIMLR
library(Matrix)
library(parallel)

# load the CIMLR R package
source("./R/SIMLR.R")
source("./R/compute.multiple.kernel.R")
source("./R/network.diffusion.R")
source("./R/utils.simlr.R")
source("./R/tsne.R")

# load the C file

# NOTE 1: we make use of an external C program during the computations of SIMLR.
# The code is located in the R directory in the file projsplx_R.c. In order to 
# use SIMLR one needs to compite the program. To do so, one needs to run on the 
# shell the command R CMD SHLIB -c projsplx_R.c. 
# The precompiled projsplx_R.so is already provided for MAC OS X only. 
# If one wants to use SIMLR on other operative systems, the file projsplx_R.so 
# needs to be deleted, and re-compiled. 

# NOTE 2: for Windows, the command dyn.load("./R/projsplx_R.so") needs to be 
# substituted with the command dyn.load("./R/projsplx_R.dll"). 

dyn.load("./R/projsplx_R.so")

# load the datasets
load(file="./data/Test_1_mECS.RData")

# test SIMLR.R on example 1
set.seed(11111)
cat("Performing analysis for Test_1_mECS","\n")
res_example1 = SIMLR(X=Test_1_mECS$in_X,c=Test_1_mECS$n_clust)
nmi_1 = compare(Test_1_mECS$true_labs[,1],res_example1$y$cluster,method="nmi")
