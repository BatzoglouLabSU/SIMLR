# required external packages for CIMLR
library(Matrix)
library(parallel)

# load the CIMLR R package
source("./R/CIMLR.R")
source("./R/compute.multiple.kernel.cimlr.R")
source("./R/network.diffusion.R")
source("./R/utils.simlr.R")
source("./R/tsne.R")

# load the C file

# NOTE 1: we make use of an external C program during the computations of CIMLR.
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
load(file="./data/Test_6_gliomas_multi_omic_data.RData")

# test CIMLR.R on a dataset of lower grade gliomas
set.seed(75584)
cat("Performing analysis for Test_6_gliomas_multi_omic_data","\n")
res_llg = CIMLR(X=Test_6_gliomas_multi_omic_data$in_X,c=Test_6_gliomas_multi_omic_data$n_clust)
