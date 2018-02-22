# required external packages for SIMLR
library(Matrix)
library(parallel)

# load the igraph package to compute the NMI
library(igraph)

# load the palettes for the plots
library(grDevices)

# load the SIMLR R package
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
load(file="./data/Test_2_Kolod.RData")
load(file="./data/Test_3_Pollen.RData")
load(file="./data/Test_4_Usoskin.RData")

# test SIMLR.R on example 1
set.seed(11111)
cat("Performing analysis for Test_1_mECS","\n")
res_example1 = SIMLR(X=Test_1_mECS$in_X,c=Test_1_mECS$n_clust)
nmi_1 = compare(Test_1_mECS$true_labs[,1],res_example1$y$cluster,method="nmi")

# test SIMLR.R on example 2
set.seed(22222)
cat("Performing analysis for Test_2_Kolod","\n")
res_example2 = SIMLR(X=Test_2_Kolod$in_X,c=Test_2_Kolod$n_clust)
nmi_2 = compare(Test_2_Kolod$true_labs[,1],res_example2$y$cluster,method="nmi")

# test SIMLR.R on example 3
set.seed(33333)
cat("Performing analysis for Test_3_Pollen","\n")
res_example3 = SIMLR(X=Test_3_Pollen$in_X,c=Test_3_Pollen$n_clust)
nmi_3 = compare(Test_3_Pollen$true_labs[,1],res_example3$y$cluster,method="nmi")

# test SIMLR.R on example 4
set.seed(44444)
cat("Performing analysis for Test_4_Usoskin","\n")
res_example4 = SIMLR(X=Test_4_Usoskin$in_X,c=Test_4_Usoskin$n_clust)
nmi_4 = compare(Test_4_Usoskin$true_labs[,1],res_example4$y$cluster,method="nmi")

# make the scatterd plots
plot(res_example1$ydata,col=c(topo.colors(Test_1_mECS$n_clust))[Test_1_mECS$true_labs[,1]],xlab="SIMLR component 1", ylab="SIMLR component 2",pch=20,main="SIMILR 2D visualization for Test_1_mECS")

plot(res_example2$ydata,col=c(topo.colors(Test_2_Kolod$n_clust))[Test_2_Kolod$true_labs[,1]],xlab="SIMLR component 1", ylab="SIMLR component 2",pch=20,main="SIMILR 2D visualization for Test_2_Kolod")

plot(res_example3$ydata,col=c(topo.colors(Test_3_Pollen$n_clust))[Test_3_Pollen$true_labs[,1]],xlab="SIMLR component 1", ylab="SIMLR component 2",pch=20,main="SIMILR 2D visualization for Test_3_Pollen")

plot(res_example4$ydata,col=c(topo.colors(Test_4_Usoskin$n_clust))[Test_4_Usoskin$true_labs[,1]],xlab="SIMLR component 1", ylab="SIMLR component 2",pch=20,main="SIMILR 2D visualization for Test_4_Usoskin")
