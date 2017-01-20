# required external packages for SIMLR
library(Matrix)
library(parallel)
library(rsvd)
library(largeVis)
library(RSpectra)

# load the igraph package to compute the NMI
library(igraph)

# load the palettes for the plots
library(grDevices)

# load the SIMLR R package
source("./R/SIMLR_Large_Scale.R")
source("./R/utils.simlr.large.scale.R")
source("./R/utils.simlr.R")

# load the C file
dyn.load("./R/projsplx_R.so")

# load the datasets
load(file="./data/Zelsel.RData")

# test SIMLR.R on the large scale dataset
set.seed(12345)
res_large_scale = SIMLR_Large_Scale(X=Zelsel$in_X,c=Zelsel$n_clust)
nmi_large_scale = compare(Zelsel$true_labs[,1],res_large_scale$y$cluster,method="nmi")
