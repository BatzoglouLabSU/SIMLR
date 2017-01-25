# required external packages for SIMLR
library(Rcpp)
library(Matrix)
library(rsvd)
library(RANN)
library(rARPACK)

# load the igraph package to compute the NMI
library(igraph)

# load the palettes for the plots
library(grDevices)

# load the SIMLR R package
source("./R/SIMLR_Large_Scale.R")
source("./R/utils.simlr.large.scale.R")
source("./R/utils.simlr.R")
source("./R/SIMLR.Rtsne.R")

# load the C files
dyn.load("./R/projsplx_R.so")
sourceCpp("./src/Rtsne.cpp")

# load the datasets
load(file="./data/Zelsel.RData")

# test SIMLR.R on the large scale dataset
set.seed(12345)
res_large_scale = SIMLR_Large_Scale(X=Zelsel$in_X,c=Zelsel$n_clust)
nmi_large_scale = compare(Zelsel$true_labs[,1],res_large_scale$y$cluster,method="nmi")

# make the scatterd plots
plot(res_large_scale$ydata,col=c(topo.colors(Zelsel$n_clust))[Zelsel$true_labs[,1]],xlab="SIMLR component 1", ylab="SIMLR component 2",pch=20,main="SIMILR 2D visualization for Zelsel")
