# required external packages for SIMLR large scale
library(Rcpp)
library(Matrix)
library(pracma)
library(RcppAnnoy)
library(RSpectra)

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
load(file="./data/Marcos.RData")

# test SIMLR.R on the large scale dataset of Zelsel
set.seed(11111)
cat("Performing analysis for Zelsel","\n")
res_large_scale_1 = SIMLR_Large_Scale(X=Zelsel$in_X,c=Zelsel$n_clust,k=30,kk=200)
nmi_large_scale_1 = compare(Zelsel$true_labs[,1],res_large_scale_1$y$cluster,method="nmi")

# test SIMLR.R on the large scale dataset of Marcos
set.seed(22222)
cat("Performing analysis for Marcos","\n")
res_large_scale_2 = SIMLR_Large_Scale(X=Marcos$in_X,c=Marcos$n_clust,k=30,kk=500)
nmi_large_scale_2 = compare(Marcos$true_labs[,1],res_large_scale_2$y$cluster,method="nmi")

# make the scatterd plots
plot(res_large_scale_1$ydata,col=c(topo.colors(Zelsel$n_clust))[Zelsel$true_labs[,1]],xlab="SIMLR component 1", ylab="SIMLR component 2",pch=20,main="SIMILR 2D visualization for Zelsel")

plot(res_large_scale_2$ydata,col=c(topo.colors(Marcos$n_clust))[Marcos$true_labs[,1]],xlab="SIMLR component 1", ylab="SIMLR component 2",pch=20,main="SIMILR 2D visualization for Marcos")
