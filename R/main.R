# required external Packages
library(Matrix)
library(parallel)
library(Rtsne)

# load the c files
dyn.load("./../src/projsplx_R.so")

# test compute.multiple.kernel.R
load(file="multiple.kernel.data.reduced.RData")
source("compute.multiple.kernel.R")
res1 = multiple.kernel(multiple.kernel.data.reduced)

# test network.diffusion.R
load(file="network.diffiusion.data.reduced.RData")
source("network.diffusion.R")
res2 = network.diffusion(network.diffiusion.data.reduced,K=10)

# test SIMLR.R on reduced data
load(file="my_data_reduced.Rdata")
source("SIMLR.R")
source("compute.multiple.kernel.R")
source("network.diffusion.R")
source("utils.simlr.R")
res3 = SIMLR(my_data_reduced,c=10)

# test SIMLR.R on a bigger data input
load(file="my_data.Rdata")
source("SIMLR.R")
source("compute.multiple.kernel.R")
source("network.diffusion.R")
source("utils.simlr.R")
res4 = SIMLR(my_data[1:70,1:1000],c=10)
