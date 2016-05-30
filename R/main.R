# test compute.multiple.kernel.R
load(file="multiple.kernel.data.reduced.RData")
source("compute.multiple.kernel.R")
res1 = multiple.kernel(multiple.kernel.data.reduced)

# test network.diffusion.R
load(file="network.diffiusion.data.reduced.RData")
source("network.diffusion.R")
res2 = network.diffusion(network.diffiusion.data.reduced,K=10)

# test SIMLR.R
load(file="my_data_reduced.Rdata")
source("SIMLR.R")
source("compute.multiple.kernel.R")
source("network.diffusion.R")
source("utils.simlr.R")
res3 = SIMLR(my_data_reduced,c=90)
