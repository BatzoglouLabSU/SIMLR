# test compute.multiple.kernel.R
load(file="multiple.kernel.data.reduced.RData")
source("compute.multiple.kernel.R")
res1 = multiple.kernel(multiple.kernel.data.reduced)

# test network.diffusion.R
load(file="multiple.kernel.data.reduced.RData")
source("network.diffusion.R")
res2 = network.diffusion(multiple.kernel.data.reduced,K=10)
