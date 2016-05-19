# load the data
load(file="my_data.RData")

# run network diffusion
source("network.diffusion.R")
res_network_diffusion = network.diffusion(my_data,100)
