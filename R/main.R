# load the data
load(file="multiple.kernel.data.reduced.RData")

# test compute.multiple.kernel.R
source("compute.multiple.kernel.R")
res = multiple.kernel(multiple.kernel.data.reduced)
