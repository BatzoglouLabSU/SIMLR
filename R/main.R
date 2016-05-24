# load the data
load(file="my_data_reduced.RData")

# test compute.multiple.kernel.R
source("compute.multiple.kernel.R")
res = multiple.kernel(t(my_data_reduced))
