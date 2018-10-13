Sys.setenv("R_TESTS" = "")

library("testthat")
library("SIMLR")

test_check("SIMLR")
