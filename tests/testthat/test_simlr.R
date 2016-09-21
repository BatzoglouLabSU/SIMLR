library(igraph)
set.seed(11111)
example = SIMLR(X=test$in_X, c=test$n_clust)


context("SIMLR")
test_that("structure of output", {
	expect_equal(names(example), c("y", "S", "F", "ydata",
		"alphaK", "execution.time", "converge", "LF"))
})

test_that("test2", {
	expect_equal(
		round(compare(test$true_labs[,1], example$y$cluster, method="nmi"), 2), 
		round(0.888298, 2))
})
