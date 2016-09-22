library(igraph)
set.seed(11111)
example = SIMLR(X = test$in_X, c = test$n_clust, cores.ratio = 0)
ranks = SIMLR_Feature_Ranking(A = test$results$S, X = test$in_X)

context("SIMLR")
test_that("structure of output is compliant", {
	expect_equal(names(example), c("y", "S", "F", "ydata",
		"alphaK", "execution.time", "converge", "LF"))
})

test_that("correct approximation of output", {
	expect_equal(
		round(compare(test$true_labs[,1], example$y$cluster, method="nmi"), 2), 
		round(0.888298, 2))
})

context("SIMLR ranking")
test_that("structure of output is compliant", {
	expect_equal(names(ranks), c("pval", "aggR"))
})
