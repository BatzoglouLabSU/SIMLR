library(igraph)
set.seed(11111)

ranks = SIMLR_Feature_Ranking(A = test$results$S, X = test$in_X)
normal = SIMLR(X = test$in_X, c = test$n_clust, cores.ratio = 0)
if.impute = SIMLR(X = test$in_X, c = test$n_clust, cores.ratio = 0, if.impute = TRUE)
normalise = SIMLR(X = test$in_X, c = test$n_clust, cores.ratio = 0, normalize = TRUE)

context("SIMLR")
test_that("structure of output is compliant", {
	expect_equal(names(normal), c("y", "S", "F", "ydata",
		"alphaK", "execution.time", "converge", "LF"))
	expect_equal(names(if.impute), c("y", "S", "F", "ydata",
		"alphaK", "execution.time", "converge", "LF"))
	expect_equal(names(normalise), c("y", "S", "F", "ydata",
		"alphaK", "execution.time", "converge", "LF"))
})

context("SIMLR ranking")
test_that("structure of output is compliant", {
	expect_equal(names(ranks), c("pval", "aggR"))
})
