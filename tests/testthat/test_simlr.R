normal = SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0)
if.impute = SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0, if.impute = TRUE)
normalise = SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0, normalize = TRUE)
ranks = SIMLR_Feature_Ranking(A = BuettnerFlorian$results$S, X = BuettnerFlorian$in_X)

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