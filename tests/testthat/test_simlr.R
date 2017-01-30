library(igraph)
set.seed(11111)

normal = SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0)
if.impute = SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0, if.impute = TRUE)
normalise = SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0, normalize = TRUE)
ranks = SIMLR_Feature_Ranking(A = BuettnerFlorian$results$S, X = BuettnerFlorian$in_X)
normal_large_scale = SIMLR_Large_Scale(X = ZeiselAmit$in_X, c = ZeiselAmit$n_clust)

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

context("SIMLR large scale")
test_that("structure of output is compliant", {
    expect_equal(names(normal_large_scale), c("y", "S0", "F", "ydata",
        "alphaK", "val", "ind", "execution.time"))
})

library(scran)
ncells = 100
ngenes = 50
mu <- 2^runif(ngenes, 3, 10)
gene.counts <- matrix(rnbinom(ngenes*ncells, mu=mu, size=2), nrow=ngenes)
rownames(gene.counts) = paste0("X", seq_len(ngenes))
sce = newSCESet(countData=data.frame(gene.counts))
output = SIMLR(X = sce, c = 8, cores.ratio = 0)

context("SIMRL SCESet")
test_that("structure of output is compliant", {
    expect_equal(names(output), c("y", "S", "F", "ydata",
        "alphaK", "execution.time", "converge", "LF"))
})
