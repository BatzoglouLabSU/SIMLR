# Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding 
# 
# Wrapper for the C++ implementation of Barnes-Hut t-Distributed Stochastic Neighbor Embedding. t-SNE is a method for constructing a low dimensional embedding of high-dimensional data, distances or similarities. Exact t-SNE can be computed by setting theta=0.0. 
# 
# Given a distance matrix \eqn{D} between input objects (which by default, is the euclidean distances between two objects), we calculate a similarity score in the original space p_ij. \deqn{ p_{j | i} = \frac{\exp(-\|D_{ij}\|^2 / 2 \sigma_i^2)}{\sum_{k \neq i} \exp(-\|D_{ij}\|^2 / 2 \sigma_i^2)} } which is then symmetrized using: \deqn{ p_{i j}=\frac{p_{j|i} + p_{i|j}}{2n}}. The \eqn{\sigma} for each object is chosen in such a way that the perplexity of p_{j|i} has a value that is close to the user defined perplexity. This value effectively controls how many nearest neighbours are taken into account when constructing the embedding in the low-dimensional space.
# For the low-dimensional space we use the Cauchy distribution (t-distribution with one degree of freedom) as the distribution of the distances to neighbouring objects:
# \deqn{ q_{i j} = \frac{(1+ \| y_i-y_j\|^2)^{-1}}{\sum_{k \neq l} 1+ \| y_k-y_l\|^2)^{-1}}}. 
# By changing the location of the objects y in the embedding to minimize the Kullback-Leibler divergence between these two distributions \eqn{ q_{i j}} and \eqn{ p_{i j}}, we create a map that focusses on small-scale structure, due to the assymetry of the KL-divergence. The t-distribution is chosen to avoid the crowding problem: in the original high dimensional space, there are potentially many equidistant objects with moderate distance from a particular object, more than can be accounted for in the low dimensional representation. The t-distribution makes sure that these objects are more spread out in the new representation.
# 
# For larger datasets, a problem with the a simple gradient descent to minimize the Kullback-Leibler divergence is the computational complexity of each gradient step (which is \eqn{O(n^2)}). The Barnes-Hut implementation of the algorithm attempts to mitigate this problem using two tricks: (1) approximating small similarities by 0 in the \eqn{p_{ij}} distribution, where the non-zero entries are computed by finding 3*perplexity nearest neighbours using an efficient tree search. (2) Using the Barnes-Hut algorithm in the computation of the gradient which approximates large distance similarities using a quadtree. This approximation is controlled by the \code{theta} parameter, with smaller values leading to more exact approximations. When \code{theta=0.0}, the implementation uses a standard t-SNE implementation. The Barnes-Hut approximation leads to a \eqn{O(n log(n))} computational complexity for each iteration.
# 
# During the minimization of the KL-divergence, the implementation uses a trick known as early exaggeration, which multiplies the \eqn{p_{ij}}'s by 12 during the first 250 iterations. This leads to tighter clustering and more distance between clusters of objects. This early exaggeration is not used when the user gives an initialization of the objects in the embedding by setting \code{Y_init}. During the early exaggeration phase, a momentum term of 0.5 is used while this is changed to 0.8 after the first 250 iterations.
# 
# After checking the correctness of the input, the \code{Rtsne} function (optionally) does an initial reduction of the feature space using \code{\link{prcomp}}, before calling the C++ TSNE implementation. Since R's random number generator is used, use \code{\link{set.seed}} before the function call to get reproducible results.
# 
SIMLR.Rtsne <- function( I, J, V, dims = 2, perplexity = 30, theta = 0.5, verbose = FALSE, max_iter = 300, Y_init = NULL, ...) {
  
  if (is.null(Y_init)) {
    init <- FALSE
    Y_init <- matrix()
  } else {
    init <- TRUE
  }
  
  return(Rtsne_cpp(I,J,V,dims,perplexity,theta,verbose,max_iter,Y_init,init))
  
}
