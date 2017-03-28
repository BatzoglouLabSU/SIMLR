# Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding modified for SIMLR large scale 
Rtsne <- function( I, J, V, dims = 2, perplexity = 30, theta = 0.5, verbose = FALSE, max_iter = 300, Y_init = matrix(), init=FALSE, ... ) {
  return(.Rtsne_cpp(I,J,V,dims,perplexity,theta,verbose,max_iter,Y_init,init))
}
