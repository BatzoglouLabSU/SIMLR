#include <Rcpp.h>
#include "tsne.h"
using namespace Rcpp;

// Function that runs the Barnes-Hut implementation of t-SNE
// [[Rcpp::export]]
Rcpp::List Rtsne_cpp(NumericVector I, NumericVector J, NumericVector V, int no_dims_in, double perplexity_in, double theta_in, bool verbose, int max_iter, NumericMatrix Y_in, bool init) {

  int origN, N, D, no_dims = no_dims_in;

  int  *dataI;
  int  *dataJ;
  double  *dataV;
  TSNE* tsne = new TSNE();
  double perplexity = perplexity_in;
  double theta = theta_in;

  N = I.size()-1;
  origN = J.size();
  D = 1; //I.size();
    
  dataI = (int*) calloc(D * N, sizeof(int));
  dataJ = (int*) calloc(D * origN, sizeof(int));
  dataV = (double*) calloc(D * origN, sizeof(double));
  if(dataI == NULL || dataJ == NULL || dataV == NULL) { Rcpp::stop("Memory allocation failed!\n"); }
  for (int i = 0; i < N+1; i++){
      dataI[i*D] = I(i)-1;
  }
  for (int i = 0; i < origN; i++){
      dataJ[i*D] = J(i)-1;
      dataV[i*D] = V(i);
  }
    
  // Make dummy landmarks
  if (verbose) Rprintf("Read the %i x %i data matrix successfully!\n", N, D);
  int* landmarks = (int*) malloc(N * sizeof(int));
  if(landmarks == NULL) { Rcpp::stop("Memory allocation failed!\n"); }
  for(int n = 0; n < N; n++) landmarks[n] = n;

	double* Y = (double*) malloc(N * no_dims * sizeof(double));
	double* costs = (double*) calloc(N, sizeof(double));
	double* itercosts = (double*) calloc((int)(ceil(max_iter/50.0)), sizeof(double));
  if(Y == NULL || costs == NULL) { Rcpp::stop("Memory allocation failed!\n"); }
    
  // Initialize solution (randomly)
  if (init) {
    for (int i = 0; i < N; i++){
      for (int j = 0; j < no_dims; j++){
        Y[i*no_dims+j] = Y_in(i,j);
      }
    }
    if (verbose) Rprintf("Using user supplied starting positions\n");
  }
    
  // Run tsne
  tsne->run(dataI, dataJ, dataV, N, Y, no_dims, theta, verbose, max_iter, costs, itercosts, init);

  // Save the results
  Rcpp::NumericMatrix Yr(N, no_dims);
  for (int i = 0; i < N; i++){
      for (int j = 0; j < no_dims; j++){
          Yr(i,j) = Y[i*no_dims+j];
      }
  }
    
  Rcpp::NumericVector costsr(N);
  for (int i = 0; i < N; i++){
    costsr(i) = costs[i];
  }
  Rcpp::NumericVector itercostsr((int)(ceil(max_iter/50.0)));
  for (int i = 0; i < (int)(ceil(max_iter/50.0)); i++) {
    itercostsr(i) = itercosts[i];
  }
  
  //free(dataI); dataI = NULL;
  dataI = NULL;
  //free(dataJ); dataJ = NULL;
  dataJ = NULL;
  //free(dataV); dataV = NULL;
  dataV = NULL;
  free(Y); Y = NULL;
  free(costs); costs = NULL;
  free(landmarks); landmarks = NULL;
  delete(tsne);
  
  Rcpp::List output = Rcpp::List::create(Rcpp::_["theta"]=theta, Rcpp::_["perplexity"]=perplexity, Rcpp::_["N"]=N,Rcpp::_["origD"]=D,Rcpp::_["Y"]=Yr, Rcpp::_["costs"]=costsr, Rcpp::_["itercosts"]=itercostsr);
  return output;
}
