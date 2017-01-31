//
//  mex_top_eig.cpp
//  
//
//  Created by Bo_Royce on 8/17/16.
//
//
#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/Core>
#include <SymEigsSolver.h>
#include <MatOp/SparseSymMatProd.h>
#include <time.h>
using namespace Eigen;
using namespace Spectra;
void make_top_eigenvectors(double *val,double *ind, int KK, int NN, int NK, double *eigenvectors,double *eigenvalues){
    clock_t begin = clock();
    Eigen::SparseMatrix<double> mat((const int) NN,(const int) NN);         // default is column major
    mat.reserve(Eigen::VectorXi::Constant((const int) NN, (const int) NK));
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve((const int) NN*NK);
    for(int i=0; i<NN; i++){
        for (int j = 0; j< NK; j++){
            tripletList.push_back(T((int) ind[i+NN*j]-1,i,val[i+j*NN]));
        }
    }
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    mat += Eigen::SparseMatrix<double>(mat.transpose());
    clock_t end = clock();
    //printf("Elapsed time in initialization is %f seconds\n", (double)(end - begin)/CLOCKS_PER_SEC);
    
    SparseSymMatProd<double> op(mat);
    begin = clock();
    // Construct eigen solver object, requesting the largest KK eigenvalues
    SymEigsSolver< double, LARGEST_ALGE, SparseSymMatProd<double> > eigs(&op, KK, 2*KK);
    
    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute();
    
    // Retrieve results
    
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evectors;
    if(eigs.info() == SUCCESSFUL){
        evalues = eigs.eigenvalues();
        evectors = eigs.eigenvectors();
    }
    //std::cout << "Eigenvalues found:\n" << evalues << std::endl;
    end = clock();
    printf("Elapsed time in eigen-decomposition is %f seconds\n", (double)(end - begin)/CLOCKS_PER_SEC);
    ///
    begin = clock();
    for (int j = 0; j< KK; j++){
        eigenvalues[j] = evalues[j];
        for (int i = 0; i < NN; i++){
            eigenvectors[i+NN*j] = evectors.col(j)[i];
        }
    }
   end = clock();
    printf("Elapsed time in copying eigenvectors is %f seconds\n", (double)(end - begin)/CLOCKS_PER_SEC);
    ///
}

/// usage: eigenvectors = mex_top_eig(val, ind, KK);
/// input: val of size NxK, the value of transition matrix 
///        ind of size NxK, the index of the values (Note ind starts with 0);
///        KK , the number of eigenvectors    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    double *val, *ind, *eigenvectors,*eigenvalues;
    int KK, NN, NK;
    val = mxGetPr(prhs[0]);
    ind = mxGetPr(prhs[1]); //
    KK = (int) mxGetScalar(prhs[2]); //number of eigenvalues
    NN = mxGetM(prhs[0]);
    NK = mxGetN(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(NN,KK,mxREAL);
    eigenvectors = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(KK,1,mxREAL);
    eigenvalues = mxGetPr(plhs[1]);
    make_top_eigenvectors(val,ind, KK,  NN, NK, eigenvectors,eigenvalues);
}
