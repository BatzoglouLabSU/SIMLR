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
#include "annoylib.h"
#include "kissrandom.h"
#include <time.h>
#include <map>

void make_top_neighbors_annoy(double *X, int NK, int NN,int NF,double *ind,double *val){
    clock_t begin = clock();
    
    AnnoyIndex<int, double, Euclidean, Kiss64Random> t = AnnoyIndex<int, double, Euclidean, Kiss64Random>(NF);
    for(int i = 0; i<NN; i++){
        double *vec = (double *) malloc( NF * sizeof(double) );
        for(int z=0; z<NF; ++z){
			vec[z] = X[ i + NN*z];
		}
        t.add_item(i, vec);
    }
    clock_t end = clock();
    printf("Elapsed time in copying data is %f seconds\n", (double)(end - begin)/CLOCKS_PER_SEC);
    begin = clock();
    t.build(100);
    end = clock();
    printf("Elapsed time in building trees is %f seconds\n", (double)(end - begin)/CLOCKS_PER_SEC);
    
    printf("initializatioin works\n");
    
    std::vector<int> closest;
    begin = clock();
    for (int i = 0; i< NN; i++){
       
        //std::vector<double> distances;
        t.get_nns_by_item(i, NK, -1, &closest, nullptr);
        //printf("get nn search for item %d works\n", i);
            
        for (int j = 0; j < NK; j++){
            ind[i + j*NN] = double(closest[j]);
            val[i + j*NN] = t.get_distance(i,closest[j] )   ;
        }
        closest.clear(); vector<int>().swap(closest);
        //distances.clear(); 
    }
   end = clock();
   printf("Elapsed time in assign values is %f seconds\n", (double)(end - begin)/CLOCKS_PER_SEC);
    ///
}

/// usage: eigenvectors = mex_top_eig(val, ind, KK);
/// input: val of size NxK, the value of transition matrix 
///        ind of size NxK, the index of the values (Note ind starts with 0);
///        KK , the number of eigenvectors    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    double *val, *ind, *X;
    int KK, NN,NF;
    X = mxGetPr(prhs[0]);
    KK = (int) mxGetScalar(prhs[1]); //number of eigenvalues
    printf("KK is %d \n", KK);
    NN = mxGetM(prhs[0]);
    NF = mxGetN(prhs[0]);
    printf("NN is %d \n", NN);
    printf("NF is %d \n", NF);
    plhs[0] = mxCreateDoubleMatrix(NN,KK,mxREAL);
    ind = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(NN,KK,mxREAL);
    val = mxGetPr(plhs[1]);
    make_top_neighbors_annoy(X, KK, NN, NF, ind,val);
}
