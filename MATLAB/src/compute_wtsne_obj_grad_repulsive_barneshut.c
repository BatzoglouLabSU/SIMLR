/*
 * Compute the repulsive part of objective and gradient of  t-SNE
 *
 * All rights reserved.
 *
 */

#include "mex.h"
#include <math.h>
#include "barnes_hut.h"

void getRepulsiveObjGradi(int id, double *pos, QuadTree* tree, double farness_factor, int nArgOut, double *pqsum, double *grad) {
    double dist2, dist, diff, decoef, q;
    int i,d;
    double tmp;
    
    if (tree==NULL) return;
    if (tree->node!=NULL && tree->node->id==id) return;
    
    dist2 = 0.0;
    for (d=0;d<2;d++) {
        diff = pos[d]-tree->position[d];
        dist2 += diff*diff;
    }
    
    tmp = farness_factor*getTreeWidth(tree);
    if (tree->childCount>0 && dist2<tmp*tmp) {
        for(i=0;i<tree->childrenLength;i++) {
            if (tree->children[i]!=NULL) {
                getRepulsiveObjGradi(id, pos, tree->children[i], farness_factor, nArgOut, pqsum, grad);
            }
        }
    } else {
        q = 1.0 / (1+dist2);
        (*pqsum) += q * tree->weight;
        
        if (nArgOut>1) {
            tmp = tree->weight * q * q;
            for (d=0;d<2;d++) {
                grad[d] += (tree->position[d] - pos[d]) * tmp;
            }
        }
    }
}

void getRepulsiveObjGrad(double *Y, double *weights, int n, double eps, double farness_factor, int nArgOut, double *pobj, double *grad) {
    double pos[2], gradi[2], d2, diff;
    double *coef;
    mwIndex i, j, d;
    QuadTree *tree;
    double qsumi, qsum;
    
    tree = buildQuadTree(Y,weights,n);

    qsum = 0.0;
    for (i=0;i<n;i++) {
        for(d=0;d<2;d++) {
            gradi[d] = 0.0;
            pos[d] = Y[d*n+i];
        }
        qsumi = 0.0;
        getRepulsiveObjGradi(i, pos, tree, farness_factor, nArgOut, &qsumi, gradi);
        qsum += qsumi*weights[i];
        if (nArgOut>1) {
            for(d=0;d<2;d++)
                grad[d*n+i] = 4 * gradi[d] * weights[i];
        }
    }
    
    if (nArgOut>1) {
        for (i=0;i<n;i++) {
            for(d=0;d<2;d++)
                grad[d*n+i] /= qsum;
        }
    }

    (*pobj) = log(qsum+eps);
    
    destroyQuadTree(tree);
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *Y, eps, *pobj, *grad, farness_factor, *weights;
    int n, nArgOut;
    
    Y = mxGetPr(prhs[0]); /* note that Y is a vector */
    weights = mxGetPr(prhs[1]);
    farness_factor = mxGetScalar(prhs[2]);
    nArgOut = (int)mxGetScalar(prhs[3]);
    
    n = mxGetM(prhs[0]);
    eps = mxGetEps();
    
    plhs[0] = mxCreateDoubleScalar(0.0);
    pobj = mxGetPr(plhs[0]);
    if (nArgOut>1) {
        plhs[1] = mxCreateDoubleMatrix(n,2,mxREAL);
        grad = mxGetPr(plhs[1]);
    }
    else
        grad = NULL;

    getRepulsiveObjGrad(Y, weights, n, eps, farness_factor, nArgOut, pobj, grad);
}