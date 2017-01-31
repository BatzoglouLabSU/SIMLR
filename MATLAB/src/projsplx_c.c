#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    int numDims, m, n, k, d, j,npos,ft;
    const mwSize *dims;
    double *y, *s, *x, *vs;
    double sumResult = -1, tmpValue, tmax, f,lambda_m;
    
    
    dims = mxGetDimensions(prhs[0]);
    numDims = mxGetNumberOfDimensions(prhs[0]);

    m = dims[0];
    n=dims[1];
    
    y  = mxGetPr(prhs[0]);
    
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
    x = mxGetPr(plhs[0]);
    s = (double*) calloc (m,sizeof(double));
    vs = (double*) calloc (m,sizeof(double));
    for (k=0;k<n;k++){
        /* s = sort(y,'ascend'); */
        
        double means = 0;
        double mins = 100000;
        for(j = 0; j < m; j++ ){
            s[j] = y[j+k*m];
            means += s[j];
            mins = (mins > s[j])? s[j]:mins;
        }

        for(j = 0; j < m; j++ ){
            s[j] -= (means-1)/m;
        }
        ft=1;
        if(mins<0){
            f = 1;
            lambda_m=0;
            while(fabs(f)>1e-10){
                npos = 0;
                   f = 0;
                for(j = 0; j < m; j++ ){
                    vs[j] = s[j]-lambda_m;
                 
                    if (vs[j]>0){
                        npos+=1;
                        f+=vs[j];
                    }
                }
                lambda_m += (f-1)/npos;
                if(ft>100){
                    for(j = 0; j <= m-1; j++){
                        x[j+k*m] = (vs[j] > 0)? vs[j]:0;
                    }
                    break;
                }
                ft+=1;
            }
            for(j = 0; j <= m-1; j++){
                x[j+k*m] = (vs[j] > 0)? vs[j]:0;;
            }
            
        }
        else{
            for(j = 0; j <= m-1; j++){
                x[j+k*m] = s[j];
            }
            
        }
    }
    
}

