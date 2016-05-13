#include "mex.h"
#include <math.h>

/*
 * This function computes
 * A = 0;
 * for i=1:length(beta), A = A + K(:,:,i)*beta(i); end
 * for symmetric and not-symmetric K
 */    
void mexFunction(int , mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	if (nrhs<2 || nrhs>3)
		mexErrMsgTxt("Wrong number of input arguments");

	if (mxGetNumberOfDimensions(prhs[1])!=2)
		mexErrMsgTxt("Second input argument must be 2D");

	unsigned int nDims = mxGetNumberOfDimensions(prhs[0]);
	const mwSize *dims = mxGetDimensions(prhs[0]);

	if (nDims<2 || nDims>3)
		mexErrMsgTxt("First input arg must be 2D or 3D");

	double *K = mxGetPr(prhs[0]);
	double *beta = mxGetPr(prhs[1]);

	unsigned int n1 = dims[0];
	unsigned int n2 = dims[1];
	unsigned int m = (nDims==3)?dims[2]:1;
	unsigned int nn = n1*n2;

	if (mxGetM(prhs[1]) != m)
		mexErrMsgTxt("Input dimensions mismatch");

	bool symmetric = false;
	if (n1==n2)
		if (nrhs>2)
			symmetric = (bool)(*(mxGetPr(prhs[2]))>0);

	// output arguments
	plhs[0]=mxCreateDoubleMatrix(n1,n2,mxREAL); 
	double *Kbeta = mxGetPr(plhs[0]); 

	if (symmetric)
	{
		for ( unsigned k=0,knn=0 ; k<m ; k++,knn+=nn )
		{
			if ( beta[k] <= 1e-8)
				continue;

			for ( unsigned int i=0,in = 0 ; i<n1 ; i++,in+=n1 )
				for ( unsigned int j=i ; j<n1 ; j++ )
					Kbeta[in+j] += beta[k] * K[knn+in+j];
		}

		/* symmetrize the output */
		for ( unsigned int i=0,in=0 ; i<n1 ;i++,in+=n1 )
			for ( unsigned int j=i+1 ; j<n1 ; j++)
				Kbeta[j*n1+i] = Kbeta[in+j];

	}
	else
	{
		for ( unsigned k=0,knn=0 ; k<m ; k++,knn+=nn )
		{
			if (beta[k]<=1e-8)
				continue;
			
			for ( unsigned int i=0,in=0 ; i<n2 ; i++,in+=n1 )
				for ( unsigned int j=0 ; j<n1 ; j++ )
					Kbeta[in+j] += beta[k] * K[knn+in+j];
		}
	}

}



