#include <R.h>
#include <Rmath.h>

/*
 * This function computes
 * A = 0;
 * for i=1:length(beta), A = A + K(:,:,i)*beta(i); end
 * for symmetric and not-symmetric K
 */    
void KbetaR(double *K, double *beta, double *Kbeta)
{

	

	unsigned int n1 = sizeof(Kbeta);
	unsigned int n2 = sizeof(Kbeta)/sizeof(Kbeta[0]);
	unsigned int m = sizeof(beta);
	unsigned int nn = n1*n2;

	

	bool symmetric = true;

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



