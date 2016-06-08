#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

SEXP projsplx_R(SEXP y, SEXP x)
{

    int m,n,j,npos,ft;
    SEXP s, vs;
    double f,lambda_m;
    SEXP Rdim = getAttrib(y, R_DimSymbol);
    m = INTEGER(Rdim)[0];
    n = INTEGER(Rdim)[1];
    
    s = allocMatrix(REALSXP, m, 1);
    vs = allocMatrix(REALSXP, m, 1);
    
    for (int k=0;k<n;k++){
        
        double means = 0;
        double mins = 100000;
        for(j = 0; j < m; j++ ){
            REAL(s)[j] = REAL(y)[j+k*m];
            means += REAL(s)[j];
            mins = (mins > REAL(s)[j])? REAL(s)[j]:mins;
        }

        for(j = 0; j < m; j++ ){
            REAL(s)[j] -= (means-1)/m;
        }
        ft=1;
        if(mins<0){
            f = 1;
            lambda_m=0;
            while(fabs(f)>1e-10){
                npos = 0;
                f = 0;
                for(j = 0; j < m; j++ ){
                    REAL(vs)[j] = REAL(s)[j]-lambda_m;
                 
                    if (REAL(vs)[j]>0){
                        npos+=1;
                        f+=REAL(vs)[j];
                    }
                }
                lambda_m += (f-1)/npos;
                if(ft>100){
                    for(j = 0; j <= m-1; j++){
                        REAL(x)[j+k*m] = (REAL(vs)[j] > 0)? REAL(vs)[j]:0;
                    }
                    break;
                }
                ft+=1;
            }
            for(j = 0; j <= m-1; j++){
                REAL(x)[j+k*m] = (REAL(vs)[j] > 0)? REAL(vs)[j]:0;
            }
            
        }
        else{
            for(j = 0; j <= m-1; j++){
                REAL(x)[j+k*m] = REAL(s)[j];
            }
            
        }
    }

    return x;
}
