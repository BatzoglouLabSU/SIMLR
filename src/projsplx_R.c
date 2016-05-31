#include <Rmath.h>
#include <R.h>

void projsplx_R(double *y, double *x)
{
    
    int m, n, d, j,npos,ft;
    double *s, *vs;
    double sumResult = -1, tmpValue, tmax, f,lambda_m;
    
    m = sizeof(y);
    n=sizeof(y)/sizeof(y[0]);
    
    
    /*  set the output pointer to the output matrix */
    s = (double*) calloc (m,sizeof(double));
    vs = (double*) calloc (m,sizeof(double));
    for (int k=0;k<n;k++){
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

