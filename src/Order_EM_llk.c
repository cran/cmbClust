#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "matrix.h"
#include "Order_EM.h"
#define pi 3.141593



double mixLLK_(int n, int p, int p0, int K, int m, double **y,  double **beta1, double **s2_1, double ***s2_2, double ***mu, double * tau, double ** class_prob){

    int i, k;
    double llk = 0, f;

    for(i =0; i<n; i++){
        f = 0;
        for(k = 0; k<K; k++){
           f += tau[k]*density_(i, p, p0, k, m, y, beta1, s2_1, s2_2, mu);
        }

        llk += log(f);
     }

    return llk;
}


double BIC(int n, int p, int p1, int K, int m, double llk){

    double bic;
    int df1, df2;
    df1 = (p1+p1*(p1-1)*m/2 + p1)*K ;
    df2 = K*(p1+1)*(p-p1)+K*(p-p1)*(p-p1+1)/2;

    bic = -2*llk + (K-1 + df1+ df2)*log(n);

    return bic;
}

int check_tol_(double llk, double prev_llk, double eps){
   int stop;
   if((llk-prev_llk)/fabs(llk)<eps) {
    stop = 1;
   }else      stop = 0;

   return stop;
}

