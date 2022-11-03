#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "array.h"
#include "matrix.h"
#include "em.h"
#define pi 3.141593



double mixLLK(int n, int p, int K, int m, double **y, double **sd, int **indicator, double **beta, double *tau){

    int i, k;
    double llk = 0, f;

    for(i =0; i<n; i++){
        f = 0;
        for(k = 0; k<K; k++){
           f += tau[k]*density(n, p, m, i, k, y, sd, indicator, beta);
        }

        llk += log(f);
     }

    return llk;
}


int check_tol(double llk, double prev_llk, double eps){
   int stop;
   if((llk-prev_llk)/fabs(llk)<eps) {
    stop = 1;
   }else      stop = 0;

   return stop;
}

