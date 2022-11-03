
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "array.h"
#include "matrix.h"
#include "em.h"
#define pi 3.141593



void EM(int n, int p, int K, int m, int nbeta, double tol, int niter, double **y, double **sd, int **indicator, double **class_prob, double **beta, double *tau, int *id, double *ll){

    int i, j, ii, k, z, sum_indicator = 0, t = 0, zeo = 0;
    int **sub_indicator;
    double prev_llk, llk = 0, x = 0;
    MAKE_2ARRAY(sub_indicator, K, p);




    for(k = 0; k < K; k++){
        t = 0;
        for(j = 0; j < p; j++){
            sub_indicator[k][j] = 0;
            for(z = 0; z < 1+j*m; z++){
                sub_indicator[k][j] += indicator[k][t];
                    t++;
           }
           sum_indicator += sub_indicator[k][j];
          }
    }




    prev_llk = llk;


    /*clock_t st;
    st = clock();*/
    update_parameters(n, p, K, m, y, sd, indicator, sub_indicator, beta, tau, class_prob);
    /*st = clock() - st;
    double time_taken = ((double)st)/CLOCKS_PER_SEC; // in seconds
    printf("Mstep takes %f seconds\n",  time_taken);*/

    llk = mixLLK(n, p, K, m, y, sd, indicator, beta, tau);
    /*printf("llk is %f\n", llk);*/

    ii = 1;
    do{
         ii++;
         /*st = clock();*/
         update_class_prob(n, p, K, m, y, sd, indicator, beta, tau, class_prob);
         /*st = clock() - st;
         time_taken = ((double)st)/CLOCKS_PER_SEC; // in seconds
         printf("Estep takes %f seconds\n", time_taken);*/

         /*st = clock();*/
         update_parameters(n, p, K, m, y, sd, indicator, sub_indicator, beta, tau, class_prob);
         /*st = clock() - st;
         time_taken = ((double)st)/CLOCKS_PER_SEC; // in seconds
         printf("Mstep takes %f seconds\n",  time_taken);*/


         prev_llk = llk;
         llk = mixLLK(n, p, K, m, y, sd, indicator, beta, tau);
        /*printf("%d, ", ii);
         printf(" llk is %f\n", llk);*/
         if(ii == niter) break;
      }while (check_tol(llk, prev_llk, tol)==0&&(!isnan(llk)));


    ll[0] = llk;
    ll[1] = -2*llk + ((K-1) + sum_indicator + p*K)*log(n);
    /*printf("number of iteration is %d\n", ii);
    printf("llk is %f\n", llk);
    printf("bic is %f\n", ll[1]);*/
    ll[2] = ((K-1) + sum_indicator + p*K);

    for(i = 0; i<n; i++){
        x = class_prob[i][0];
        id[i] = zeo;
        for(k=1; k<K; k++){
            if(x < class_prob[i][k]){
                    id[i] = k;
                    x = class_prob[i][k];
                    }

        }
    }



     FREE_2ARRAY(sub_indicator);

}
