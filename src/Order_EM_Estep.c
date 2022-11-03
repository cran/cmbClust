#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "array.h"
#include "matrix.h"
#include "em.h"
#define pi 3.14159265358979323846





double dmvnorm(int n1, int p1, int K1, double **y, double ***s2_2, double ***mu){

   double f, *detS, **sigma, **inv_sigma;
   double f1, f2 = 0;
   MAKE_1ARRAY(detS, 1);
   MAKE_2ARRAY(sigma, p1, p1);
   MAKE_2ARRAY(inv_sigma, p1, p1);

   for(int j1= 0; j1<p1; j1++){
      for(int j2 = 0; j2<p1; j2++){

          sigma[j1][j2] = s2_2[j1][j2][K1];
      }
   }


   inverse_det(p1, sigma, inv_sigma, detS);
   double detss = *detS;


   f = 1/sqrt(pow(2*pi, p1)*detss);

   for(int j1 = 0; j1< p1; j1++ ){
        f1 = 0;
        for(int j2 = 0; j2< p1; j2++){
             f1 += mu[n1][j2][K1]*inv_sigma[j1][j2];
        }
        f2 += f1*mu[n1][j1][K1];
   }

   f *= exp(-0.5*f2);


   FREE_2ARRAY(sigma);
   FREE_2ARRAY(inv_sigma);
   FREE_1ARRAY(detS);

   return f;

}


double density_(int n1, int p, int p0, int K1, int m, double **y, double **beta1, double **s2_1, double ***s2_2, double ***mu){
   int i, j=0, z, t = 0;

   double mu1 = beta1[K1][0] ;
   double f = dnorm(n1, j, K1, y, mu1, s2_1);

   for(j= 1; j < p0; j++){
        t++;
        mu1 = beta1[K1][t];
        for(i = 0; i < j; i++){
            for(z = 0; z < m; z++){
                  t++;
                  mu1 += pow(y[n1][i], (z + 1))*beta1[K1][t];
            }
        }
       f *=   dnorm(n1, j, K1, y, mu1, s2_1);
    }


    f *= dmvnorm(n1, p-p0, K1, y, s2_2, mu);


   return f;

 }


void classprob_(int n1, int p, int p0, int K, int m, double **y, double **beta1, double **s2_1, double ***s2_2, double ***mu, double *tau, double ** class_prob){
   int k;
   double sum_pi = 0.0;

   for(k = 0; k< K; k++){
        class_prob[n1][k]  = tau[k]*density_(n1, p, p0, k, m, y, beta1, s2_1, s2_2, mu);
        sum_pi += class_prob[n1][k];
        }

   for(k = 0; k< K; k++){
        class_prob[n1][k]  /= sum_pi;

    }


}


void update_class_prob_(int n, int p, int p0, int K, int m, double **y, double **beta1, double **s2_1, double ***s2_2, double ***mu, double * tau, double ** class_prob){
    int i;

     for (i = 0; i < n; i++){
            classprob_(i, p, p0, K, m, y, beta1, s2_1, s2_2, mu, tau, class_prob);
     }


}
