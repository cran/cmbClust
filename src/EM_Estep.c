#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "array.h"
#include "matrix.h"

#define pi 3.14159265358979323846

double dnorm(int n1, int p1, int K1, double **y, double mean, double **sd){
   double f;
   f = (1/sqrt(2*pi*sd[K1][p1]))*exp(-pow((y[n1][p1]-mean),2)/(2*sd[K1][p1]));
   return f;

}

double density(int n, int p, int m, int n1, int K1, double **y, double **sd, int **indicator, double **beta){
   int i, j=0, z, t = 0;
   double mu=0.0, f = dnorm(n1, j, K1, y, beta[K1][0], sd);
   for(j= 1; j<p; j++){
        t++;
        mu = beta[K1][t];
            for(i = 0; i < j; i++){
                for(z = 0; z < m; z++){
                  t++;
                  mu += pow(y[n1][i], (z + 1))*beta[K1][t];
            }
        }
       f *=   dnorm( n1, j, K1, y, mu, sd);
    }
   return f;

 }


void classprob(int n, int p, int K, int m, int n1, double **y, double **sd, int **indicator, double **beta, double * tau, double **class_prob){
   int k;
   double sum_pi = 0.0;

   for(k = 0; k< K; k++){
        class_prob[n1][k]  = tau[k]*density(n, p, m, n1, k, y, sd, indicator, beta);
        sum_pi += class_prob[n1][k];
        }

   for(k = 0; k< K; k++){
        class_prob[n1][k]  /= sum_pi;
    }


}


void update_class_prob(int n, int p, int K, int m, double **y, double **sd, int **indicator, double **beta, double * tau, double ** class_prob){
    int i;

     for (i = 0; i < n; i++){
            classprob(n, p, K, m, i, y, sd, indicator, beta, tau, class_prob);
     }

}
