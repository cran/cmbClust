#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "array.h"
#include "matrix.h"

#define pi 3.141593


void updata_tau(int n, int K, double *tau,  double ** class_prob){

    int k, i;
    double t;
    for (k = 0; k < K; k++){
          t = 0.0;
          for (i = 0; i < n; i++) {
                t += class_prob[i][k];
          }
          t /= n;
          tau[k] = t;
      }

}


double updata_sd(int n, int K1, int p1, double **y, double *mu, double *tau, double **class_prob){

   int i;
   double v = 0.0;

   for(i = 0; i<n; i++){

    v += pow((y[i][p1] - mu[i]), 2.0)*class_prob[i][K1]/(n*tau[K1]) ;

   }

   return v;
}

/*=====================================================================================================================*/


void update_beta_sd1(int n, int K1, int p1, int m, double **y, double **sd, double **beta, double * tau, double **class_prob){
    int i;
    double beta1 = 0.0;
    double *mu;
    MAKE_1ARRAY(mu, n);

    for ( i = 0; i < n; i++) beta1 += y[i][p1] * class_prob[i][K1];
    beta1 /= n * tau[K1];
    beta[K1][p1+(p1-1)*p1*m/2] = beta1;
    for(i = 0; i < n; i++) mu[i] = beta1;
    sd[K1][p1] = updata_sd(n, K1, p1, y, mu, tau, class_prob);

    FREE_1ARRAY(mu);

}
/*=====================================================================================================================*/

void xy_matrix(int n, int m, int K1, int p1, double **y, int ** indicator, double **class_prob, double **x_matrix, double ** x_m, double *y_vector){

    int i, z, j, t, tt;

    for(i = 0; i < n; i++){
         t = 0;
         tt = 0;
         x_matrix[i][t] =  pow(class_prob[i][K1], 0.5);
         x_m[i][t] =  1;

         for(j = 0; j < p1; j++){
             for(z = 0; z < m; z++){
                   t++;
                   if(indicator[K1][p1 + (p1 -1)*p1*m/2 + t] == 1){
                         tt++;
                         x_matrix[i][tt] = pow(y[i][j], (z+1.0)) * pow(class_prob[i][K1], 0.5);
                         x_m[i][tt] = pow(y[i][j], (z+1.0));
                       }


               }

        }
        y_vector[i] = y[i][p1]*pow(class_prob[i][K1], 0.5);

    }

}


void update_beta(int ncol_x, int n, double **x_matrix, double *y_vector, double* beta2){

    double **xx_matrix, **inv_matrix, **ax;
    MAKE_2ARRAY(xx_matrix, ncol_x, ncol_x);
    MAKE_2ARRAY(inv_matrix, ncol_x, ncol_x);
    MAKE_2ARRAY(ax, ncol_x, n);

    /*clock_t st = clock();*/
    xx_product(x_matrix, ncol_x, n, xx_matrix);
    /*st = clock() - st;
    double time_taken = ((double)st)/CLOCKS_PER_SEC; // in seconds
    printf("XXproduct takes %f seconds\n",  time_taken);*/

    /*st = clock();*/
    inverse(ncol_x, xx_matrix, inv_matrix);
    /*st = clock() - st;
    time_taken = ((double)st)/CLOCKS_PER_SEC; // in seconds
    printf("XXinverse takes %f seconds\n",  time_taken);*/




    ax_product(inv_matrix, ncol_x, x_matrix, n, ax);
    ay_product(ax, ncol_x, n, y_vector, beta2);



    FREE_2ARRAY(xx_matrix);
    FREE_2ARRAY(inv_matrix);
    FREE_2ARRAY(ax);


}

void updata_mu(int n, int num, double **x_m, double *beta2, double *mu){
   int i, j;

   for(i = 0; i < n; i++){
      mu[i] = 0;
      for(j = 0; j < num; j++){
        mu[i] += x_m[i][j]*beta2[j];
      }
   }
}



void update_beta_sd2(int n, int m, int K1, int p1, double **y, double ** sd, int **indicator, int **sub_indicator, double **beta, double * tau, double **class_prob){

    int i, t, tt=0;
    double **x_matrix, **x_m, *y_vector, *beta2, *mu;

    MAKE_2ARRAY(x_matrix, n, sub_indicator[K1][p1]);
    MAKE_2ARRAY(x_m, n, sub_indicator[K1][p1]);

    MAKE_1ARRAY(y_vector, n);
    MAKE_1ARRAY(beta2, sub_indicator[K1][p1]);
    MAKE_1ARRAY(mu, n);


    xy_matrix(n, m, K1, p1, y, indicator, class_prob, x_matrix, x_m, y_vector);
    update_beta(sub_indicator[K1][p1], n, x_matrix, y_vector, beta2);


    updata_mu(n, sub_indicator[K1][p1], x_m, beta2, mu);
    sd[K1][p1] = updata_sd(n, K1, p1, y, mu, tau, class_prob);


    for(t = 0; t <(1+p1*m); t++){
          i = p1 + (p1-1)*p1*m/2 + t;
          if(indicator[K1][i] == 1){
                    beta[K1][i] = beta2[tt];
                    tt++;
          }

    }





    FREE_2ARRAY(x_matrix);
    FREE_1ARRAY(y_vector);
    FREE_1ARRAY(beta2);
    FREE_1ARRAY(mu);
    FREE_2ARRAY(x_m);

    }





/*=====================================================================================================================*/

void update_parameters (int n, int p, int K, int m, double **y, double ** sd, int **indicator, int ** sub_indicator, double **beta, double * tau,  double ** class_prob){

    int j, k;

    updata_tau(n, K, tau, class_prob);

    for(k = 0; k < K; k++){
        for(j = 0; j < p; j++){
               if(sub_indicator[k][j] ==1){
                    update_beta_sd1(n, k, j, m, y, sd, beta, tau, class_prob);
               }else{
                    update_beta_sd2(n, m, k, j, y, sd, indicator, sub_indicator, beta, tau, class_prob);
               }

        }
    }



}








