
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "array.h"
#include "matrix.h"


#define pi 3.141593


/*=====================================================================================================================*/


void update_beta_sd1_(int n, int p_1, int K1, int m, double **y, double **s2_1, double **beta1, double *tau, double **class_prob){
    int i;
    double  zero = 0;
    beta1[K1][p_1] = zero;
    s2_1[K1][p_1] = zero;

    for ( i = 0; i < n; i++) beta1[K1][p_1] += y[i][p_1] * class_prob[i][K1];
    beta1[K1][0] /= (n * tau[K1]);


    for(i = 0; i<n; i++) s2_1[K1][p_1] += (pow((y[i][p_1] - beta1[K1][p_1]), 2.0)*class_prob[i][K1]/(n*tau[K1]));



    /*printf("%d and %f \n", K1,  beta[K1][p_1+(p_1-1)*p_1*m/2]);*/
}

/*=====================================================================================================================*/


void x_matrix1(int n, int p_1, int K1, int m, double **x, double **y){

    int i, z, j, t;
    double one = 1;

    for(i = 0; i < n; i++){
         t = 0;
         x[i][t] =  one;
         for(j = 0; j < p_1; j++){
             for(z = 0; z < m; z++){
                   t++;
                   x[i][t] = pow(y[i][j], (z+1.0));
               }
        }
    }


}


void update_beta2_(int n, int p_1, int K1, int ncol_x, double **x, double **y, double **class_prob, double *beta_2){

    int i, j1, j2;
    double **xx, *xy, **inv_xx, zero;
    zero = 0.0;
    MAKE_2ARRAY(xx, ncol_x, ncol_x);
    MAKE_2ARRAY(inv_xx, ncol_x, ncol_x);
    MAKE_1ARRAY(xy, ncol_x);

    for(j1 = 0; j1< ncol_x; j1++){
        xy[j1] = zero;
        for(j2 = j1; j2 < ncol_x; j2++){
            xx[j1][j2] = zero;
            for(i = 0; i < n; i++){
               xx[j1][j2] += (class_prob[i][K1]*x[i][j1]*x[i][j2]);
               if(j1 == j2) xy[j1] += (class_prob[i][K1]*x[i][j1]*y[i][p_1]);
            }
            xx[j2][j1] = xx[j1][j2];
       }
    }

   inverse(ncol_x, xx, inv_xx);


   for(j1 = 0; j1< ncol_x; j1++){
         beta_2[j1] = zero;
        for(j2 = 0; j2 < ncol_x; j2++){
                beta_2[j1] += (inv_xx[j1][j2]*xy[j2]);
            }

    }


   FREE_1ARRAY(xy);
   FREE_2ARRAY(xx);
   FREE_2ARRAY(inv_xx);
}


double updata_sd2_(int n, int p_1, int K1, int ncol_x, double **x, double **y, double *beta_2, double *tau, double **class_prob){

   int i ,j;
   double mu, v = 0.0;

   for(i = 0; i<n; i++){
      mu = 0.0;
      for(j = 0; j < ncol_x; j++){
        mu += x[i][j]*beta_2[j];
      }
      v += (pow((y[i][p_1] - mu), 2.0)*class_prob[i][K1]/(n*tau[K1]));
   }

   return v;
}



void update_beta_sd2_(int n, int p_1, int K1, int m, double **y, double **s2_1, double **beta1, double *tau, double **class_prob){

    int i, t, nbeta = 1 + m*p_1;
    double **x, *beta_2;
    MAKE_2ARRAY(x, n, nbeta);
    MAKE_1ARRAY(beta_2, nbeta);

    x_matrix1(n, p_1, K1, m, x, y);
    update_beta2_(n, p_1, K1, nbeta, x, y, class_prob, beta_2);
    s2_1[K1][p_1] = updata_sd2_(n, p_1, K1, nbeta, x, y, beta_2, tau, class_prob);

    for(t = 0; t < nbeta; t++){
          i = p_1 + (p_1-1)*p_1*m/2 + t;
          beta1[K1][i] = beta_2[t];
    }

    FREE_2ARRAY(x);
    FREE_1ARRAY(beta_2);

}



/*=====================================================================================================================*/
void update_parameters1(int n, int p, int K, int m, double **y, double **s2_1, double **beta1, double *tau, double **class_prob){

    int j, k;



    for(k = 0; k < K; k++){
        update_beta_sd1_(n, 0, k, m, y, s2_1, beta1, tau, class_prob);
        for(j = 1; j < p; j++){
           update_beta_sd2_(n, j, k, m, y, s2_1, beta1, tau, class_prob);
        }
    }


}



/*=====================================================================================================================*/





