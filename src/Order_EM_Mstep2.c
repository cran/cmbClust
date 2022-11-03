#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "array.h"
#include "matrix.h"


#define pi 3.141593


/*=====================================================================================================================*/


void x_matrix2(int n, int p, int p0, int K1, int m, double **x1, double **x2, double **y){

    int i, z, j, t;
    double one = 1;

    for(i = 0; i < n; i++){
         t = 0;
         x1[i][t] =  one;

         for(j = 0; j < p0; j++){
            for(z = 0; z < m; z++){
                   t++;
                   x1[i][t] = pow(y[i][j], (z+1.0));

               }
        }

        for(j = 0; j<(p-p0); j++){
               x2[i][j] = y[i][j+p0];

        }

    }


}



void update_beta_(int n, int p, int p0, int K1, int nbeta, double **x1, double **x2, double ***beta2, double **class_prob){

    int i, j1, j2;
    double **xx, **xy, **inv_xx, zero;
    zero = 0.0;
    MAKE_2ARRAY(xx, nbeta, nbeta);
    MAKE_2ARRAY(inv_xx, nbeta, nbeta);
    MAKE_2ARRAY(xy, nbeta, p-p0);


    for(j1 = 0; j1< nbeta; j1++){
        for(j2 = j1; j2 < nbeta; j2++){
            xx[j1][j2] = zero;
            for(i = 0; i < n; i++){
               xx[j1][j2] += (class_prob[i][K1]*x1[i][j1]*x1[i][j2]);
            }
            xx[j2][j1] = xx[j1][j2];
       }

        for(j2 = 0; j2 < p-p0; j2++){
            xy[j1][j2] = zero;
            for(i = 0; i < n; i++){
                xy[j1][j2] += (class_prob[i][K1]*x1[i][j1]*x2[i][j2]);
             }
        }


    }

   inverse(nbeta, xx, inv_xx);

    for (int j = 0; j < p-p0; j++){
        for (i = 0; i< nbeta; i++){
            beta2[j][i][K1] = zero;
            for (int k=0; k < nbeta; k++)  {beta2[j][i][K1] += inv_xx[i][k] * xy[k][j];}
        }
    }


   FREE_2ARRAY(xy);
   FREE_2ARRAY(xx);
   FREE_2ARRAY(inv_xx);
}



void updata_sd_(int n, int p, int p0, int K1, int nbeta, double **x1, double **x2, double ***beta2, double ***s2_2, double ***mu, double *tau, double **class_prob){

   int i ;
   double  zero = 0.0;




   for(i = 0; i<n; i++){
        for(int j1 = 0; j1 < p-p0; j1++){
            mu[i][j1][K1] = zero;
            for(int j2 = 0; j2 < nbeta; j2++){
                mu[i][j1][K1] += beta2[j1][j2][K1]*x1[i][j2];
             }
            mu[i][j1][K1] = x2[i][j1]-mu[i][j1][K1];
        }
    }

    for(int j1 = 0; j1< p - p0; j1++){
        for(int j2 = j1; j2 < p- p0; j2++){
             s2_2[j1][j2][K1] = zero;
             for(i = 0; i<n; i++){
                    s2_2[j1][j2][K1] += (class_prob[i][K1]*mu[i][j1][K1]*mu[i][j2][K1]/(n*tau[K1]));
                  }
             s2_2[j2][j1][K1] = s2_2[j1][j2][K1];

             }
    }




}



void update_beta_sd_(int n, int p, int p0, int K1, int m, double **y, double ***s2_2, double ***beta2, double ***mu, double *tau, double **class_prob){


    int nbeta = m*p0+1;
    double **x1, **x2;
    MAKE_2ARRAY(x1, n, nbeta);
    MAKE_2ARRAY(x2, n, p-p0);

    x_matrix2(n, p, p0, K1, m, x1, x2, y);

    update_beta_(n, p, p0, K1, nbeta, x1, x2, beta2, class_prob);
    updata_sd_(n, p, p0, K1, nbeta, x1, x2, beta2, s2_2, mu, tau, class_prob);


    FREE_2ARRAY(x1);
    FREE_2ARRAY(x2);
}



/*=====================================================================================================================*/
void update_parameters2(int n, int p, int p0, int K, int m, double **y, double ***s2_2, double ***beta2, double ***mu, double *tau,  double **class_prob){

    int k;

    for(k = 0; k < K; k++){
       update_beta_sd_(n, p, p0, k,  m, y, s2_2, beta2, mu, tau, class_prob);
    }

}







/*=====================================================================================================================*/





