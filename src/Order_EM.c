#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "array.h"
#include "Order_EM.h"
#include "em.h"
#include "matrix.h"



void EM_all(int n, int p, int p1, int K, int m, double tol, int EM_iter, double **y, int *id, double *ll){



    /*printf("%d -th iteration \n", p1);*/
    double **s2_1, ***s2_2, **beta1, ***beta2, ***mu, *tau, **class_prob, llk, prev_llk;
    int ii = 1, stop = 0;

    MAKE_2ARRAY(s2_1, K, p1+1);
    MAKE_2ARRAY(beta1, K, p1+1 + (p1+1)*p1*m/2);
    MAKE_3ARRAY(s2_2,  p-p1-1, p-p1-1, K);
    MAKE_3ARRAY(beta2, p-p1-1, m*(p1+1)+1, K);
    MAKE_3ARRAY(mu, n, p-p1-1, K);
    MAKE_1ARRAY(tau, K);
    MAKE_2ARRAY(class_prob, n, K);


   idTOclassprob(n, K,  id, class_prob);
   updata_tau(n, K, tau, class_prob);
   update_parameters1(n, p1+1, K, m, y, s2_1, beta1, tau, class_prob);
   update_parameters2(n, p, p1+1, K, m, y, s2_2, beta2, mu, tau, class_prob);
   llk = mixLLK_(n, p, p1+1, K, m, y,  beta1, s2_1, s2_2, mu, tau, class_prob);
   /*printf("%d, ", 1);
   printf("llk is %f \n", llk);*/


   do{
         ii++;

         update_class_prob_(n, p, p1+1, K, m, y, beta1, s2_1, s2_2, mu, tau, class_prob);
        for(int i = 0; i<n; i++){
            for(int k = 0; k<K; k++){
                if(isnan(class_prob[i][k])) stop = 1;
            }
         }
         if(stop ==1) break;
         updata_tau(n, K, tau, class_prob);

         update_parameters1(n, p1+1, K, m, y, s2_1, beta1, tau, class_prob);

         update_parameters2(n, p, p1+1, K, m, y, s2_2, beta2, mu, tau, class_prob);

         prev_llk = llk;
         llk = mixLLK_(n, p, p1+1, K, m, y,  beta1, s2_1, s2_2, mu, tau, class_prob);

         /*printf("%d, ", ii);
         printf(" llk is %f\n", llk);*/
         if(ii == EM_iter) break;
          }while (check_tol_(llk, prev_llk, tol)==0&&(!isnan(llk)));


    ll[0] = llk;
    ll[1] = BIC(n, p, p1+1, K, m, llk);
    classprobTOid(n, K, class_prob, id);



    FREE_2ARRAY(s2_1);
    FREE_2ARRAY(beta1);
    FREE_3ARRAY(s2_2);
    FREE_3ARRAY(beta2);
    FREE_3ARRAY(mu);
    FREE_1ARRAY(tau);
    FREE_2ARRAY(class_prob);



}




