#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "array.h"
#include "matrix.h"
#include "Order_EM.h"


void order(double *y1, int (*n1), int (*p1), int (*K1),  int (*m1), int *id, int *orders, double *ll, int (*n_em1),int (*em_iter1), int (*nk_min1), int (*EM_iter1), double (*tol1)){


    int n, K, p, m, EM_iter, n_em, em_iter, nk_min;
	double tol, a;
	double **y;
	int *id1;


	n = (*n1);
	K = (*K1);
	p = (*p1);
	m = (*m1);
	tol = (*tol1);
	EM_iter = (*EM_iter1);
    n_em = (*n_em1);
    em_iter = (*em_iter1);
    nk_min = (*nk_min1);
	MAKE_2ARRAY(y, n, p);
	MAKE_1ARRAY(id1, n);
	array1to2(n, p, y1, y);



    for(int j1 = 0; j1 < (p-1); j1++){
       double bic[p-j1];
       double **id_all;
       MAKE_2ARRAY(id_all, n, p-j1);
       if(j1 == 0) {EMEM(n, p, 0, K, m, tol, n_em, em_iter, nk_min, y, id, ll);}
       for(int i = 0; i<n; i++) id1[i] = id[i];
       EM_all(n, p, j1, K, m, tol, EM_iter, y, id1, ll);
       bic[0] = ll[1];
       for(int i = 0; i<n; i++) id_all[i][0] = id1[i];

       for(int j2 = 1; j2 < (p-j1); j2++){
          ItoJ(n, j1, j2+j1, y);
          if(j1 == 0) {EMEM(n, p, 0, K, m, tol, n_em, em_iter, nk_min, y, id, ll);}
          for(int i = 0; i<n; i++) id1[i] = id[i];
          EM_all(n, p, j1, K, m, tol, EM_iter, y, id1, ll);
          for(int i = 0; i<n; i++) id_all[i][j2] = id1[i];
          bic[j2] = ll[1];
          ItoJ(n, j1, j2+j1, y);

       }

       int z = which_min(p-j1, bic);
       for(int i = 0; i<n; i++) id[i] = id_all[i][z];

       if(z!= 0)
        {a = orders[z+j1];
         orders[z+j1]  = orders[j1];
         orders[j1] = a;
         ItoJ(n, j1, j1+z, y);
        }

      FREE_2ARRAY(id_all);
    }

      /*for(int i = 0; i<p; i++){
        printf("%d  ", orders[i]);
      }
      printf("\n");*/

    FREE_2ARRAY(y);
    FREE_1ARRAY(id1);

}
