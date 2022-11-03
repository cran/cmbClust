#include<math.h>
#include <stdio.h>
#include <stdlib.h>
/*#include <time.h>*/

#include "array.h"
#include "matrix.h"
#include "Order_EM.h"

	#include <Rmath.h>

/* sampling without replacement */

void srswor(int M, int n, int *y){

	int flag;
	int k, v;
	int *indy;

	MAKE_VECTOR(indy, n);

    for (k=0; k<n; k++) indy[k] = 0;

	for (k=0; k<M; k++){

		flag = 0;

		while (flag == 0){

			v = floor(runif(0.0, n));
			if (indy[v] == 0){
				y[k] = v;
				indy[v] = 1;
                flag = 1;
            }
        }

	}

	FREE_VECTOR(indy);

	return;

}



void kmeans1(int n, int p, int K, double **y, int *id){

    int *r_id;
    double dist[K];
    MAKE_1ARRAY(r_id, K);
    // Use current time as
    // seed for random generator
    /*srand(time(0));*/
    srswor(K, n, r_id);



    for(int i = 0; i<n; i++){
            for(int k = 0; k<K; k++){
                 dist[k] = 0;
                 for(int j = 0; j<p; j++) dist[k] += pow(y[i][j]-y[r_id[k]][j], 2);
                 dist[k] = sqrt(dist[k]);
            }
           id[i] = which_min(K, dist);
        }


    FREE_1ARRAY(r_id);


}



void EMEM(int n, int p, int p2, int K, int m, double tol, int n_em, int em_iter, int nk_min, double **y, int *id, double *ll){


   double *bic, **nid;
   int *ascend_order, num_id[K];

   MAKE_1ARRAY(ascend_order, n_em);
   MAKE_1ARRAY(bic, n_em);
   MAKE_2ARRAY(nid, n, n_em);



    for(int i = 0; i<n_em; i++){

       kmeans1(n, p, K, y, id);
       EM_all(n, p, p2, K, m, tol, em_iter, y, id, ll);
       bic[i] = ll[1];
       for(int j = 0; j<n; j++) nid[j][i] = id[j];

    }


    for(int i = 0; i<n_em; i++) ascend_order[i] = i;
    ordervector(n_em, bic, ascend_order);



    for(int j = 0; j<n_em; j++){
        for(int k=0; k <K; k++){
           num_id[k] = 0;
           for(int i = 0; i<n; i++){
            if(nid[i][ascend_order[j]]==k) num_id[k] += 1;
           }
           /*printf("%d ", num_id[k]);*/
        }
        /*printf("%d \n", miny(K, num_id));*/

        if(miny(K, num_id) > nk_min){
            for(int i = 0; i<n; i++) {id[i] = nid[i][ascend_order[j]];}
            break;
        }

    }



    FREE_2ARRAY(nid);
    FREE_1ARRAY(bic);
    FREE_1ARRAY(ascend_order);
}


