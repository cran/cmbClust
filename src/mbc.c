#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "array.h"
#include "matrix.h"
#include "em.h"
#define pi 3.141593

void mbc(double *y1, int (*n1), int (*K1), int (*p1), int (*m1), int *id, double *ll, double *tau, int *indicator1, double *beta1, double *sd1, double *class_prob1, int (*niter1), double (*tol1)){

    int  j, k, t = 0;
    int n, K, p, m, niter;
	double tol;
	int **indicator;
	double **y, **sd, **beta, **class_prob;


	n = (*n1);
	K = (*K1);
	p = (*p1);
	m = (*m1);
	tol = (*tol1);
	niter = (*niter1);
	int nbeta = p + (p-1)*p*m/2;


    MAKE_2ARRAY(indicator, K, nbeta);
	MAKE_2ARRAY(y, n, p);
	MAKE_2ARRAY(sd, K, p);
	MAKE_2ARRAY(beta, K, nbeta);
	MAKE_2ARRAY(class_prob, n, K);


    array1to2i(K, nbeta, indicator1, indicator);
	array1to2(n, p, y1, y);
	array1to2(K, p, sd1, sd);
	array1to2(n, K, class_prob1, class_prob);



    for(k = 0; k < K; k++){
        for(j = 0; j < nbeta; j++){
            if(indicator[k][j] ==1){
             beta[k][j] = beta1[t];
            }else{beta[k][j] = 0;
            }
            t++;
        }
    }



    clock_t st;
    st = clock();
    EM(n, p, K, m, nbeta, tol, niter, y, sd, indicator, class_prob, beta, tau, id, ll);
    st = clock() - st;
    /*double time_taken = ((double)st)/CLOCKS_PER_SEC; // in seconds
    printf("EM takes %f seconds\n",  time_taken);*/


    array2to1i(K, nbeta, indicator1, indicator);
	array2to1(n, p, y1, y);
	array2to1(K, p, sd1, sd);
    array2to1(K, nbeta, beta1, beta);
    array2to1(n, K, class_prob1, class_prob);

    FREE_2ARRAY(indicator);
    FREE_2ARRAY(y);
    FREE_2ARRAY(sd);
    FREE_2ARRAY(beta);
    FREE_2ARRAY(class_prob);

}
