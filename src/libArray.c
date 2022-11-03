#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.141593

void array1to2(int a, int b, double *y, double **x){

	int i, j, k;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){

			x[i][j] = y[k];
			k++;

		}
	}


}



void array2to1(int a, int b, double *y, double **x){

	int i, j, k;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){

			y[k] = x[i][j];
			k++;

		}
	}


}


void array1to2i(int a, int b, int *y, int **x){

	int i, j, k;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){

			x[i][j] = y[k];
			k++;

		}
	}


}


void array2to1i(int a, int b, int *y, int **x){

	int i, j, k;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){

			y[k] = x[i][j];
			k++;

		}
	}


}
