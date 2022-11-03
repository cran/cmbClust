#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "array.h"
#include "cephes_eigens.h"
#define pi 3.141593

/*  matrix product  x' %*% x*/
void xx_product(double **x, double ncol, double nrow, double **xx){
   int i, j1, j2;
   double zeo = 0.0;

      for(j1 = 0; j1 < ncol; j1++){
            for(j2 = 0; j2 < ncol; j2++){
                  xx[j1][j2] = zeo;
                  for(i = 0; i < nrow; i++){
                    xx[j1][j2] += (x[i][j1] * x[i][j2]);
                  }
      }
   }
}


/*  matrix product  a %*% x' kbyk nbyk*/
void ax_product(double **a, double nrow_a, double **x, double nrow_x, double **ax){
   int i, j1, j2;
   double zeo = 0;

       for(j1 = 0; j1 < nrow_a; j1++){
            for(j2 = 0; j2 < nrow_x; j2++){
                 ax[j1][j2] = zeo;
                  for(i = 0; i < nrow_a; i++){
                        ax[j1][j2] += a[j1][i]*x[j2][i];
                  }

      }

   }

}

/*  matrix product  a %*% y   */
void ay_product(double **a, double nrow_a, double ncol_a, double *y_vector, double *ay){
   int i, j1;
   double zeo = 0;

       for(j1 = 0; j1 < nrow_a; j1++){
            ay[j1] = zeo;
            for(i = 0; i < ncol_a; i++){
                 ay[j1] += a[j1][i]*y_vector[i];
            }
      }

}


void tA(double **A, int a, int b, double **Res){

	int i,j;

   	for (i=0; i<a; i++){
		for (j=0; j<b; j++){
			Res[i][j] = A[j][i];
		}
	}

}



void multiply(double **a, int arows, int acols, double **b, int brows, int bcols, double **c){
  int i, j, k;

  for (i=0; i<arows; i++)
    for (j=0; j<bcols; j++) {
      c[i][j] = 0;
      for (k=0; k<acols; k++)
	c[i][j] += a[i][k] * b[k][j];
    }

}


void XAXt(double **X, int p, double **A, double **Res){

	double **Res1, **Res2;

	MAKE_MATRIX(Res1, p, p);
	MAKE_MATRIX(Res2, p, p);

	tA(X, p, p, Res2);

	multiply(X, p, p, A, p, p, Res1);
	multiply(Res1, p, p, Res2, p, p, Res);

	FREE_MATRIX(Res1);
 	FREE_MATRIX(Res2);

}


/* switch j1 and j2 columns of matrix y */
void ItoJ(int n, int j1, int j2, double **y){

    double a;
    for(int i = 0; i < n; i++){
              a = y[i][j1];
              y[i][j1] = y[i][j2];
              y[i][j2] = a;
     }

}




/* provides 0-matrix */

void Anull(double **X, int ax, int bx){

     int i, j;

     for (i=0; i<ax; i++){
         for (j=0; j<bx; j++){
		X[i][j] = 0.0;
	 }
     }
}


void inverse(int p, double **S, double **Sinv){


int i, j;
double  detS;
double *Eig, **L;


MAKE_MATRIX(L, p, p);
MAKE_VECTOR(Eig, p);


double **S_eigV;
MAKE_2ARRAY(S_eigV, p, p);

for(i=0; i<p; i++){
    for(j= 0; j<p; j++){
        S_eigV[i][j] = S[i][j];
    }
}

cephes_symmeigens_down(p, Eig, S_eigV, &detS);


Anull(L, p, p);

for (j=0; j<p; j++){
L[j][j] = 1 / Eig[j];
}

XAXt(S_eigV, p, L, Sinv);



FREE_MATRIX(L);
FREE_VECTOR(Eig);
FREE_2ARRAY(S_eigV);
}



void inverse_det(int p, double **S, double **Sinv, double (*detS)){


int i, j;
double *Eig, **L;


MAKE_MATRIX(L, p, p);
MAKE_VECTOR(Eig, p);


double **S_eigV;
MAKE_2ARRAY(S_eigV, p, p);

for(i=0; i<p; i++){
    for(j= 0; j<p; j++){
        S_eigV[i][j] = S[i][j];
    }
}

cephes_symmeigens_down(p, Eig, S_eigV, detS);


Anull(L, p, p);

for (j=0; j<p; j++){
L[j][j] = 1 / Eig[j];
}

XAXt(S_eigV, p, L, Sinv);



FREE_MATRIX(L);
FREE_VECTOR(Eig);
FREE_2ARRAY(S_eigV);
}







double det(int p, double **S){

int i, j;
double detS;

double *Eig;
MAKE_VECTOR(Eig, p);


double **S_eigV;
MAKE_2ARRAY(S_eigV, p, p);

for(i=0; i<p; i++){
    for(j= 0; j<p; j++){
        S_eigV[i][j] = S[i][j];
    }
}

cephes_symmeigens_down(p, Eig, S_eigV, &detS);


FREE_VECTOR(Eig);
FREE_2ARRAY(S_eigV);

return detS;

}

/*which minimum of a vector*/
int which_min(int p, double *y){

    double y1 = y[0];
    int z = 0;

    for(int i = 1; i<p; i++){
        if(y1 > y[i]&&(!isnan(y[i]))){
            y1 = y[i];
            z = i;
        }
    }
 return z;

}

/*minimum of a vector*/
int miny(int p, int *y){

    double y1 = y[0];

    for(int i = 1; i<p; i++){
        if(y1 > y[i]&&(!isnan((float)y[i]))){
            y1 = y[i];
        }
    }

 return y1;

}


/* the ascending order of a vector y */
void ordervector(int p, double *y, int *r){

   double y1;
   int z;

   for(int i = 0; i<p; i++){
      y1 = y[i];
      z = i;
      for(int j = i; j<p; j++){
         if(y1 > y[j]&&(!isnan(y[j]))){
            y1 = y[j];
            z = j;
          }
      }
      int r1 = r[i];
      r[i] = r[z];
      r[z] = r1;
      y[z] = y[i];
      y[i] = y1;


   }

   /*for(int i = 0; i<p; i++) printf("%f ", y[i]);
   printf("\n");
   for(int i = 0; i<p; i++) printf("%d ", r[i]);

    printf("\n");*/

}


void classprobTOid(int n, int K, double **class_prob, int *id){
    double x;
    int zero = 0;

    for(int i = 0; i<n; i++){
        x = class_prob[i][0];
        id[i] = zero;
        for(int k=1; k<K; k++){
            if(x < class_prob[i][k]){
                    id[i] = k;
                    x = class_prob[i][k];
                    }

        }
    }

}

void idTOclassprob(int n, int K, int *id, double **class_prob){

    double one = 1, zero = 0;

    for(int i = 0; i<n; i++){
	      for(int k = 0; k<K; k++){
              if(k == id[i]){
              class_prob[i][k] = one;
            }else{class_prob[i][k] = zero;}
        }
    }


}


int Factorial(int a){
    int i;
    int res;

    res=1;
    for (i=1; i<(a+1); i++){
        res=res*i;
    }

    return res;
}
