#ifndef MATRIX_H
#define MATRIX_H

void xx_product(double **x, double ncol, double nrow, double **xx);
void ax_product(double **a, double nrow_a, double **x, double nrow_x, double **ax);
void ay_product(double **a, double ncol_a, double nrow_a, double *y_vector, double *ay);
void array1to2(int a, int b, double *y, double **x);
void array2to1(int a, int b, double *y, double **x);
void array1to2i(int a, int b, int *y, int **x);
void array2to1i(int a, int b, int *y, int **x);


void XAXt(double **X, int p, double **A, double **Res);
void multiply(double **a, int arows, int acols, double **b, int brows, int bcols, double **c);
void tA(double **A, int a, int b, double **Res);
void ItoJ(int n, int j1, int j2, double **y);
void Anull(double **X, int ax, int bx);
double det(int p, double **S);
void inverse(int p, double **S, double **Sinv);
void inverse_det(int p, double **S, double **Sinv, double (*detS));
int which_min(int p, double *y);
int miny(int p, int *y);
void ordervector(int p, double *y, int *r);
void classprobTOid(int n, int K,  double **class_prob, int *id);
void idTOclassprob(int n, int K,  int *id, double **class_prob);
int Factorial(int a);

#endif /* MATRIX_H */
