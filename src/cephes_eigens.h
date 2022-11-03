#ifndef CEPHES_EIGENS_H
#define CEPHES_EIGENS_H


void cephes_eigens(double *A, double *EV, double *E, int N);
void cephes_symmeigens_down(int p, double *eval, double **A, double (*determinant));

#endif /* CEPHES_EIGENS_H */
