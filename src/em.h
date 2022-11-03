#ifndef EM_H
#define EM_H
void mbc(double *y1, int (*n1), int (*K1), int (*p1), int (*m1), int *id, double *ll, double *tau, int *indicator1, double *beta1, double *sd1, double *class_prob1, int (*niter1), double (*tol1));
void EM(int n, int p, int K, int m, int nbeta, double tol, int niter, double **y, double **sd, int **indicator, double **class_prob, double **beta, double *tau, int *id, double *ll);
int check_tol(double llk, double prev_llk, double eps);
double mixLLK(int n, int p, int k, int m, double **y, double **sd, int **indicator, double **beta, double *tau);

/*
void updata_tau(int K, int n, double * tau,  double ** class_prob);
double updata_sd(int n, int K1, int p1, double **y, double *mu, double *tau, double **class_prob);
void update_parameters (int n, int p, int K, int m, double **y, double ** sd, int **indicator, int ** sub_indicator, double **beta, double * tau,  double ** class_prob);
void update_beta_sd1(int n, int K1, int p1, int m, double **y, double **sd, double **beta, double * tau, double **class_prob);
void xy_matrix(int n, int m, int K1, int p1, double **y, int ** indicator, double **class_prob, double **x_matrix, double ** x_m, double *y_vector);
void update_beta(int ncol_x, int n, double **x_matrix, double *y_vector, double* beta_vector);
void updata_mu(int n, int num, double **x_matrix, double *beta2, double *mu);
void update_beta_sd2(int n, int m, int K1, int p1, double **y, double ** sd, int **indicator, int **sub_indicator, double **beta, double * tau, double **class_prob);
*/


void update_parameters(int n, int p, int K, int m, double **y, double ** sd, int **indicator, int ** sub_indicator, double **beta, double * tau,  double ** class_prob);
void update_beta_sd2(int n, int m, int K1, int p1, double **y, double ** sd, int **indicator, int **sub_indicator, double **beta, double * tau, double **class_prob);
double updata_sd2(int n, int K1, int p1, int ncol_x, double **x, double **y, double *beta2, double *tau, double **class_prob);
void update_beta2(int n, int K1, int p1, int ncol_x, double **x, double **y, double **class_prob, double *beta2);
void x_matrix(int n, int m, int K1, int p1, double **x, double **y, int ** indicator);
void update_beta_sd1(int n, int K1, int p1, int m, double **y, double **sd, double **beta, double * tau, double **class_prob);
void updata_tau(int n, int K, double *tau, double ** class_prob);



void update_class_prob(int n, int p, int K, int m, double **y, double **sd, int **indicator, double **beta, double * tau, double ** class_prob);
void classprob(int n, int p, int K, int m, int n1, double **y, double **sd, int **indicator, double **beta, double * tau, double **class_prob);
double density(int n, int p, int m, int n1, int K1, double **y, double **sd, int **indicator, double **beta);
double dnorm(int n1, int p1, int K1, double **y, double mean, double **sd);

void inverse(int p, double **S, double **Sinv);


#endif/* EM_H */
