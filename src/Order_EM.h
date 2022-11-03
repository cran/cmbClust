
#ifndef Order_EM_H
#define Order_EM_H


void EM_all(int n, int p, int p_1, int K, int m, double tol, int EM_iter, double **y, int *id, double *ll);


void update_parameters1(int n, int p, int K, int m, double **y, double **s2_1, double **beta1, double *tau, double **class_prob);
void update_beta_sd2_(int n, int p_1, int K1, int m, double **y, double **s2_1, double **beta1, double *tau, double **class_prob);
double updata_sd2_(int n, int p_1, int K1, int ncol_x, double **x, double **y, double *beta_2, double *tau, double **class_prob);
void update_beta2_(int n, int p_1, int K1, int ncol_x, double **x, double **y, double **class_prob, double *beta_2);
void x_matrix1(int n, int p_1, int K1, int m, double **x, double **y);
void update_beta_sd1_(int n, int p_1, int K1, int m, double **y, double **s2_1, double **beta1, double *tau, double **class_prob);
/*void updata_tau(int n, int K, double *tau,  double ** class_prob);*/



void update_parameters2(int n, int p, int p0, int K, int m, double **y, double ***s2_2, double ***beta2, double ***mu, double *tau,  double **class_prob);
void update_beta_sd_(int n, int p, int p0, int K1, int m, double **y, double ***s2_2, double ***beta2, double ***mu, double *tau, double **class_prob);
void updata_sd_(int n, int p, int p0, int K1, int nbeta, double **x1, double **x2, double ***beta2, double ***s2_2, double ***mu, double *tau, double **class_prob);
void update_beta_(int n, int p, int p0, int K1, int nbeta, double **x1, double **x2, double ***beta2, double **class_prob);
void x_matrix2(int n, int p, int p0, int K1, int m, double **x1, double **x2, double **y);



double dmvnorm(int n1, int p1, int K1, double **y, double ***s2_2, double ***mu);
double density_(int n1, int p, int p0, int K1, int m, double **y, double **beta1, double **s2_1, double ***s2_2, double ***mu);
void classprob_(int n1, int p, int p0, int K, int m, double **y, double **beta1,  double **s2_1, double ***s2_2, double ***mu, double * tau, double ** class_prob);
void update_class_prob_(int n, int p, int p0, int K, int m, double **y, double **beta1, double **s2_1, double ***s2_2, double ***mu, double * tau, double ** class_prob);



double mixLLK_(int n, int p, int p0, int K, int m, double **y,  double **beta1, double **s2_1, double ***s2_2, double ***mu, double * tau, double ** class_prob);
double BIC(int n, int p, int p1, int K, int m, double llk);
int check_tol_(double llk, double prev_llk, double eps);


void srswor(int M, int n, int *y);
void kmeans1(int n, int p, int K, double **y, int *id);
void EMEM(int n, int p, int p2, int K, int m, double tol, int n_em, int em_iter, int nk_min, double **y, int *id, double *ll);
#endif /* Order_EM_H */
