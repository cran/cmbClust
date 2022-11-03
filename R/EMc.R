EM_c = function(x, m, K, s.id, id0, iter, tol)
{

  n =nrow(x)
  p = ncol(x)

  nbeta = (p + (p-1)*p*m/2)
  indicator = rep(1, K*nbeta)
  if(m >0&&p>1){
        ind = NULL
        for(k in 1:K){
              ind = c(ind, 1)
              for(j1 in 2:p){
                 ind = c(ind, 1)
                     for(j2 in 1:(j1-1)){
                       for(z in 1:m){
                        ind = c(ind, s.id[[k]][j1,j2,z])
                     }
                 }
              }
       }

       indicator = ind
  }

  beta = indicator


  Pi = matrix(0, n, K)
  for (k in 1:K){Pi[which(id0 == k),k] = 1}
  Pi = c(t(Pi))

  id0 = id0 -1



  #dyn.load("C:/Users/yang wang/Desktop/third project/C_mbc_EM/mbc.dll")
  obj = .C("mbc", y1 = as.double(as.vector(t(x))), n1 = as.integer(n), K1 = as.integer(K), p1 = as.integer(p), m1 = as.integer(m), id = as.integer(id0), ll = as.double(c(0,0,0)), tau= as.double(rep(1/K,K)), indicator1 = as.integer(indicator), beta1 = as.double(beta), sd1 = as.double(rep(1,K*p)), class_prob1 = as.double(Pi), niter1 = as.integer(iter), tol1 = as.double(tol))
  #dyn.unload("C:/Users/yang wang/Desktop/third project/C_mbc_EM/mbc.dll")
  

  Pi = matrix(obj$class_prob1,n,k, byrow = "TRUE")

  return(list(BIC = obj$ll[2], llk = obj$ll[1], id = obj$id+1, Pi = Pi,s.id = s.id, beta = obj$beta1, s2 = obj$sd1, tau = obj$tau, n_pars = obj$ll[3]))

}




