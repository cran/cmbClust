backward = function(x, orders, m, K, id0 = id0, n.em, em.iter, EM.iter, nk.min, tol, silent, Parallel, n.cores)
{
  p = ncol(x)
  n = nrow(x)
  s.id1 = inital_bac_sid(p,m,K)
  s.id1.new = s.id1
  stop = numeric(K)
  x = x[, orders]

  #options(warn = -1)
  if(!is.null(id0)){
    obj =   EM_c(x, m, K, s.id = s.id1, id0 = id0, EM.iter, tol)
    BIC.all = obj$BIC
    id0 = obj$id
  }
  if(is.null(id0) && K > 1){
    obj = EMEM(x = x, m, K, s.id = s.id1, n.em, em.iter, EM.iter, nk.min, tol)
    id0 = obj$id
    BIC.all = obj$BIC
  }
  if(is.null(id0) && K == 1){
    id0 = rep(1, n)
    obj =   EM_c(x, m, K, s.id = s.id1, id0 = id0, EM.iter, tol)
    BIC.all = obj$BIC
    id0 = obj$id
  }

  


  i = 0
  while(m >0&&length(table(id0))==K){
    obj = select_drop(x, m, K, id0 = id0, s.id = s.id1, BIC.all, EM.iter, nk.min, tol, silent, orders = orders, Parallel, n.cores)
    s.id1.new = obj$s.id
    id0 = obj$id
    BIC.all = obj$BIC
    for (k in 1:K) {
      stop[k] = all(s.id1[[k]] == s.id1.new[[k]])
    }
    if(all(stop==1)) break
    
    
    s.id1 = s.id1.new
    #print(i)
  }

  model = NULL
  for (k in 1:K) model = append(model,list(expr( s.id1.new[[k]],p, m, orders = orders)))

  return(list(model = model, BIC = BIC.all, llk = obj$llk, id = id0, Pi = obj$Pi, beta = obj$beta, s2 = obj$s2, tau = obj$tau, n_pars = obj$n_pars))
}
