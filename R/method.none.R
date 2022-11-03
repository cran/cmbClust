none = function(x, orders, m, K, id0, n.em, em.iter, EM.iter, nk.min, tol){

  p = ncol(x)
  n = nrow(x)
  s.id1 = inital_bac_sid(p,m,K)
  x = as.matrix(x[, orders])
  #options(warn = -1)
  
  if(!is.null(id0)){
    obj =   EM_c(x, m, K, s.id = s.id1, id0 = id0, EM.iter, tol)
  }
  
  
 
  
  if(is.null(id0) && K == 1){
    id0 = rep(1, n)
    obj =   EM_c(x, m, K, s.id = s.id1, id0 = id0, EM.iter, tol)
  }

  if(is.null(id0) && K > 1){
    obj = EMEM(x = x, m, K, s.id = s.id1, n.em, em.iter, EM.iter, nk.min, tol)
    
  }


  
  model = NULL
  for (k in 1:K) model = append(model,list(expr( s.id1[[k]],p,m, orders = orders)))


  return(list(model = model, BIC = obj$BIC, llk = obj$llk, id = obj$id, Pi = obj$Pi, beta = obj$beta, s2 = obj$s2, tau = obj$tau, n_pars = obj$n_pars))
}

