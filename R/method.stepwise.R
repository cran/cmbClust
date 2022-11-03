stepwise = function(x, orders, m, K, id0 = id0, n.em, em.iter, EM.iter, nk.min, tol, silent, Parallel, n.cores)
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
  
  if(is.null(id0) && K == 1){
    id0 = rep(1, n)
    obj =   EM_c(x, m, K, s.id = s.id1, id0 = id0, EM.iter, tol)
    BIC.all = obj$BIC
    id0 = obj$id
  }
  
  
  if(is.null(id0) && K > 1){
    obj = EMEM(x = x, m, K, s.id = s.id1, n.em, em.iter, EM.iter, nk.min, tol)
    id0 = obj$id
    BIC.all = obj$BIC
  }
  


  i = 0
  while(m > 0&&length(table(id0))==K){
    obj1 = select_drop(x, m, K, id0 = id0, s.id = s.id1, BIC.all, EM.iter, nk.min, tol, silent, orders = orders, Parallel, n.cores)
    s.id1.new = obj1$s.id
    id0 = obj1$id
    s.id2 = s.id1.new
    BIC.all = obj1$BIC

    for (k in 1:K){
      stop[k] = all(s.id1[[k]] == s.id1.new[[k]])
    }
    if(all(stop==1)) {
      obj = obj1
      break
    }

    
    
    obj2 = select_add(x, m, K, id0 = id0, s.id = s.id2, BIC.all, EM.iter, nk.min, tol, silent, orders = orders, Parallel, n.cores)
    s.id2.new = obj2$s.id
    id0 = obj2$id
    BIC.all = obj2$BIC
    obj = obj2

    for (k in 1:K) {
      stop[k] = all(s.id1[[k]] == s.id2.new[[k]])
    }
    if(all(stop==1)) {
      obj = obj2
       break
    }
    
    
    s.id1 =  s.id2.new

    i = i+1
    #print(i)
    

  }


   model = NULL

   
   if(!is.null(obj$s.id)) for (k in 1:K) model = append(model,list(expr( obj$s.id[[k]],p,m, orders = orders)))



  return(list(model = model, BIC = BIC.all, llk = obj$llk, id = obj$id, Pi = obj$Pi, beta = obj$beta, s2 = obj$s2, tau = obj$tau, n_pars = obj$n_pars))
}
