cmb.em = function(x, order = NULL, l, K, method = "stepwise", id0 = NULL, n.em = 200, em.iter = 5, EM.iter = 200, nk.min = NULL, max.spur = 5, tol = 1e-06, silent = FALSE, Parallel = FALSE, n.cores = 4){

  
  m = l
  
  x = as.matrix(x)
  p = ncol(x)
  if(is.null(nk.min)) nk.min = (m*(p-1)+1)*2
  if(is.null(order)) order = c(1:p)

  for (ii in 1:max.spur) {
    
      if(p >1){
          if(method == "stepwise") obj = stepwise(x, orders = order, m, K, id0 = id0, n.em, em.iter, EM.iter, nk.min, tol, silent, Parallel, n.cores)
          if(method == "forward")  obj = forward(x, orders = order, m, K,  id0 = id0, n.em, em.iter, EM.iter, nk.min, tol, silent, Parallel, n.cores)
          if(method == "backward") obj = backward(x, orders = order, m, K, id0 = id0, n.em, em.iter, EM.iter, nk.min, tol, silent, Parallel, n.cores)
          if(method == "none")     obj = none(x, orders = order, m, K, id0 = id0, n.em, em.iter, EM.iter, nk.min, tol)
        
          }
      
      if(p ==1){
        obj = none(x, orders = order, m, K, id0 = id0, n.em, em.iter, EM.iter, nk.min, tol)
      }
      

      if((!is.na(obj$llk))&&(!is.nan(obj$llk))&&(!is.infinite(obj$llk))&&(min(obj$s2, na.rm = TRUE)> 1e-07)&&(min(table(obj$id))>= nk.min)&&(length(table( obj$id))==K)){
        break}else{
          if(ii < max.spur)  warning("Spurious solution has been detected, trying a new random seed")
          if(ii == max.spur)  
            {  obj = NULL
               warning("The maximum number of attempts has been reached. No solution has been found.")
            }
        }

  }
  

  if(is.null(obj))
    {return(obj)
  }else{
     p = ncol(x)
     nbeta = (p + (p-1)*p*m/2)
     beta = matrix(obj$beta, K, nbeta, byrow = TRUE)
     invisible(colnames(beta, do.NULL = FALSE))
     colnames(beta) = beta.name(p, m, orders = order)

     s2 = matrix(obj$s2, K, p, byrow = TRUE)
     model = unlist(obj$model)


  return(list(data = x, model = model, id = obj$id, loglik = obj$llk,BIC = obj$BIC,  Pi = obj$tau, tau = obj$Pi, beta = beta, s2 = s2, order = order, n_pars = obj$n_pars))
  }
  
}






