EMEM = function(x, m, K, s.id, n.em, em.iter, EM.iter, nk.min, tol)
{
  n = nrow(x)
  p = ncol(x)
  nbic = numeric(n.em)
  nid = matrix(0, n.em, n)

  for(j in 1: n.em){
    #id0 = kmeans(x, K, iter.max = 1, nstart = 1, algorithm = "Forgy")$cluster
    
    id0 = initial_id0(x, K)
    em_result = tryCatch(EM_c(x, m, K, s.id, id0 = id0, em.iter, tol), error=function(e){ })
    nbic[j] = ifelse(is.null(em_result), NA, em_result$BIC)
    if(!is.null(em_result)){
      nid[j,] = em_result$id
    }else{ nid[j,] = NA}
  }
  
  obj = NULL
  for (i in 1:n.em) {
    id = nid[which(rank(nbic,na.last = TRUE, ties.method = "first")==i),]
    if(min(table(id))> nk.min){
        obj = tryCatch(EM_c(x, m, K, id0 = id, s.id, EM.iter, tol), error=function(e){ })
    }
    if(!is.null(obj)){
      break
    }
  }


  return(list(BIC = obj$BIC, llk = obj$llk, id = obj$id, Pi = obj$Pi, beta = obj$beta, s2 = obj$s2, tau = obj$tau, n_pars = obj$n_pars))

}
