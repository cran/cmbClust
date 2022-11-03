#selection
select_drop = function(x, m, K, id0, s.id, BIC.all, EM.iter, nk.min, tol, silent, orders, Parallel, n.cores)
{
  p = ncol(x)
  n = nrow(x)
  BIC = s.id
  for (k in 1:K) {
    BIC[[k]][which(!BIC[[k]])] = NA
  }

  #cat("begin dropping, current BIC are ", BIC.all,"\n")

  id.true = NULL
  flag = NULL
  id.true1 = NULL
  
  if((!Parallel)||K ==1)
    {for (k in 1:K) {
      Id.true = which(s.id[[k]])
      id.true = append(id.true,list(Id.true))
      for (i in 1:length(Id.true)){
              j = Id.true[i]%/%(p^2)+1
              j1 = ifelse((Id.true[i]%%(p^2))%%p == 0, p, (Id.true[i]%%(p^2))%%p)
              j2 = ifelse((Id.true[i]%%(p^2))%%p == 0,(Id.true[i]%%(p^2))%/%p, (Id.true[i]%%(p^2))%/%p+1 )
              s.id.new = s.id
              s.id.new[[k]][j1,j2,j] = FALSE

              obj.new = EM_c(x = x, m, K, s.id =  s.id.new, id0 = id0, EM.iter, tol)
              BIC[[k]][j1,j2,j]  = obj.new$BIC
              if(min(table( obj.new$id)) < nk.min || length(table( obj.new$id))!=K)  BIC[[k]][j1,j2,j] = BIC.all+10
              else if(is.na(obj.new$llk)||is.infinite(obj.new$llk))  BIC[[k]][j1,j2,j]= BIC.all+10
              else if(min(obj.new$s2, na.rm = TRUE)< 1e-07 || any(is.na(obj.new$s2))) BIC[[k]][j1,j2,j]= BIC.all+10
                  #cat("BIC are ", BIC[[k]][j1,j2,j] ,"\n")
          }
     flag = append(flag, all(is.na(BIC[[k]])))
    }
   }else{
    
     chk = Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
     if(nzchar(chk) && chk == "TRUE"){
        n.cores = 2L
     }else{
        n.cores = ifelse(n.cores > detectCores(), detectCores(),  n.cores)
     }
    
     cl = makeCluster(n.cores)
     clusterSetRNGStream(cl, 123)
     clusterExport(cl, varlist=c("x", "m","K", "p", "s.id","nk.min","id0","BIC","BIC.all","id.true1", "EM.iter", "tol","EM_c"), envir=environment())
     result = parLapply(cl, 1:K, function(k) 
        { 
          Id.true = which(s.id[[k]])
          id.true1 = append(id.true1,list(Id.true))
          BIC1 = BIC[[k]]
          for (i in 1:length(Id.true)){
                 j = Id.true[i]%/%(p^2)+1
                 j1 = ifelse((Id.true[i]%%(p^2))%%p == 0, p, (Id.true[i]%%(p^2))%%p)
                 j2 = ifelse((Id.true[i]%%(p^2))%%p == 0,(Id.true[i]%%(p^2))%/%p, (Id.true[i]%%(p^2))%/%p+1 )
                 s.id.new = s.id
                 s.id.new[[k]][j1,j2,j] = FALSE
                 obj.new = EM_c(x = x, m, K, s.id =  s.id.new, id0 = id0, EM.iter, tol)
                 BIC1[j1,j2,j]  = obj.new$BIC
                 if(min(table( obj.new$id)) < nk.min || length(table( obj.new$id))!=K)   BIC1[j1,j2,j] = BIC.all+10
                 else if(is.na(obj.new$llk)||is.infinite(obj.new$llk))  BIC1[j1,j2,j]= BIC.all+10
                 else if(min(obj.new$s2, na.rm = TRUE)< 1e-07 || any(is.na(obj.new$s2))) BIC1[j1,j2,j]= BIC.all+10
                 #cat("BIC are ", BIC[[k]][j1,j2,j] ,"\n")
             }
          return(list(BIC1 = BIC1, id.true1 = id.true1))
        })
    stopCluster(cl)
    for (k in 1:K){
      BIC[[k]] = result[[k]]$BIC1
      id.true = append(id.true, result[[k]]$id.true1)
      flag = append(flag, all(is.na(BIC[[k]])))
    }
   }
  
  
  

  if(!all(flag)){
     id = numeric(K)
     bic = numeric(K)

     for (k in 1:K) {
         id[k] = which.min(c(BIC.all, BIC[[k]][which(!is.na(BIC[[k]]))]))
          bic[k] = min(c(BIC.all, BIC[[k]][which(!is.na(BIC[[k]]))]),  na.rm = TRUE)
     }


    if(!max(id) == 1){
        k = which.min(bic)
        BIC.all = min(bic)
        i = id.true[[k]][id[k]-1]
        j = i%/%(p^2)+1
        j1 = ifelse((i%%(p^2))%%p == 0, p, (i%%(p^2))%%p)
        j2 = ifelse((i%%(p^2))%%p == 0, (i%%(p^2))%/%p, (i%%(p^2))%/%p+1 )
        s.id[[k]][j1,j2,j] = FALSE
        obj = EM_c(x = x, m, K, s.id =  s.id, id0 = id0, EM.iter, tol)
        id0 = obj$id
        if(!silent) cat(paste0("Current BIC is ", BIC[[k]][j1,j2,j], "."), paste0("In cluster ", k, ", variable"),  paste0("x", orders[j2],"^",j), "for response", paste0("x", orders[j1]), "is excluded", "\n")
    }else{
      if(!silent) cat(paste0("Current BIC is ", BIC.all,". No variable is excluded"), "\n")
       obj = EM_c(x = x, m, K, s.id =  s.id, id0 = id0, EM.iter, tol)
       id0 = obj$id
    }
  }else{
    if(!silent) cat(paste0("Current BIC is ", BIC.all,". No variable is excluded"), "\n")
    obj = EM_c(x = x, m, K, s.id =  s.id, id0 = id0, EM.iter, tol)
    id0 = obj$id
  }

  return(list(BIC = BIC.all, llk = obj$llk, id = id0, Pi = obj$Pi, s.id = s.id, beta = obj$beta, s2 = obj$s2, tau = obj$tau, n_pars = obj$n_pars ))
}





select_add = function(x, m, K, id0, s.id, BIC.all, EM.iter, nk.min, tol, silent, orders, Parallel, n.cores){

  n = nrow(x)
  p = ncol(x)
  n = nrow(x)
  BIC = s.id

  for (k in 1:K) {
    BIC[[k]][which(!BIC[[k]])] = NA
    BIC[[k]][which(BIC[[k]])] = NA
  }


  #cat("begin adding BIC are ", BIC.all, "\n")

  s.id.1 = array(FALSE, c(p, p, m))
  for (i in 1:m) {
    s.id.1[lower.tri(s.id.1[,,i])] = TRUE
  }

  id.true = NULL
  id.true1 = NULL
  flag = NULL
  if((!Parallel)||K ==1)
     {for (k in 1:K)
         {
           Id.true = which(s.id.1)[!which(s.id.1)%in% which(s.id[[k]])]
           id.true = append(id.true,list(Id.true))
           for (i in 1:length(Id.true)){
               j = Id.true[i]%/%(p^2)+1
               j1 = ifelse((Id.true[i]%%(p^2))%%p == 0, p, (Id.true[i]%%(p^2))%%p)
               j2 = ifelse((Id.true[i]%%(p^2))%%p == 0,(Id.true[i]%%(p^2))%/%p, (Id.true[i]%%(p^2))%/%p+1 )
               s.id.new = s.id
               s.id.new[[k]][j1,j2,j] = TRUE
               obj.new = EM_c(x = x, m, K, s.id =  s.id.new, id0 = id0, EM.iter, tol)
               BIC[[k]][j1,j2,j]  = obj.new$BIC
               if(min(table( obj.new$id)) < nk.min || length(table( obj.new$id))!=K) BIC[[k]][j1,j2,j] = BIC.all+10
               else if(is.na(obj.new$llk)||is.infinite(obj.new$llk))  BIC[[k]][j1,j2,j]= BIC.all+10
               else if(min(obj.new$s2, na.rm = TRUE)< 1e-07 || any(is.na(obj.new$s2))) BIC[[k]][j1,j2,j]= BIC.all+10
               # cat(" BIC are ",   BIC[[k]][j1,j2,j], "\n")
              }
          flag = append(flag, all(is.na(BIC[[k]])))
        }
  }else{
         chk = Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
         if(nzchar(chk) && chk == "TRUE"){
              n.cores = 2L
          }else{
              n.cores = ifelse(n.cores > detectCores(), detectCores(),  n.cores)
          }
         
         cl = makeCluster(n.cores)
         clusterSetRNGStream(cl, 123)
         clusterExport(cl, varlist=c("x", "m","n","K", "p", "s.id","s.id.1","nk.min","id0","BIC","BIC.all","id.true1", "EM.iter", "tol","EM_c"), envir=environment())
         result = parLapply(cl, 1:K, function(k){
            BIC1 = BIC[[k]]
            Id.true = which(s.id.1)[!which(s.id.1)%in% which(s.id[[k]])]
            id.true1 = append(id.true1,list(Id.true))
            for (i in 1:length(Id.true)){
                    j = Id.true[i]%/%(p^2)+1
                    j1 = ifelse((Id.true[i]%%(p^2))%%p == 0, p, (Id.true[i]%%(p^2))%%p)
                    j2 = ifelse((Id.true[i]%%(p^2))%%p == 0,(Id.true[i]%%(p^2))%/%p, (Id.true[i]%%(p^2))%/%p+1 )
                    s.id.new = s.id
                    s.id.new[[k]][j1,j2,j] = TRUE
                    obj.new = EM_c(x = x, m, K, s.id =  s.id.new, id0 = id0, EM.iter, tol)
                    BIC1[j1,j2,j] = obj.new$BIC
                    if(min(table( obj.new$id)) < nk.min || length(table( obj.new$id))!=K) BIC1[j1,j2,j]= BIC.all+10
                    else if(is.na(obj.new$llk)||is.infinite(obj.new$llk))  BIC1[j1,j2,j]= BIC.all+10
                    else if(min(obj.new$s2, na.rm = TRUE)< 1e-07 || any(is.na(obj.new$s2))) BIC1[j1,j2,j]= BIC.all+10
                     # cat(" BIC are ",   BIC[[k]][j1,j2,j], "\n")
             }
          return(list(BIC1 = BIC1, id.true1 = id.true1))
        }) 
      stopCluster(cl)
      for (k in 1:K){
        BIC[[k]] = result[[k]]$BIC1
        id.true = append(id.true, result[[k]]$id.true1)
        flag = append(flag, all(is.na(BIC[[k]])))
      }
   }


  if(!all(flag)){
  id = numeric(K)
  bic = numeric(K)

  for (k in 1:K) {
    id[k] = which.min(c(BIC.all, BIC[[k]][which(!is.na(BIC[[k]]))]))
    bic[k] = min(c(BIC.all, BIC[[k]][which(!is.na(BIC[[k]]))]), na.rm = TRUE)
  }

  if(!max(id) == 1){
    k = which.min(bic)
    BIC.all = min(bic)
    i = id.true[[k]][id[k]-1]
    j = i%/%(p^2)+1
    j1 = ifelse((i%%(p^2))%%p == 0, p, (i%%(p^2))%%p)
    j2 = ifelse((i%%(p^2))%%p == 0, (i%%(p^2))%/%p, (i%%(p^2))%/%p+1 )
    s.id[[k]][j1,j2,j] = TRUE
    obj = EM_c(x = x, m, K, s.id =  s.id, id0 = id0, EM.iter, tol)
    id0 = obj$id
    if(!silent) cat(paste0("Current BIC is ", BIC[[k]][j1,j2,j], "."), paste0("In cluster ", k, ", variable"), paste0("x", orders[j2],"^",j), "for response", paste0("x", orders[j1]), "is included", "\n")
  }else{
    if(!silent) cat(paste0("Current BIC is ", BIC.all,". No variable is included"), "\n")
    obj = EM_c(x = x, m, K, s.id =  s.id, id0 = id0, EM.iter, tol)
    id0 = obj$id
  }
  }else{
    if(!silent) cat(paste0("Current BIC is ", BIC.all,". No variable is included"), "\n")
    obj = EM_c(x = x, m, K, s.id =  s.id, id0 = id0, EM.iter, tol)
    id0 = obj$id
  }

  return(list(BIC = BIC.all, llk = obj$llk, id = id0, Pi = obj$Pi, s.id = s.id, beta = obj$beta, s2 = obj$s2, tau = obj$tau, n_pars = obj$n_pars ))
}
