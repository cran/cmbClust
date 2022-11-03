cmb.search = function(x, l, K, method = "stepwise", all.perms  = TRUE, id0 = NULL, n.em = 200, em.iter = 5, EM.iter = 200, nk.min = NULL, max.spur = 5, tol = 1e-06, silent = FALSE, Parallel = TRUE, n.cores = 4)
{
  
  x = as.matrix(x)
  p = ncol(x)
  n = nrow(x)
  if(is.null(nk.min)) nk.min = (l*p+1)*2
 

  
  if((!all.perms)&&missing(Parallel)) Parallel = FALSE
    
  if(all.perms&&(l>=1)&&(Parallel)&&p>1){
      orders = permutations(p)
      id.all = matrix(0, n, nrow(orders))
      BIC.all = numeric(nrow(orders))
      llk.all = numeric(nrow(orders))
      Pi.all = matrix(0, K, nrow(orders))
      beta.all = NULL
      s2.all = NULL
      tau.all = NULL
      model.all = NULL
      n_pars.all = numeric(nrow(orders))
      
      chk = Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if(nzchar(chk) && chk == "TRUE"){
               n.cores = 2L
      }else{
               n.cores = ifelse(n.cores > detectCores(), detectCores(),  n.cores)
               n.cores = ifelse(n.cores > nrow(orders), nrow(orders),  n.cores)
      }
  
  
      cl = makeCluster(n.cores)
      clusterSetRNGStream(cl, 123)
      clusterExport(cl, varlist=c("x", "orders", "l","K", "method", "id0", "n.em", "em.iter", "EM.iter", "nk.min", "max.spur","tol","silent","cmb.em","none","stepwise","forward", "backward", "initial_id0","inital_bac_sid", "inital_for_sid","select_drop","select_add", "EMEM","EM_c","expr", "beta.name"), envir=environment())
      result = parLapply(cl, 1:nrow(orders), function(i) {cmb.em(x = x,  order = orders[i,], l, K, method, id0 = id0, n.em, em.iter, EM.iter, nk.min, max.spur, tol, silent)})
      stopCluster(cl)



      for (i in 1:nrow(orders)){
           if(!is.null(result[[i]]))
            {   id.all[, i] = result[[i]]$id
                BIC.all[i] = result[[i]]$BIC
                llk.all[i] = result[[i]]$loglik
                Pi.all[,i]  = result[[i]]$Pi
                n_pars.all[i] = result[[i]]$n_pars
                beta.all = append(beta.all,list(result[[i]]$beta))
                s2.all = append(s2.all,list(result[[i]]$s2))
                tau.all= append(tau.all,list(result[[i]]$tau))
                model.all = append(model.all,list(result[[i]]$model))
           }else{
                id.all[, i] = NA
                BIC.all[i] = NA
                llk.all[i] = NA
                Pi.all[,i]  = NA
                n_pars.all[i] = NA
                beta.all = append(beta.all, NA)
                s2.all = append(s2.all, NA)
                tau.all = append(tau.all, NA)
                model.all = append(model.all, NA)
          }
      }

      if(!all(is.na(BIC.all)))
        {for(i in 1: nrow(orders)){
             obj1 = result[[which(rank(BIC.all, ties.method = "random")==i)]]
             if(min(table(obj1$id))> nk.min) break
        }
      }else{
         obj1 = NA
         warning("The maximum number of attempts has been reached. No solution has been found.")
        }


      obj2 = list(model = model.all, id = id.all, loglik = llk.all, BIC = BIC.all,  Pi = Pi.all, tau = tau.all, beta = beta.all, s2 = s2.all, order = orders, n_pars = n_pars.all)
      return(list(best.model = obj1, models = obj2))
  }
  
  
  if(all.perms&&(l>0)&&(!Parallel)&&p>1){
     orders = permutations(p)
     id.all = matrix(0, n, nrow(orders))
     BIC.all = numeric(nrow(orders))
     llk.all = numeric(nrow(orders))
     Pi.all = matrix(0, K, nrow(orders))
     beta.all = NULL
     s2.all = NULL
     tau.all = NULL
     model.all = NULL
     n_pars.all = numeric(nrow(orders))
     
     
     for (i in 1:nrow(orders)) {
       if(!silent){
         epr = paste("================== Conditional order is", paste0("x",orders[i, 1],","))
         for(j1 in 2:p){
           epr = paste0(epr, " ",paste0("x",orders[i, j1]), "|")
           for(j2 in 1:(j1-1)){
            epr = paste0(epr, paste0("x",orders[i, j2]), ",")
           }
         }
         cat(epr, "==================","\n")
       }
       obj = cmb.em(x = x, order = orders[i, ], l, K, method, id0, n.em, em.iter, EM.iter, nk.min, max.spur, tol, silent)
       id.all[, i] =  obj$id
       BIC.all[i] = obj$BIC
       llk.all[i] = obj$loglik
       Pi.all[,i]  = obj$Pi
       n_pars.all[i] = obj$n_pars
       beta.all = append(beta.all,list(obj$beta))
       s2.all = append(s2.all,list(obj$s2))
       tau.all= append(tau.all,list(obj$tau))
       model.all = append(model.all,list(obj$model))
     }
    
     for(i in 1: nrow(orders)){
       id_best = which(rank(BIC.all, ties.method = "random")==i)
       if(min(table(id.all[, id_best]))> nk.min) break
     }
    
    obj1 = list(data = x, model = model.all[[id_best]], id = id.all[,id_best],loglik = llk.all[id_best], BIC = BIC.all[id_best],  Pi = Pi.all[,id_best], tau = tau.all[[id_best]], beta = beta.all[[id_best]], s2 = s2.all[[id_best]], order = orders[id_best, ], n_pars = n_pars.all[id_best])
    obj2 = list(model = model.all, id = id.all, loglik = llk.all, BIC = BIC.all,  Pi = Pi.all, tau = tau.all, beta = beta.all, s2 = s2.all, order = orders, n_pars = n_pars.all)
    return(list(best.model = obj1, models = obj2))
    
  }
  
  
 
  if(all.perms&&l == 0&&p>1){
    
    if(!silent){
      epr = paste("================== Conditional order is", paste0("x1",","))
      for(j1 in 2:p){
        epr = paste0(epr, " ",paste0("x",j1), "|")
        for(j2 in 1:(j1-1)){
          epr = paste0(epr, paste0("x",j2), ",")
        }
      }
      cat(epr, "==================","\n")
    }
    obj = cmb.em(x = x, order = c(1:p), l, K, method, id0, n.em, em.iter, EM.iter, nk.min, max.spur, tol, silent)
    return(list(best.model = obj))
  }
  
  
  if(p==1){
    if(!silent){
      epr = paste("================== Conditional order is", "x1~`1")
      cat(epr, "==================","\n")
    }
    obj = cmb.em(x = x, order = c(1:p), l, K, method, id0, n.em, em.iter, EM.iter, nk.min, max.spur, tol, silent)
    return(list(best.model = obj))
  }
  
  
 
  
  
  if(!all.perms&&p>1){
    obj =  OrderSch(x = x, l, K, n.em, em.iter, nk.min, EM.iter, tol)
    orders = obj$orders
    
    
    if(!silent){
      epr = paste("================== Conditional order is", paste0("x",orders[1],","))
      for(j1 in 2:p){
        epr = paste0(epr, " ",paste0("x",orders[j1]), "|")
        for(j2 in 1:(j1-1)){
          epr = paste0(epr, paste0("x",orders[j2]), ",")
        }
      }
      cat(epr, "==================","\n")
    }
  
    obj1 = cmb.em(x = x, order = orders, l, K, method, id0, n.em, em.iter, EM.iter, nk.min,  max.spur, tol, silent, Parallel = FALSE, n.cores)
    return(list(best.model = obj1))
  }

  
  
  
}
