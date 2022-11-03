cmb.density = function(x, orders, m, K, beta, s2, Pi){
  
  x = as.matrix(x)
  x = as.matrix(x[, orders])
  p =ncol(x)
  n = nrow(x)
  f = numeric(n)
  mu = list(NULL)
  
  for(k in 1: K)
   {
      fk = matrix(0,n,p)
      mu[[k]] = matrix(0, n, p) 
      mu[[k]][, 1] = beta[k, 1]
      fk[,1] = dnorm(x[,1], mu[[k]][, 1], rep(s2[k, 1]^0.5,n))
      t = 1
      if(p>1){
        for (j1 in 2:p){
           xx = rep(1, n)
           if(j1 > 1 ){
              for(j2 in 1: (j1-1)){
                 for(i in 1:m){
                   xx = cbind(xx, x[, j2]^i)
                 }
              }
            }
        
           mu[[k]][, j1] = xx%*%beta[k, c((t+1):(t+ 1 + (j1-1)*m))]
           fk[,j1] = dnorm(x[,j1], mu[[k]][, j1], rep(s2[k, j1]^0.5,n))
           t =t + 1 + (j1-1)*m
        }
      }
      f = f + Pi[k]*apply(fk, 1, prod) 
  }
  

  return(list(mu = mu, f = f ))
  
}



cmb.mu = function(x, orders, m, K, beta){
  
  x = as.matrix(x)
  x = x[, orders]
  p =ncol(x)
  n = nrow(x)
  mu = list(NULL)
  
  for(k in 1: K)
  {
    mu[[k]] = matrix(0, n, p) 
    mu[[k]][, 1] = beta[k, 1]
    t = 1
    for (j1 in 2:p){
      xx = rep(1, n)
      if(j1 > 1 ){
        for(j2 in 1: (j1-1)){
          for(i in 1:m){
            xx = cbind(xx, x[, j2]^i)
          }
        }
      }
      
      mu[[k]][, j1] = xx%*%beta[k, c((t+1):(t+ 1 + (j1-1)*m))]
      t =t + 1 + (j1-1)*m
    }
  }
  return(mu)
  
}



xx_matrix = function(x, j, m){
  
    n = nrow(x)
    xx = rep(1, n)
    if(j > 1){
      for(j1 in 1: (j-1)){
        for(z in 1:m){
          xx = cbind(xx, x[, j1]^z)
        }
      }
    }

  
  return(xx)
}



cmb.plot.density = function(x, orders, pair, m, K, beta, s2, Pi){
  
  
  x = as.matrix(x)
  x = x[, orders]
  pair = sort(c(which(orders == pair[1]), which(orders == pair[2])))
  p =ncol(x)
  n = nrow(x)
  fk = matrix(0, n, p)
  f = numeric(n)
  
  for(k in 1: K)
  {
    x.new = x
    t = 1
    if(pair[1] == 1 )
      { fk[,pair[1]] = dnorm(x[,1], beta[k, 1], rep(s2[k, 1]^0.5,n))
    }else{
        x.new[, 1] = beta[k, 1]
        if(pair[1] > 2){
            for(j in 2: (pair[1]-1)){
              x.new[, j] = xx_matrix(x.new, j, m)%*%beta[k, c((t+1):(t+ 1 + (j-1)*m))]
              t =t + 1 + (j-1)*m
            }
        }
       mu = xx_matrix(x.new, pair[1], m)%*%beta[k, c((t+1):(t+ 1 + (pair[1]-1)*m))]
       fk[,pair[1]] = dnorm(x.new[, pair[1]], mu, rep(s2[k, pair[1]]^0.5,n))
       t =t + 1 + (pair[1]-1)*m
    }
    
    
    
    
    
    if((pair[2] -pair[1]) ==1)
    { 
      mu = xx_matrix(x.new, pair[2], m)%*%beta[k, c((t+1):(t+ 1 + (pair[2]-1)*m))]
      fk[,pair[2]] = dnorm(x.new[, pair[2]], mu, rep(s2[k, pair[2]]^0.5,n))
      t =t + 1 + (pair[2]-1)*m
    }else{
      for(j in (pair[1]+1):(pair[2]-1)){
         x.new[, j] = xx_matrix(x.new, j, m)%*%beta[k, c((t+1):(t+ 1 + (j-1)*m))]
         t =t + 1 + (j-1)*m
        }
      mu = xx_matrix(x.new, pair[2], m)%*%beta[k, c((t+1):(t+ 1 + (pair[2]-1)*m))]
      fk[,pair[2]] = dnorm(x.new[, pair[2]], mu, rep(s2[k, pair[2]]^0.5,n))
      t =t + 1 + (pair[2]-1)*m
    }
    
    
   
    

    f = f + Pi[k]*fk[,pair[1]]*fk[,pair[2]]
  }
  
  return(f)
}



cmb.plot.mean = function(x, orders, pair, m, k, beta){
  
  
  x = as.matrix(x)
  x = x[, orders]
  p = ncol(x)
  n = nrow(x)
 
  
  x.new = x
  t = 1
  if(pair[1] > 1 )
     { x.new[, 1] = beta[k, 1]
       if(pair[1] > 2){
        for(j in 2: (pair[1]-1)){
          x.new[, j] = xx_matrix(x.new, j, m)%*%beta[k, c((t+1):(t+ 1 + (j-1)*m))]
          t =t + 1 + (j-1)*m
        }
       }
      t =t + 1 + (pair[1]-1)*m
    }
    
  
    if((pair[2] -pair[1]) ==1)
    { 
      mu = xx_matrix(x.new, pair[2], m)%*%beta[k, c((t+1):(t+ 1 + (pair[2]-1)*m))]
    }else{
      for(j in (pair[1]+1):(pair[2]-1)){
        x.new[, j] = xx_matrix(x.new, j, m)%*%beta[k, c((t+1):(t+ 1 + (j-1)*m))]
        t =t + 1 + (j-1)*m
      }
      mu = xx_matrix(x.new, pair[2], m)%*%beta[k, c((t+1):(t+ 1 + (pair[2]-1)*m))]
    }
    

  
  return(mu)
}


