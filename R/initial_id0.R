initial_id0 = function(x, K){
  n = nrow(x)
  p = sample(n, K)
  id = matrix(0, n, K)
  for(k in 1:K){
    id[,k] = apply(x, 1, function(i) dist(rbind(i,x[p[k],]), method =  "euclidean"))
  }
  
  id = apply(id, 1, which.min)
  id = as.integer(id)
  return(id)
  
}
