permutations = function(n){

  x = numeric(n*factorial(n))
  
 
  obj = .C("AllPerms", size1 = as.integer(n), perms1 = as.integer(x))
 
  
  x = matrix(obj$perms1+1, factorial(n),n, byrow = TRUE)
  
  return(x)

}



  