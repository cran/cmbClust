OrderSch = function(x, m, K, n.em = 200, em.iter = 5, nk.min = NULL, EM.iter = 200, tol = 1e-06)
{
  
  n =nrow(x)
  p = ncol(x)
  id0 = numeric(n)
  nk.min = ifelse(is.null(nk.min), 2*(p-1)*(1+m), nk.min)
 
  #dyn.load("/Users/wangy4/Desktop/third\ project/C_mbc_ordersch/Order.dll")
  obj = .C("order", y1 = as.double(as.vector(t(x))), n1 = as.integer(n), p1 = as.integer(p), K1 = as.integer(K), m1 = as.integer(m), id = as.integer(id0), orders = c(0:(p-1)), ll = as.double(c(0,0)), n_em1 =  as.integer(n.em), em_iter1 = as.integer(em.iter), nk_min1 = as.integer(nk.min), EM_iter1 = as.integer(EM.iter),  tol1 = as.double(tol))
  #dyn.unload("/Users/wangy4/Desktop/third\ project/C_mbc_ordersch/Order.dll")
  

  return(list(BIC = obj$ll[2], llk = obj$ll[1], id = obj$id+1, orders = obj$orders+1))
  
}
















