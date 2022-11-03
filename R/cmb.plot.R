
fnctns = function(j2, m, orders, beta){
  
  
  K = nrow(beta)
  
  st1 = (j2 + (j2-2)*(j2-1)*m/2) 
  st2 = st1 + m*(j2-1)
  b = matrix(beta[, st1:st2], K)
  
  expr = numeric()
  for(k in 1: K){
    y = paste0("X", orders[j2], " = ", round(b[k, 1],3))
    if(!all(round(b[k,],3)==0)) 
       {for(j in 1: (j2-1)){
          for(z in 1:m)
            { if(z ==1){
                  i = (z+(j-1)*m)+1
                  if(round(b[k, i],3)!=0) 
                    {if(b[k, i]<0) y = paste0(y, round(b[k, i],3), "X", orders[j]) else y = paste0(y, "+", round(b[k, i],3), "X", orders[j]) 
                    }
              }else{
                 i = (z+(j-1)*m)+1
                 if(round(b[k, i],3)!=0) 
                   {if(b[k, i]<0) y = paste0(y, round(b[k, i],3), "X", orders[j],"^", z) else y = paste0(y, "+", round(b[k, i],3), "X", orders[j],"^", z)
                   } 
            } 
        }
       }
    }
    expr[k] = y
  }
  
  
  
  
  return(expr)
}








xlabels = function(x, orders){
  
  p = length(orders)
  xlabels = numeric(p)
  for(i in 1:p){
    xlabels[i] = paste0("X", i, " - ", x[i])
  }
  
  xlabel = xlabels[1]
  if(p>1) for (i in 2:p) {xlabel = paste(xlabel, xlabels[i], sep = "\n")}
  
  return(xlabel)
}






cmb.all.plot = function(x, id.true, orders, m, K, id, beta, s2, Pi, allcolors,allpch, lwd, cex.text, cex.point, mar, oma, nlevels){
  
  
  x = as.matrix(x)
  p = ncol(x)
  n = nrow(x)
  if(is.null(allcolors))
  {allcolors = c( "dodgerblue2", "tomato4" ,"darkseagreen", "violet", "khaki", "darkorange", "skyblue1",  "violetred4", "forestgreen" ,  
                  "steelblue4", "slategrey",  "darkgoldenrod3",  "olivedrab",
                  "royalblue" , "brown","cyan2", "springgreen2")}
  
  if(is.null(allpch)){
    allpch = c(c(15,16,17,18,3,4),7:14)
    
  }
  
  opar <- par("mfrow", "oma", "mar")
  on.exit(par(opar))
  
  z = matrix(0, n, n)
  #if(p ==1) warning("for multi-dimensional data only")
  
  if(p > 1){
    par(mfrow=c(p-1, p-1), oma = oma)
    for(j2 in 2 : p)
    { 
      for(j1 in 1: (j2-1)) 
      { 
        par(mar = mar)
        x1 = seq(min(x[, orders[j1]]), max(x[, orders[j1]]), length = n)
        x2 = seq(min(x[, orders[j2]]), max(x[, orders[j2]]), length = n) 
        x.new = x
        x.new[, orders[j1]] = x1
        for(i in 1:n){
          x.new[, orders[j2]] = x2[i]
          z[, i] = cmb.plot.density(x.new, orders = orders, pair = c(orders[j1], orders[j2]), m, K, beta, s2, Pi)
        }
        if(j2 != p&&j1!=1)  {contour(x1, x2, z, nlevels =  nlevels,lwd = lwd*(0.7), col = "gray64", drawlabels = FALSE, xaxt='n', yaxt = 'n')
        } else if(j2 == p && j1!=1)   {contour(x1, x2, z, nlevels =  nlevels, lwd = lwd*(0.7), col = "gray64", drawlabels = FALSE, yaxt = 'n')
        } else if(j2 != p && j1==1)   {contour(x1, x2, z, nlevels =  nlevels, lwd = lwd*(0.7), col = "gray64", drawlabels = FALSE, xaxt='n')
        } else {contour(x1, x2, z, nlevels =  nlevels,lwd = lwd*(0.7), col = "gray64", drawlabels = FALSE)}
        
        
        if(j2 == p) mtext(text = colnames(x)[orders[j1]],side = 1,line = 2.2, cex = cex.text)
        if(j1 == 1) mtext(text = colnames(x)[orders[j2]],side = 2,line = 2.2, cex = cex.text)
        
        
        for(k in 1:K){
          x1 = seq(min(x[which(id ==k), orders[j1]]), max(x[which(id ==k), orders[j1]]), length = n) 
          x.new = x
          x.new[, orders[j1]] = x1
          mu = cmb.plot.mean(x.new, orders = orders, pair = c(j1, j2), m, k, beta)
          points(x1, mu, col = allcolors[k], type = "l", lwd = lwd)
          
          if(is.null(id.true)) 
          {points(x[which(id ==k),orders[j1]], x[which(id ==k), orders[j2]], col = allcolors[k], pch = allpch[k], cex = cex.point)
          }else{
            for(ii in 1: K){
              points(x[which(id ==k&id.true==ii),orders[j1]], x[which(id ==k&id.true==ii), orders[j2]], col = allcolors[k], pch = allpch[ii], cex = cex.point)
            }
          }
        }
        expr = fnctns(j2, m, orders, beta)
        if(j1 == j2-1) legend("topleft", expr,  col = allcolors[1:K], pch = allpch[1:K], cex = (0.8*cex.text), pt.cex = 1.2, inset=c(1,0), xpd = NA, bty="n", text.width = 2)
        
        if(j2 == 2 && p == 2)
        { 
          legend("topleft", xlabels(colnames(x), orders), cex = (0.8*cex.text), pt.cex = 1.2, inset = c(1.4,0.5), xpd = NA, bty="n",)
        }
      }
      
      if(j2 < p)
      { for(j in j2:(p-1)) 
      { par(mar = mar)
        plot(x = 1:5,y = 1:5, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
        if(j2 == 2 && j == (p-1))legend("topleft", xlabels(colnames(x), orders), cex = (0.8*cex.text), pt.cex = 1.2, inset = c(1,0),xpd = NA,bty="n",)
      }
      }
      
    }
    
  }
}







cmb.plot = function(obj, allcolors = NULL,allpch = NULL, lwd = 1,  cex.text = 1, cex.point = 0.6, mar = c(0.6,0.6,0.6,0.6), oma = c(3.5,3.5,2.5,14) , nlevels = 30) 
{
  
  id.true = NULL
  if(!is.null(obj$best.model)) obj = obj$best.model
  data = obj$data
  p = ncol(as.matrix(data))
  if(p>1) m = (ncol(obj$beta)-p)*2/(p*(p-1))
  if(p ==1) m = NULL
  opar <- par("mfrow", "oma")
  on.exit(par(opar))
  
  if(p>1){cmb.all.plot(data, id.true = id.true, orders = obj$order, m = m,  K = nrow(obj$s2), id = obj$id, beta = obj$beta, s2 = obj$s2, Pi = obj$Pi, allcolors = allcolors, allpch = allpch, lwd = lwd, cex.text = cex.text, cex.point = cex.point, mar = mar, oma = oma , nlevels = nlevels)
  }else{
    
    x.new = seq(min(data), max(data), length = 2*nrow(as.matrix(data)))
    f =  cmb.density(x.new, 1, m = m, K = nrow(obj$s2),beta = obj$beta, s2 = obj$s2, Pi = obj$Pi)$f
    par(mfrow=c(1, 1), oma = c(1,1,1,1))
    plot(x.new, f,  type = "l", main = NULL, xlab = "x", ylab = "density")
    
  }
  
}





#
#obj = cmb.em(iris[,-5], m = 2, K =3)
#pdf("C:/Users/yang wang/Desktop/plots.pdf", width = 8, height = 6, family = "Helvetica")
#cmb.plot( obj) 
#dev.off()


#obj1 = cmb.search(iris[,-5], 1, 3)
#pdf("C:/Users/yang wang/Desktop/plot_m1.pdf", width = 8, height = 6, family = "Helvetica")
#cmb.plot( obj1) 
#dev.off()


