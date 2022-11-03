expr = function(s.id, p, m, orders)
{
  y = NULL

  for (j in 1:p){
    y_f = paste0("x",orders[j])
    y_f = paste(y_f, "~ 1")
    if(m > 0){
      for (d in 1:p) {
        for (i in 1:m) {
          if(s.id[j,d,i]) {
            x = paste0("x",orders[d])
            x = paste0(x,"^",i)
            y_f = paste(y_f, "+",x)
          }
        }
      }
    }
    if(!is.null(y)) {y = paste(y,",",y_f)} else{
      y = y_f
    }
  }

  return(y)
}


