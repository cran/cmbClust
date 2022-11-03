beta.name = function(p, m, orders){

   name = paste0("x", orders[1], "~1")

   if(p>1){
      for (j in 2:p){
        s = paste0("x", orders[j], "~1")
        name = c(name, s)
        if(m > 0){
            for (d in 1:(j-1)) {
                  for (i in 1:m) {
                     s = paste0("x", orders[j], "~x", orders[d], "^",i)
                     name = c(name, s)

                   }
             }

         }
      }
   }
 return(name)

}

