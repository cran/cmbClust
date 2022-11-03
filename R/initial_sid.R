inital_bac_sid = function(p,m,K)
{
  if(m == 0)
    {s.id.new = NULL
  }else{
         s.id = array(FALSE, c(p, p, m))
         for (i in 1:m) {
                s.id[lower.tri(s.id[,,i])] = TRUE
               }
         s.id.new = NULL

         for (k in 1:K) {
               s.id.new = append(s.id.new,list(s.id))
              }
  }

  return(s.id.new)
}



inital_for_sid = function(p,m,K)
{
    s.id = array(FALSE, c(p, p, m))
    s.id.new = NULL
    if(m > 0){
       for (k in 1:K) {
           s.id.new = append(s.id.new,list(s.id))
       }
    }
  return(s.id.new)
    
}
