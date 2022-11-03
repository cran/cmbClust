set.seed(0)
K <- 3
l <- 2
x <- as.matrix(iris[, -5])
obj <- cmb.search(x = x, l, K, method = "stepwise", 
                     all.perms = TRUE, Parallel = TRUE, silent = TRUE)
# print all objects returned by the list $best.model
print(obj$best.model)	
#===============================================================
# print all objects returned by the list $models
#R> print(obj$models) 