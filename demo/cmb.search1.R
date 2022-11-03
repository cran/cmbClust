set.seed(3)
K <- 3
l <- 2
x <- as.matrix(iris[, -5])
obj <- cmb.search(x = x, l, K, method = "stepwise",
                     all.perms = FALSE, Parallel = TRUE, silent = TRUE)
# print all objects returned by the list $best.model
print(obj$best.model)	
