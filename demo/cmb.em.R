set.seed(1)
K <- 3
l <- 2
x <- as.matrix(iris[,-5])
obj <- cmb.em(x = x, order = c(1,2,3,4), l, K, method = "stepwise", silent = TRUE, Parallel = FALSE)
print(obj)


