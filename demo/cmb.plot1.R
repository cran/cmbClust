set.seed(1)
data(smltn)
K <- 2
l <- 2
order <- c(1, 2)
x <- as.matrix(smltn1[,2:3])
obj <- cmb.em(x = x, order, l, K, method = "stepwise",
                 Parallel = TRUE, silent = TRUE)
cmb.plot(obj, oma = c(3.5,3.5,2.5,15), nlevels = 25)
