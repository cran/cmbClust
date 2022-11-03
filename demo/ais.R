#Experiments on the dataset ais
#Require packages "MixSim","mclust","HDclassif"
library(MixSim)
library(mclust)
library(HDclassif)



data(ais)
x <- as.matrix(ais[,c(-1,-2)])
id.ais <- as.character(ais[,1])
table(id.ais)



#=====mclust==========================================================
obj1 <- Mclust(x, G = 1:5)
summary(obj1)
table(id.ais,obj1$classification)
#calculate ARI and # parameters for the best model
(ARI.Mclust <- RandIndex(id.ais,obj1$classification)$AR)
(npr.Mclust <- obj1$df)
#calculate all ARI and parameters for the models with 1:5 clusters
ALLBIC.Mclust <- numeric(5)
ALLARI.Mclust <- numeric(5)
for (k in 1:5) {
  obj1 <- Mclust(x, G = k)
  ALLBIC.Mclust[k] <- -obj1$bic
  ALLARI.Mclust[k] <- RandIndex(id.ais, obj1$class)$AR
}
ALLBIC.Mclust
ALLARI.Mclust


#==========HDclassif============================================================

hddc.modelnames <- c("AkjBkQkDk", "AkBkQkDk", "ABkQkDk", "AkjBQkDk", 
                     "AkBQkDk", "ABQkDk", "AkjBkQkD", "AkBkQkD", 
                     "ABkQkD", "AkjBQkD", "AkBQkD", "ABQkD", "AjBQD",
                     "ABQD")
obj2 <- hddc(x, K = 1:5, model = hddc.modelnames)
table(id.ais, obj2$class)
(obj2$K)
(BIC.hddc <- -obj2$BIC) 
(ARI.hddc <- RandIndex(id.ais, obj2$class)$AR)

set.seed(111)
ALLBIC.hddc <- numeric(5)
ALLARI.hddc <- numeric(5)
for (k in 1:5) {
  obj2 <- hddc(x, K = k, model = hddc.modelnames)
  ALLBIC.hddc[k] <- -obj2$BIC
  ALLARI.hddc[k] <- RandIndex(id.ais, obj2$class)$AR
}
ALLBIC.hddc
ALLARI.hddc


# ==============cmbClust=====================
set.seed(3)
ALLBIC.cmb <- numeric(5)
ALLARI.cmb <- numeric(5)
for (k in 1:5) {
  obj3 <- cmb.search(x = x, l = 2, K = k, all.perms = FALSE, method = "stepwise",silent = TRUE,n.em = 600, nk.min = 16, max.spur = 25, tol = 1e-09, em.iter = 5, EM.iter = 500)
  ALLBIC.cmb[k] <- obj3$best.model$BIC
  ALLARI.cmb[k] <- RandIndex(id.ais, obj3$best.model$id)$AR
}
ALLBIC.cmb
ALLARI.cmb
















