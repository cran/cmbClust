
.onLoad <- function(libname, pkgname){
  library.dynam("cmbClust", pkgname, libname)
  set.seed(NULL)
} # End of .onLoad()

.onUnload <- function(libpath){
  library.dynam.unload("cmbClust", libpath)
} # End of .onUnload()

