.First.lib <- function(libname = .libPaths(), pkgname = "aod"){
  cat("Package aod, version", 
      read.dcf(file = system.file("DESCRIPTION", package = "aod"), fields = "Version"), "\n")
  }
