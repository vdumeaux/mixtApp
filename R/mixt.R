.onLoad <- function(libname, pkgname) {
  data("modules", "geneModuleOverview", "cc.biopsy_msigdb", package=pkgname, envir=parent.env(environment()))
}