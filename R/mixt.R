.onLoad <- function(libname, pkgname) {
  data("modules", "geneModuleOverview", "cc.biopsy_msigdb", "goterms", package=pkgname, envir=parent.env(environment()))
}