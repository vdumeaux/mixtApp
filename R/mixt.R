.onLoad <- function(libname, pkgname) {
  data("merged_bresat","merged_moduleColors", "goterms", "combat_data","msigdb", package=pkgname, envir=parent.env(environment()))
}