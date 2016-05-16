.onLoad <- function(libname, pkgname) {
  data("bresat","moduleColors", "goterms", "combat_data","msigdb", "TOM-net",package=pkgname, envir=parent.env(environment()))
}