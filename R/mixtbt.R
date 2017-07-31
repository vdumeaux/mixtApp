# Loads all datasets into memory when the package gets loaded. This is to
# speed up the Kvik Compute Service and OpenCPU calls.

.onLoad <- function(libname, pkgname) {
  utils::data("bresat", "moduleColors", "goterms", "dat", "msigdb.enrichment",
              "perm_cor_p", "mod_clinical_fdr", "net", package=pkgname,
              envir=parent.env(environment()))
}

