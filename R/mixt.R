.onLoad <- function(libname, pkgname) {
  utils::data("bresat", "moduleColors", "goterms", "combat_data", "msigdb", "perm_cor_p", "mod_clinical_fdr","TOM-net", package=pkgname, envir=parent.env(environment()))
  #datasets = c("merged_bresat.RData", "merged_moduleColors.RData", "topGO_merged_mod.RData", "go_common.RData", "msigdb.RData", "combat_data.RData", "merged_TOM.RData")
}
