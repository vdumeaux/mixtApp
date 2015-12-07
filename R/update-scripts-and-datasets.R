update_scripts_and_datasets <- function() {
  datasets = c("merged_bresat.RData", "topGO_merged_mod.RData")
  scripts = c("bresat.R", "common.R", "huc.R", "modules.R", "pathway_analyses")
  
  datadir = "/home/mixt/data/mixt/"
  scriptdir = "/home/mixt/src/"
  
  setwd("/home/rstudio/mixt/")
  
  for (script in scripts) {
    filename = paste0("R/", script)
    file.copy(paste0(scriptdir, script), filename, overwrite = TRUE, copy.mode = TRUE)
    Sys.chmod(filename, "0666") 
  }
  
  for (dataset in datasets) {
    filename = paste0("data/", dataset)
    file.copy(paste0(datadir, dataset), filename, overwrite = TRUE, copy.mode = TRUE)
    Sys.chmod(filename, "0666") 
  }
}