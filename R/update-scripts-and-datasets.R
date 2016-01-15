update_scripts_and_datasets <- function() {
  datasets = c("merged_bresat.RData", "merged_moduleColors.RData", "topGO_merged_mod.RData", "msigdb.RData", "combat_data.RData")
  scripts = c("bresat.R", "common.R", "huc.R", "modules.R", "pathway_analyses.R")
  
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
  
  saveGOTerms() 
  
  
}

saveGOTerms <- function(){
  load("data/topGO_merged_mod.RData") 
  goterms <- all.single
  for (tissue in names(goterms)){
    for(module in names(goterms[[tissue]])){
      goterms[[tissue]][[module]]$GO.data <- NULL
      goterms[[tissue]][[module]]$resultFisher <- NULL
      goterms[[tissue]][[module]]$common <- NULL
      #goterms[[tissue]][[module]]$common <- go.common[[tissue]][[module]]
    }
  }
  save(goterms, file="data/goterms.RData") 
}

initparallel <- function(){
  library(doParallel)
  library(plyr)
  
  nodes <- detectCores() 
  cl <- makeCluster(nodes)
  registerDoParallel(cl)
}

initDev <- function(){ 
  datasets = c("merged_bresat.RData", "merged_moduleColors.RData", "topGO_merged_mod.RData", "msigdb.RData", "combat_data.RData")
  scripts = c("bresat.R", "common.R", "huc.R", "heatmap.R" , "modules.R", "pathway_analyses.R")
  datadir = "/home/mixt/data/mixt/"
  scriptdir = "/home/mixt/src/"
  
  setwd(datadir)
  for(dataset in datasets){
    load(dataset)
  }
  
  setwd(scriptdir)
  for(script in scripts){
    source(script) 
  }
  
  
}

