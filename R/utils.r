load.modules <- function(dat, mod)
{
  modules <- NULL
  modules$biopsy$exprs <- dat$biopsy$matched.tumor$exprs
  modules$blood$exprs <- dat$blood$matched.tumor$exprs
  modules$biopsy$clinical <- dat$biopsy$matched.tumor$clinical.heatmap
  modules$blood$clinical <- dat$blood$matched.tumor$clinical.heatmap
  modules$biopsy$modules <- mod$biopsy$modules.GeneList
  modules$blood$modules <- mod$blood$modules.GeneList
  ## make sure the grey module always comes first
  biopsy.grey.idx <- which(names(modules$biopsy$modules) == "grey")
  blood.grey.idx <- which(names(modules$blood$modules) == "grey")
  modules$biopsy$modules <- c(modules$biopsy$modules[biopsy.grey.idx], modules$biopsy$modules[-biopsy.grey.idx])
  modules$blood$modules <- c(modules$blood$modules[blood.grey.idx], modules$blood$modules[-blood.grey.idx])
  ## remove periods and '_' from clinical names
  names(modules$biopsy$clinical) <- sub("_", " ", sub("\\.", " ", names(modules$biopsy$clinical)))
  names(modules$blood$clinical) <- sub("_", " ", sub("\\.", " ", names(modules$blood$clinical)))
  
  return(modules)
}


## Loads the modules from an rData file located at 'rawModulesFilename'
## Computes the ranksum + ROI and stores everything in a Rdata file 'modulesFilename'
loadModulesAndROI <- function(rawModulesFilename, exprsFilename, modulesFilename, tissues=c("blood", "biopsy")){ 
  if(file.exists(modulesFilename)){ 
    load(modulesFilename)
    return (modules)
  } else {
    ### Load datasets 
    load(rawModulesFilename) # modules
    load(exprsFilename)   # gene expression and others 
    
    names(cc.biopsy) <- tissues 
    names(cc.biopsy.modules) <- tissues
    modules <- load.modules(cc.biopsy, cc.biopsy.modules)
    
    # Ranksum
    modules$blood$bresat <- mclapply(modules$blood$modules[-1], function(mod) {
      sig.ranksum(modules$blood$exprs, ns=mod, full.return=TRUE)
    })
    modules$biopsy$bresat <- mclapply(modules$biopsy$modules[-1], function(mod) {
      sig.ranksum(modules$biopsy$exprs, ns=mod, full.return=TRUE)
    })
    
    ### roi function
    roi<-NULL
    
    for (tissue in tissues)
    {
      module.names <- names(modules[[tissue]]$bresat)
      roi[[tissue]]<- mclapply(module.names, function(module) {
        random.ranks(modules[[tissue]]$bresat[[module]])
      })
      names(roi[[tissue]])<-module.names
    }  
    
    for (tissue in tissues)
    {
      module.names <- names(modules[[tissue]]$bresat)
      for (module in module.names){
        modules[[tissue]]$bresat[[module]]$roi<-roi[[tissue]][[module]]
      }
    }
    
    ### define roi categories
    roi.cat<-NULL
    for (tissue in tissues)
    {
      module.names <- names(modules[[tissue]]$bresat)
      roi.cat[[tissue]]<- mclapply(module.names, function(module) {
        define.roi.regions(modules[[tissue]]$bresat[[module]], modules[[tissue]]$bresat[[module]]$roi )
      })
      names(roi.cat[[tissue]])<-module.names
    }  
    for (tissue in tissues)
    {
      module.names <- names(modules[[tissue]]$bresat)
      for (module in module.names){
        modules[[tissue]]$bresat[[module]]$roi.cat<-roi.cat[[tissue]][[module]]
      }
    }  
    save(modules, file=modulesFilename)
    return(modules)
  }
}
