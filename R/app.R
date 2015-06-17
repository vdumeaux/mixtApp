# Where data and scripts are stored 
datadir <- "/home/bjorn/mixt/data"
#scriptdir <- "/home/bjorn/mixt/src"

# Get helper scripts 
#source(paste0(scriptdir, "/bresat.R"), chdir=TRUE)
#source(paste0(scriptdir, "/heatmap.R"), chdir=TRUE)
#source(paste0(scriptdir, "/mixt-utils.r"), chdir=TRUE)

# Get datafiles 
#rawModulesFilename <- paste0(datadir, "/cc.blood-biopsy-Modules.RData")
#exprsFilename <-  paste0(datadir, "/CC-Biopsy-Expressions.RData")
#modulesFilename <-paste0(datadir, "/modules-complete-pepi.RData")

#modules <- loadModulesAndROI(rawModulesFilename,exprsFilename,modulesFilename)
  
### Set Kvik option so that the output is readable in Kvik 
#options(width=10000) 

### Where to store images
#imgpath <- "images"
#dir.create(imgpath,showWarnings = FALSE)
### Directory to store tables (output as csv files)
#tablePath <- "tables"
#dir.create(tablePath,showWarnings = FALSE)

#' Generate heatmap plot for the given tissue and module. If the heatmap
#' already exists, it finds the appropriate png file where it is supposed
#' to store a new one, it returns this file. This heat map function
#' generates both a png and a pdf. All plots are stored in the path
#' given by 'imgpath'
#' @param tissue we want to look at
#' @param module the module we want to generate a heetmap for
#' @export
#' @examples
#' heatmap("blood", "green")
#' heatmap("biopsy", "blue") 
#' 

heatmap <- function(tissue, module) { 
    plot.new()
    create.modules.heatmap(bs=mymodules[[tissue]]$bresat[[module]],exprs=mymodules[[tissue]]$exprs, 
                           clinical=mymodules[[tissue]]$heatmap.clinical, re.order=FALSE,
                           title=paste(module, tissue ,sep="-"))
    #dev.off() 
}

#' Returns a list of modules found for the given tissue
#' @param tissue tissue we want to retrieve modules for 
#' @export
getModules <- function(tissue) {
  return (names(mymodules[[tissue]]$modules))
}

#' Returns a list of all genes found in  all modules across all tissues. 
#' @export
getAllGenes <- function(){ 
  return (geneModuleOverview$gene)
}

#' Get all modules a specific gene is found in. 
#' @param gene the interesting gene
#' @export
getAllModules <- function(gene) { 
  id = match(gene,geneModuleOverview$gene)
  g = geneModuleOverview[id,]
  d = c(lapply(g,as.character))
  return(c(d$blood, d$biopsy))
}

#' Retrieves an overview of all genes and the modules they participate in.
#' @export  
getAllGenesAndModules <- function() {
  res <- NULL
  tissues <- getAllTissues()
  for (tissue in tissues){
    for(module in names(mymodules[[tissue]]$modules)) {
      if(module == "grey"){
        next 
      }
      gs <- mymodules[[tissue]]$bresat[[module]]$gene.order
      for(gene in gs){
        if(length(res[[gene]])==0) {
          res[[gene]] = list()
          for (ti in tissues){ 
            res[[gene]][[ti]] = NA
          }
        }
        res[[gene]][[tissue]] = c(module)
      }
    }
  }
  geneModuleOverview = matrix(unlist(res), nrow=length(names(res)))
  geneModuleOverview = cbind(names(res), geneModuleOverview)
  colnames(geneModuleOverview) <-  c("gene",tissues)
  geneModuleOverview = as.data.frame(geneModuleOverview) 
  save(geneModuleOverview, file="data/geneModuleOverview.RData")
  return (geneModuleOverview)
}

#' Get available tissues
#' @export
getAllTissues <- function() {
  return (names(mymodules))
}

#' Get a list of genes for a specific module and tissue.
#' @param tissue is the tissue we're interested in 
#' @param module is the module we want to get the genes from
#' @export  
getGeneList <- function(tissue,module){
  genes <- mymodules[[tissue]]$bresat[[module]]$gene.order
  up.dn <- mymodules[[tissue]]$bresat[[module]]$up.dn

  res <- matrix(c(genes,up.dn), nrow=length(genes))
  colnames(res) <- c("Gene", "up.dn")
  res = data.frame(res)
  
  # get correlation and join it 
  a = matrix(c(mymodules[[tissue]]$bresat[[module]]$up.cor, mymodules[[tissue]]$bresat[[module]]$dn.cor))
  colnames(a) <- c("cor")
  df = data.frame(a) 
  df$gene = names(a)
  df[match(res$Gene, df$gene),]
  res$cor = df$cor

  print(names(a)) 
  
  return(res)
} 

#' Get enrichment scores for specific module, tissue and specific
#' gene set if applicable. If no gene set is specified it will return
#' all gene sets
#' @param tissue is the tissue, e.g. "blood"
#' @param module is the module, e.g. "red"
#' @param genesets is a vector of genesets we want to retrieve scores from. unspecified will return all gene sets
#' @export  
getEnrichmentScores <- function(tissue, module,genesets=c()) {
  if(length(genesets) < 1 ){ 
    return (msigdb.enrichment[[tissue]][[module]]$table)
  }
  return (msigdb.enrichment[[tissue]][[module]]$table[genesets,])
}

#' Get gene set names available to the MIxT app
#' @export 
getGeneSetNames <- function() {
  return (names(msigdb.enrichment[[1]][[1]][[1]]))
}

#' Get enrichment scores for all modules 
#' @param tissue is the tissue we're interested in, e.g. blood
#' @param genesets is a vector of the genesets we want to look at, default is the first... 
#' @export 
getEnrichmentForTissue <- function(tissue, genesets=c(1)) {
  res = lapply(msigdb.enrichment[[tissue]], gs, sets = genesets)
  res = lapply(res, addTissue, tissue=tissue) 
  res = do.call("rbind", res)
  return(res)
}

gs <- function(d, sets){
  return(d$table[sets,])
}

addTissue <- function(d, tissue){
  e = d;
  e$tissue = tissue
  return(e)
}
