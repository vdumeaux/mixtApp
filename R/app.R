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

heatmap <- function(tissue,module) { 
    plot.new()
    create.modules.heatmap(modules[[tissue]]$bresat[[module]],modules[[tissue]]$clinical,
                           title=Hmisc::capitalize(paste(tissue, module)))
    #dev.off() 
}

#' Returns a list of modules found for the given tissue
#' @param tissue tissue we want to retrieve modules for 
#' @export
getModules <- function(tissue) {
  return (names(modules[[tissue]]$modules))
}

#' Returns a list of all genes found in  all modules across all tissues. 
#' @export
getAllGenes <- function(){ 
  genesAndModules = getAllGenesAndModules()
  return (genesAndModules$gene)
}

#' Get all modules a specific gene is found in. 
#' @param gene the interesting gene
#' @export
getAllModules <- function(gene) { 
  genesAndModules = getAllGenesAndModules() 
  id = match(gene,genesAndModules$gene)
  g = genesAndModules[id,]
  d = c(lapply(g,as.character))
  return(c(d$blood, d$biopsy))
}

#' Retrieves all genes and the modules they participate in.
#' @export  
getAllGenesAndModules <- function() {
  res <- NULL
  tissues <- c("blood", "biopsy")
  for (tissue in tissues){
    for(module in names(modules[[tissue]]$modules)) {
      if(module == "grey"){
        next 
      }
      gs <- modules[[tissue]]$bresat[[module]]$gene.order
      for(gene in gs){
        if(length(res[[gene]])==0) {
          res[[gene]] = list()
          res[[gene]][["blood"]] = NA
          res[[gene]][["biopsy"]] = NA
        }
        res[[gene]][[tissue]] = c(module)
      }
    }
  }
  genes = matrix(unlist(res), nrow=length(names(res)))
  genes = cbind(names(res), genes)
  colnames(genes) <-  c("gene",tissues)
  genes = as.data.frame(genes) 
  return (genes)
}

#' Get available tissues
#' @export
getTissues <- function() {
  return (names(modules))
}

#' Get a list of genes for a specific module and tissue.
#' @param tissue is the tissue we're interested in 
#' @param module is the module we want to get the genes from
#' @export  
getGeneList <- function(tissue,module){
  genes <- modules[[tissue]]$bresat[[module]]$gene.order
  up.dn <- modules[[tissue]]$bresat[[module]]$up.dn
  res <- matrix(c(genes,up.dn), nrow=length(genes))
  colnames(res) <- c("Gene", "up.dn")
  return(res)
} 