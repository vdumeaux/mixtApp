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
  
  # get correlation and merge  
  a = matrix(unlist(c(mymodules[[tissue]]$bresat[[module]]$up.cor, mymodules[[tissue]]$bresat[[module]]$dn.cor)))
  colnames(a) <- c("cor")
  df = data.frame(a) 
  df$Gene = names(c(mymodules[[tissue]]$bresat[[module]]$up.cor, mymodules[[tissue]]$bresat[[module]]$dn.cor))
  
  res = merge(res, df, by="Gene", sort=FALSE) 
  return(res)
} 

#' Get enrichment scores for specific module, tissue and specific
#' gene set if applicable. If no gene set is specified it will return
#' all gene sets. Only returns genesets with updn pvalue < 1. Since
#' p-value is stored as a factor we're checking for != 1. 
#' @param tissue is the tissue, e.g. "blood"
#' @param module is the module, e.g. "red"
#' @param genesets is a vector of genesets we want to retrieve scores from. unspecified will return all gene sets
#' @export  
getEnrichmentScores <- function(tissue, module,genesets=c()) {
  if(length(genesets) < 1 ){ 
    return (subset(msigdb.enrichment[[tissue]][[module]]$results, updn.pval != 1))
  }
  return (subset(msigdb.enrichment[[tissue]][[module]]$results[genesets,], updn.pval != 1))
}

#' Get gene set names available to the MIxT app
#' @export 
getGeneSetNames <- function() {
  return (names(msigdb.enrichment[[1]][[1]][[1]]$sig.set))
}

#' Get enrichment scores for all modules 
#' @param tissue is the tissue we're interested in, e.g. blood
#' @param genesets is a vector of the genesets we want to look at, default is the first... 
#' @export 
getEnrichmentForTissue <- function(tissue, genesets=c(1)) {
  res = lapply(msigdb.enrichment[[tissue]], gs, sets = genesets)
  res = lapply(res, addTissue, tissue=tissue) 
  res = do.call("rbind", res)
  res = subset(res, updn.pval != 1)
  return(res)
}

gs <- function(d, sets){
  return(d$results[sets,])
}

addTissue <- function(d, tissue){
  e = d;
  e$tissue = tissue
  return(e)
}

#' Get go terms for the given module and tissue. Possible to specify
#' which specific terms you're interested in. 
#' @export 
#' @param tissue is the tissue, e.g. blood
#' @param module is the module, e.g. blue, pink etc.
#' @param terms are a vector GO terms we're interested in, default is all, given as a vector. 
getGOTerms <- function(tissue, module, terms=c()){
  if(length(terms) < 1){
    return (subset(goterms[[tissue]][[module]]$GO.table, classicFisher != "1.00000"))
  }
  return (subset(subset(goterms[[tissue]][[module]]$GO.table, Term==terms), classicFisher != "1.00000"))
}


saveGOTerms <- function(){
  load("data/topGO.RData.latest") 
  goterms <- all.single
  for (tissue in names(goterms)){
    for(module in names(goterms[[tissue]])){
      goterms[[tissue]][[module]]$GO.data <- NULL
      goterms[[tissue]][[module]]$resultFisher <- NULL
      goterms[[tissue]][[module]]$common <- NULL
      goterms[[tissue]][[module]]$common <- go.common[[tissue]][[module]]
    }
  }
  save(goterms, file="data/goterms.RData") 
}

#' Get common genes for given tissue, module, geneset and optionally status
#' @param geneset is the geneset, from any of the msgidb sets
#' @param status is the status, up.common, dn.common or default updn.common
#' @export 
getCommonGenes <- function(tissue, module, geneset, status = "updn.common") { 
    return (msigdb.enrichment[[tissue]][[module]][[status]][[geneset]])
} 

#' Get the common genes from the go term anlysis for a specific
#' tissue, module and GO term id. 
#' @export
#' @param gotermID is the GO term id, e.g. "GO:0070848" 
getCommonGOTermGenes <- function(tissue,module,gotermID){
   return (goterms[[tissue]][[module]]$common[[gotermID]])
}

#' Get p-values from the hypergeometric test between an arbirtrary gene list
#' and the modules. User must specify tissue and genelist
#' @param tissue is the tissue we're interested in, e.g. "blood"
#' @param genelist is the genelist as a vector. 
#' @export
userEnrichmentScores <- function(tissue, genelist) {
  modules = mymodules[[tissue]]$modules[-1]
  all_genes = unlist(modules)
  genelist = genelist[genelist %in% all_genes]
  intersections = lapply(modules, function(module) intersect(module, genelist)) 
  
  q = sapply(intersections, length) - 1
  k = sapply(modules, length)
  m = length(genelist)
  n = length(all_genes) - length(genelist) 
  
  p_values = p.adjust(phyper(q, m, n, k, lower.tail=FALSE), method="BH") 
  p_values = as.data.frame(p_values) 
  names(p_values) <- c("p-value")
  p_values$module = row.names(p_values) 
  row.names(p_values) <- NULL
  
  return(p_values)
}

#' Get all go term names
#' @export 
getGOTermNames <- function(){
  return(goterms$blood$blue$GO.table$Term)
}

#' Get all scores for a specific term in a tissue 
#' @export 
getGOScoresForTissue <- function(tissue,term) {
  res <- NULL
  for (i in 1:length(goterms[[tissue]])){
    module = names(goterms[[tissue]])[i]
    score = subset(goterms[[tissue]][[module]]$GO.table, Term == term)
    score$module = module
    res[[module]] <- score
  }
  res = do.call(rbind, res) 
  return(res)
}
