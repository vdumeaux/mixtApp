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
  #return(c("biopsy","blood"))
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
  if(tissue == "nblood"){
    tissue = "blood"
  }
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
  if(tissue == "nblood"){
    tissue = "blood"
  }
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
  
  p_values$common = intersections
  
  return(p_values)
}

#' Get common genes from er analysis with genelist and module
#' @export 
commonEnrichmentScoreGenes <- function(tissue, module, genelist){
  scores = userEnrichmentScores(tissue,genelist)
  return(scores$common[[module]])
}

#' Get all go term names
#' @export 
getGOTermNames <- function(){
  return(goterms$blood$blue$GO.table$Term)
}

#' Get all scores for a specific term in a tissue 
#' @export 
getGOScoresForTissue <- function(tissue,term) {
  if(tissue == "nblood"){
    tissue = "blood"
  }
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

#' Calculates eigengene correlations between tissue A and tissue B. Returns
#' the p-values. 
#' @export 
eigengeneCorrelation <- function(tissueA,tissueB){
  if(tissueA==tissueB){
    module2Cor<-NULL
    module2Pvalue<-NULL
    
    module2Cor[[tissueA]] = cor(MEs[[tissueA]], use = "p");
    module2Pvalue[[tissueA]] = corPvalueStudent(module2Cor[[tissueA]], ncol(mymodules$blood$exprs)); 
    
    res <- NULL
    
    moduleNames = names(MEs[[tissueA]])
    
    res = matrix(unlist(module2Pvalue[[tissueA]]), ncol=length(moduleNames))
    
    #log transform 
    #res = -log10(res) 
    #res[res > 10] <- 10
    
    res = cbind(moduleNames, res)
    colnames(res) = c("module", moduleNames)
    return(res) 
  }
  
  ## Correlation analyses of eigengenes across tissues
  moduleCor = cor(MEs[[tissueA]], MEs[[tissueB]], use = "p");
  modulePvalue = corPvalueStudent(moduleCor, ncol(mymodules$blood$exprs));
  
  res = matrix(unlist(modulePvalue), ncol=length(names(MEs[[tissueB]])))
  res = cbind(names(MEs[[tissueA]]),res)
  colnames(res) = c("module", names(MEs[[tissueB]]))
  return(res) 
}

getEigengenes <- function(tissue){
  return(names(MEs[[tissue]]))
}

#' Run gene overlap test between eigengenes in tissueA and tissueB. Returns 
#' p values 
#' @param tissueA is the first tissue we are interested in, e.g. "blood"
#' @param tissueB is the other, e.g. "biopsy" 
#' @export  
eigengeneHypergeometricTest <- function(tissueA, tissueB){
  moduleCor = cor(MEs[[tissueA]], MEs[[tissueB]], use = "p");
  modulePvalue = corPvalueStudent(moduleCor, ncol(mymodules$blood$exprs));
  
  hyper <- geneOverlapTest(mymodules,tissueA,tissueB)
  #hyper <- hyper[match(rownames(modulePvalue), names(MEs[[tissueA]])),
  #                                   match(colnames(modulePvalue), names(MEs[[tissueB]]) )] 
  
  ret = as.matrix(hyper, ncol=names(MEs[[tissueA]]))
  ret = cbind(rownames(hyper), ret)
  colnames(ret) = c("module", names(MEs[[tissueA]]))  
  rownames(ret) = NULL
  return (ret)
}


geneOverlapTest <- function(modules, tissueA="biopsy", tissueB="blood"){
  
  all.genes <- intersect(unlist(modules[[tissueA]]$modules[-1]), unlist(modules[[tissueB]]$modules[-1]))
  
  pvals <- sapply(modules[[tissueA]]$modules[-1], function(biopsy.mod) {
    sapply(modules[[tissueB]]$modules[-1], function(blood.mod) {
      s <- length(intersect(blood.mod, all.genes))
      e <- length(intersect(biopsy.mod, all.genes))
      com <- length(intersect(biopsy.mod, blood.mod))
      hyper.test(s, length(all.genes) - s, e, com)
    })
  })
  ret <- p.adjust(pvals, method="BH")
  dim(ret) <- dim(pvals)
  dimnames(ret) <- dimnames(pvals)
  return(ret)
}













