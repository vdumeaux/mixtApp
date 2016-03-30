#' Generate heatmap plot for the given tissue and module. If the heatmap
#' already exists, it finds the appropriate png file where it is supposed
#' to store a new one, it returns this file. This heat map function
#' generates both a png and a pdf. All plots are stored in the path
#' given by 'imgpath'
#' @param tissue we want to look at
#' @param module the module we want to generate a heetmap for
#' @param orderByModule which module we want to order patients by. Default is order by module. 
#' @param orderByTissue the tissue where the orderByModule module is found. Default is the same tissue. 
#' @param re.order set to true if patients should be re ordered. Defaults to FALSE. 
#' @import colorspace
#' @export
#' @examples
#' heatmap("blood", "green")
#' heatmap("biopsy", "blue") 
#' 

heatmap <- function(tissue, module, re.order=FALSE, orderByModule=NULL, orderByTissue=NULL, subtypes=NULL) { 
    plot.new()
    title = "" 
    if(is.null(orderByModule)){
      orderByModule = module
    }
    if(is.null(orderByTissue)){
      orderByTissue = tissue
    }
    
    if(tissue == "blood"){
      clinical=huc.color.clinical(dat$blood$clinical)[,c(1,3:15, 17, 19,21:25, 27, 29:33)]
    }
    if(tissue =="biopsy"){
      clinical=huc.color.clinical(dat$biopsy$clinical)[,c(1,3:15,17, 21:25, 27, 29:33)]
    }
    
    if(tissue == "nblood") {
      clinical=huc.color.clinical(dat$nblood$clinical)[, c(19,22:27)]
    }
    
    if(tissue == "bnblood") {
      clinical=huc.color.clinical(dat$nblood$clinical)[, c(1,3:15,17, 19,21:25,26, 27, 29:33)]
    }
          
    if(re.order == FALSE || (orderByTissue==tissue && orderByModule==module)){
      title = paste(module," module from ",tissue, sep="") 
    }  else {
      title = paste(paste(module, " module from ",tissue ,sep=""), "ordered by", orderByModule, "module from", orderByTissue, sep=" ") 
    }
    

    
    create.modules.heatmap(bs = bresat[[tissue]][[module]], 
                           exprs = dat[[tissue]]$exprs, 
                           clinical = clinical,
                           re.order = re.order,
                           bs.order.by = bresat[[orderByTissue]][[orderByModule]],
                           title = title)
    #dev.off() 
}


#' Returns a list of modules found for the given tissue
#' @param tissue tissue we want to retrieve modules for 
#' @return vector of module names 
#' @export
getModules <- function(tissue) {
  return (names(bresat[[tissue]]))
}

#' Returns a list of all genes found in  all modules across all tissues. 
#' @return vector of gene names 
#' @export
getAllGenes <- function(){ 
  return(names(c(moduleColors$biopsy, moduleColors$blood)))
}

#' Get all modules a specific gene is found in. 
#' @param gene the interesting gene
#' @return vector with module name from blood and module name from biopsy
#' @export
getAllModules <- function(gene) { 
  id = match(gene, names(moduleColors$blood)) # same id blood or biopsy 
  bloodModule = as.character(moduleColors$blood[id])
  biopsyModule = as.character(moduleColors$biopsy[id])
  return(c(bloodModule, biopsyModule))
}

#' Retrieves an overview of all genes and the modules they participate in.
#' @export  
getAllGenesAndModules <- function() {
  res <- NULL
  tissues <- getAllTissues()
  for (tissue in tissues){
    for(module in names(bresat[[tissue]])) {
      if(module == "grey"){
        next 
      }
      gs <- bresat[[tissue]][[module]]$gene.order
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
  return (names(bresat))
}

#' Get a list of genes for a specific module and tissue.
#' @param tissue is the tissue we're interested in 
#' @param module is the module we want to get the genes from
#' @return matrix with columns: gene name, up.dn, cor 
#' @export  
getGeneList <- function(tissue,module){
  genes <- bresat[[tissue]][[module]]$gene.order
  up.dn <- bresat[[tissue]][[module]]$up.dn

  res <- matrix(c(genes,up.dn), nrow=length(genes))
  colnames(res) <- c("Gene", "up.dn")
  res = data.frame(res)
  
  # get correlation and merge  
  a = matrix(unlist(c(bresat[[tissue]][[module]]$up.cor, bresat[[tissue]][[module]]$dn.cor)))
  colnames(a) <- c("cor")
  df = data.frame(a) 
  df$Gene = names(c(bresat[[tissue]][[module]]$up.cor, bresat[[tissue]][[module]]$dn.cor))
  
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
  res = NULL
  if(length(genesets) < 1 ){ 
    res = subset(msigdb.enrichment[[tissue]][[module]]$results, updn.pval != 1)
  } else { 
    res = subset(msigdb.enrichment[[tissue]][[module]]$results[genesets,], updn.pval != 1)
  }
  # force ordering on p-val
  res$updn.pval = as.numeric(as.character(res$updn.pval))
  res = res[with(res, order(updn.pval)), ]
  return(res)
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
  
  if(tissue == "nblood" || tissue == "bnblood"){
    tissue = "blood" 
  }
  
  res = lapply(msigdb.enrichment[[tissue]], gs, sets = genesets)
  res = lapply(res, addTissue, tissue=tissue) 
  res = do.call("rbind", res)
  res = subset(res, updn.pval != 1)
  return(res)
}

#' @export
gs <- function(d, sets){
  return(d$results[sets,])
}

#' @export 
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




#' Get common genes for given tissue, module, geneset and optionally status
#' @param geneset is the geneset, from any of the msgidb sets
#' @param status is the status, up.common, dn.common or default updn.common
#' @export 
getCommonGenes <- function(tissue, module, geneset, status = "updn.common") { 
    return (msigdb.enrichment[[tissue]][[module]][[status]][[geneset]])
} 

#' Get the common genes from the go term anlysis for a specific
#' tissue, module and GO term id. 
#' @param gotermID is the GO term id, e.g. "GO:0070848" 
#' @export
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
  modules = names(bresat[[tissue]])
  all_genes = names(moduleColors[[tissue]])
  genelist = genelist[genelist %in% all_genes]
  intersections = lapply(modules, function(module) {
    intersect(bresat[[tissue]][[module]]$gene.order, genelist)
  }) 
  
  names(intersections) <- modules 
  
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
getGOScoresForTissue <- function(tissue, term) {
  if(tissue == "nblood"){
    tissue = "blood"
  }
  res <- NULL
  for (i in 1:length(goterms[[tissue]])){
    module = names(goterms[[tissue]])[i]
    score = subset(goterms[[tissue]][[module]]$GO.table, Term == term)
    
    if(dim(score)[1] < 1){
      next
    }
    
    score$module = module
    res[[module]] <- score
  }
  if(!is.null(res)){ 
    res = do.call(rbind, res) 
  } else { 
    res = 0
  }
    return(res)
}

#' Calculates eigengene correlations between tissue A and tissue B. Returns
#' the p-values. 
#' @export 
eigengeneCorrelation <- function(tissueA,tissueB){
  MEs = computeEigengenes(c(tissueA,tissueB)) 
  
  if(tissueA==tissueB){
    module2Cor<-NULL
    module2Pvalue<-NULL
    
    module2Cor[[tissueA]] = cor(MEs[[tissueA]], use = "p");
    module2Pvalue[[tissueA]] = WGCNA::corPvalueStudent(module2Cor[[tissueA]], ncol(dat$blood$exprs)); 
    
    res <- NULL
    
    moduleNames = names(bresat[[tissueA]])
    
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
  modulePvalue = WGCNA::corPvalueStudent(moduleCor, ncol(dat$blood$exprs));
  res <- NULL
  res = matrix(unlist(modulePvalue), ncol=length(names(MEs[[tissueB]])))
  rownames(res) = NULL
  res = cbind(names(bresat[[tissueA]]),res)
  colnames(res) = c("module", names(bresat[[tissueB]]))
  return(res) 
}

#' Computes module eigengenes for modules in the given tissues
#' @param tissues is a vector of tissues 
#' @return MEs object with module eigengenes for each tissue
computeEigengenes <- function(tissues) { 
  MEs = NULL 
  for (tissue in tissues)
  {
    # Recalculate MEs with color labels
    MEs[[tissue]] = WGCNA::moduleEigengenes(t(dat[[tissue]]$exprs), moduleColors[[tissue]])$eigengenes
    MEs[[tissue]] = WGCNA::orderMEs(MEs[[tissue]])
    
    #Exclude grey module
    MEs[[tissue]]<-MEs[[tissue]][-which(names(MEs[[tissue]])=="MEgrey")]
  }
  return(MEs)
}

#' Run gene overlap test between eigengenes in tissueA and tissueB. Returns 
#' p values 
#' @param tissueA is the first tissue we are interested in, e.g. "blood"
#' @param tissueB is the other, e.g. "biopsy" 
#' @export  
moduleHypergeometricTest <- function(tissueA, tissueB){
  
  MEs = computeEigengenes(c(tissueA,tissueB)) 
  
  moduleCor = cor(MEs[[tissueA]], MEs[[tissueB]], use = "p");
  modulePvalue = WGCNA::corPvalueStudent(moduleCor, ncol(dat[[tissueA]]$exprs));
  
  hyper <- geneOverlapTest(mymodules,tissueA,tissueB)
  hyper <- hyper[match(rownames(modulePvalue), paste("ME", rownames(hyper), sep="")),
                 match(colnames(modulePvalue), paste("ME", colnames(hyper), sep=""))]
  
  #hyper = t(hyper) 
  cols = colnames(hyper) 
  ret = NULL
  ret = as.matrix(hyper, ncol=length(colnames(modulePvalue)))
  ret = cbind(rownames(hyper), ret)
  colnames(ret) = c("module", cols)  
  
  rownames(ret) = NULL
  return (ret)
}

#' Compute gene overlap between genes from two tissues 
#'@export 
geneOverlapTest <- function(modules, tissueA="blood", tissueB="biopsy"){
  all.genes <- names(moduleColors$blood)
  
  pvals <- sapply(bresat[[tissueB]], function(tissueB.mod) {
    sapply(bresat[[tissueA]], function(tissueA.mod) {
      
      white <- length(intersect(tissueB.mod$gene.order, all.genes))
      black <- length(all.genes) - white
      total.drawn <- length(intersect(tissueA.mod$gene.order, all.genes))
      white.drawn <- length(intersect(tissueA.mod$gene.order, tissueB.mod$gene.order))
      
      phyper(white.drawn - 1, white, black, total.drawn, log.p = FALSE, lower.tail=FALSE)
    })
  })
  
  
  ret <- p.adjust(pvals, method="BH")
  dim(ret) <- dim(pvals)
  dimnames(ret) <- dimnames(pvals)
  return(ret)
}

#' ROI fisher's exact test
#' @export
roiTest <- function(tissueA="blood", tissueB="biopsy"){
  
  MEs = computeEigengenes(c(tissueA,tissueB)) 
  
  moduleCor = cor(MEs[[tissueA]], MEs[[tissueB]], use = "p");
  modulePvalue = WGCNA::corPvalueStudent(moduleCor, ncol(dat[[tissueA]]$exprs));
  
  # Define roi categories (from ROI.R) 
  roi.cat<-NULL
  for (tissue in c(tissueA, tissueB))
  {
    module.names <- names(bresat[[tissue]])
    roi.cat[[tissue]]<- parallel::mclapply(module.names, function(module) {
      define.roi.regions(bresat[[tissue]][[module]], bresat[[tissue]][[module]]$roi)
    })
    names(roi.cat[[tissue]])<-module.names
  }  
  
  # Fisher's exact teset between roi categories
  mod.roi <- NULL
  mod.roi <- plyr::laply(roi.cat[[tissueA]], function(x){
    plyr::laply(roi.cat[[tissueB]], function(y){
  
      fisher.test(x, y, workspace = 2e+07, hybrid=TRUE)$p
    })
  })
  
  rownames(mod.roi) <- names(roi.cat[[tissueA]])
  colnames(mod.roi) <- names(roi.cat[[tissueB]]) 
  
  mod.roi <- mod.roi[match(rownames(modulePvalue),paste("ME", rownames(mod.roi), sep="")),
                     match(colnames(modulePvalue), paste("ME", colnames(mod.roi), sep=""))]
  
  cols = colnames(mod.roi) 
  
  mod.roi = cbind(rownames(mod.roi), mod.roi)
  colnames(mod.roi) = c("module", cols)  
  rownames(mod.roi) <- NULL
  
  
  return(mod.roi)
}



#' Compute correlation between patient rank in two tissues
#' @param subtypes is a vector of subtypes
#' @export 
patientRankCorrelation <- function(tissueA="blood", tissueB="biopsy"){
  MEs = computeEigengenes(c(tissueA,tissueB)) 
  
  moduleCor = cor(MEs[[tissueA]], MEs[[tissueB]], use = "p");
  modulePvalue = WGCNA::corPvalueStudent(moduleCor, ncol(dat[[tissueA]]$exprs));
  
  rank <- NULL
  
  ### relate patient ordering across tissues
  for (tissue in c(tissueA, tissueB)){
    rank[[tissue]] <- do.call(cbind, lapply(bresat[[tissue]], '[', "pat.order"))[,]
  }
  
  rank.cor.p<-plyr::laply(rank[[tissueA]], function (x){
    plyr::laply(rank[[tissueB]], function (y){
      WGCNA::corPvalueStudent(cor(x,y, use="p"), ncol(dat[[tissueA]]$exprs))
    })
  })
  
  rownames(rank.cor.p) <- names(rank[[tissueA]])
  colnames(rank.cor.p) <- names(rank[[tissueB]])
  
  rank.cor.p<- rank.cor.p[match(rownames(modulePvalue), paste("ME", rownames(rank.cor.p), sep="")),
                          match(colnames(modulePvalue), paste("ME", colnames(rank.cor.p), sep=""))] 
  cols = colnames(rank.cor.p) 
  rank.cor.p = cbind(rownames(rank.cor.p), rank.cor.p)
  colnames(rank.cor.p) = c("module", cols)  
  rownames(rank.cor.p) <- NULL
  
  return(rank.cor.p)
}

#' Get available cohorts
#' @export
getCohorts <- function(){
  tissue="blood"
  pat.dat <- NULL  
  pat.dat[[tissue]]$cohorts <- pat.cohorts(dat[[tissue]])
  cohorts<- names(pat.dat[[tissue]]$cohorts)
  return(cohorts)
}

#' Calculate patient rank sum scores for given tissues. 
#' @export 
patientRankSum <- function(tissueA="blood",tissueB="biopsy",cohort="all") { 
    tissues = c(tissueA,tissueB)
    
    pat.dat <- NULL  
    for(tissue in tissues){
      pat.dat[[tissue]]$cohorts <- pat.cohorts(dat[[tissue]])
    }
    
    cohorts<- names(pat.dat[[tissue]]$cohorts)
    ranksum<-NULL

    for(tissue in tissues){
      ranksum$all[[tissue]]<-do.call(cbind, lapply(bresat[[tissue]], '[', "ranksum"))[,]
      ranksum[[cohort]][[tissue]]<-lapply(ranksum$all[[tissue]], function(x) x[pat.dat[[tissue]]$cohorts[[cohort]]])
    }
    
    ranksum.cor.p<-NULL
    ranksum.cor.p[[cohort]]<-plyr::laply(ranksum[[cohort]][[tissueA]], function (x){
      plyr::laply(ranksum[[cohort]][[tissueB]], function (y){
        WGCNA::corPvalueStudent(cor(x,y, use="p"), length(pat.dat[[tissue]]$cohorts[[cohort]]))
      })
    })
    rownames(ranksum.cor.p[[cohort]]) <- names(ranksum[[cohort]][[tissueA]])
    colnames(ranksum.cor.p[[cohort]]) <- names(ranksum[[cohort]][[tissueB]])
    
    
    cols = colnames(ranksum.cor.p[[cohort]]) 
    ranksum.cor.p[[cohort]] = cbind(rownames(ranksum.cor.p[[cohort]]), ranksum.cor.p[[cohort]])
    colnames(ranksum.cor.p[[cohort]]) = c("module", cols)  
    rownames(ranksum.cor.p[[cohort]]) <- NULL
    
    #vals<-NULL
    #vals[[subtype]] <- -log10(ranksum.cor.p[[subtype]])
    #vals[[subtype]][vals[[subtype]] > 10] <- 10
    
    return(ranksum.cor.p[[cohort]])
    
}

#' Compute all 4 different analyses for modules from two tissues. 
#' @param tissueA is the first tissue
#' @param tissueB is the second tissue
#' @param moduleA is a module from the first tissue
#' @param moduleB is a module from the second tissue
#' @export 
comparisonAnalyses <- function(tissueA, tissueB, moduleA, moduleB){
  
  analyses = NULL
  eigen = eigengeneCorrelation(tissueA,tissueB) 
  rank = patientRankCorrelation(tissueA,tissueB)
  overlap = moduleHypergeometricTest(tissueA, tissueB) 
  roi = roiTest(tissueA, tissueB) 
  ranksum = patientRankSum(tissueA,tissueB,"all")
  
  analyses$ranksum = as.numeric(ranksum[ranksum[,1] == moduleA, colnames(ranksum) == moduleB])
  analyses$eigen =  as.numeric(eigen[eigen[,1] == moduleA , colnames(eigen) == moduleB])
  analyses$rank =  as.numeric(rank[rank[,1] == moduleA , colnames(rank) == moduleB])
  analyses$overlap =  as.numeric(overlap[overlap[,1] == moduleA , colnames(overlap) == moduleB])
  analyses$roi =  as.numeric(roi[roi[,1] == moduleA , colnames(roi) == moduleB])
  
  analyses$common =  intersect(bresat[[tissueA]][[moduleA]]$gene.order,bresat[[tissueB]][[moduleB]]$gene.order)
  return(analyses)
}

#' Computes correlation between eigengenes and quantitative variables
#' and anova between eigengenes and qualitative variables.  
#' @param tissue is the tissue from which to get the eigengenes 
#' @export 
eigengeneClinicalRelation <- function(tissue){
  
  bc.qual.var<- c("er","her2" ,"pam50.parker","hybrid.parker","cit", "claudin.low", 
                  "lum" ,"lumN","prolif" ,"basalL",
                  "lumC" ,"lumaNormL" ,"basLmApo","lumBlumC", 
                  "lymph" ,
                  "menopause","hrt","weight.70.plus","age.55.plus", "medication","hospital")
  bc.quant.var<-c("age","weight", "MKS", "ERS", "LUMS", "HER2S")
  n.qual.var<- c("menopause","hrt","weight.70.plus","age.55.plus", "medication","orig.dataset")
  bn.qual.var<- c("cancer", bc.qual.var, "orig.dataset")
  n.quant.var<-c("age","weight")
  
  
  MEs = computeEigengenes(tissue) 
  
  ### Qualitative variables 
  cl <- plyr::laply(dat[[tissue]]$clinical[,bc.qual.var], function(y) {
    plyr::laply(MEs[[tissue]], function (x){
      anova(lm(x~y))$`Pr(>F)`[1]
    })
  })
  rownames(cl) <- bc.qual.var
  colnames(cl) <- names(MEs[[tissue]])
  
  ### quantitative variables 
  moduleTraitCor<-NULL
  moduleTraitPvalue<-NULL
  moduleTraitCor[[tissue]]= cor(MEs[[tissue]], dat[[tissue]]$clinical[, bc.quant.var], use = "p");
  moduleTraitPvalue[[tissue]]= WGCNA::corPvalueStudent(moduleTraitCor[[tissue]], nrow(dat[[tissue]]$clinical));
  cl<-rbind(cl, t(moduleTraitPvalue[[tissue]]))
  
  
  if(tissue == "blood") { 
    orig.dataset<- plyr::laply(MEs[[tissue]], function (x){
      anova(lm(x~dat[[tissue]]$clinical[,3]))$`Pr(>F)`[1]
    })
    
    cl <- rbind(cl, orig.dataset)
  }
  
  if(tissue == "nblood"){
    cl <-plyr::laply(dat[[tissue]]$clinical[ , n.qual.var], function(y) {
      plyr::laply(MEs[[tissue]], function(x){
        anova(lm(x~y))$`Pr(>F)`[1]})
    })
    
    rownames(cl) < -n.qual.var
    colnames(cl) <- names(MEs[[tissue]])
    
  }
  cols = colnames(cl) 
  cl = cbind(rownames(cl), cl)
  colnames(cl) = c("Clinical", cols) 
  #res[[tissue]] <- -log10(cl)
  #res[[tissue]][res[[tissue]] > 10] <- 10
  return (cl) 
} 

#' Computes ROI for the given tissue
#' @export 
computeROICategories <- function(tissue) {
    roi.cat<-NULL
    module.names <- names(bresat[[tissue]])
    roi.cat[[tissue]]<- parallel::mclapply(module.names, function(module) {
      define.roi.regions(bresat[[tissue]][[module]], bresat[[tissue]][[module]]$roi)
    })
    names(roi.cat[[tissue]])<-module.names
    return(roi.cat) 
}

#' ROI and clinical categories relation 
#' @export 
roiClinicalRelation <- function(tissue){
  
  roi.cat <- computeROICategories(tissue) 
  MEs = computeEigengenes(tissue) 
  
  bc.qual.var<- c("er","her2" ,"pam50.parker","hybrid.parker","cit", "claudin.low", 
                  "lum" ,"lumN","prolif" ,"basalL",
                  "lumC" ,"lumaNormL" ,"basLmApo","lumBlumC", 
                  "lymph" ,
                  "menopause","hrt","weight.70.plus","age.55.plus", "medication","hospital")
  bc.quant.var<-c("age","weight", "MKS", "ERS", "LUMS", "HER2S")
  n.qual.var<- c("menopause","hrt","weight.70.plus","age.55.plus", "medication","orig.dataset")
  bn.qual.var<- c("cancer", bc.qual.var, "orig.dataset")
  n.quant.var<-c("age","weight")
  
  
  cl.roi <- NULL
  # fisher's exact betgween roi and qualitative variables 
  cl.roi<-plyr::laply(dat[[tissue]]$clinical[,bc.qual.var], function(y) {
    plyr::laply(data.frame(roi.cat[[tissue]]), function (x) {
        #fisher.test(y,x, workspace=2e+8,hybrid=TRUE)$p
        chisq.test(table(y,x))$p.value
        }, .parallel=FALSE)
    }, .parallel=FALSE)
    
  rownames(cl.roi) <- bc.qual.var
  colnames(cl.roi) <- names(roi.cat[[tissue]])
  
  roiTraitPvalue<-NULL
  
  roiTraitPvalue[[tissue]]<- plyr::laply(dat[[tissue]]$clinical[,bc.quant.var], function(y) {
    plyr::laply(MEs[[tissue]], function (x){
      anova(lm(x~y))$`Pr(>F)`[1]
    })
  })
  
  rownames(roiTraitPvalue[[tissue]]) <- bc.quant.var
  colnames(roiTraitPvalue[[tissue]]) <- names(MEs[[tissue]])
  
  cl.roi <- rbind(cl.roi, roiTraitPvalue[[tissue]])

  
  if(tissue == "blood") { 
    orig.dataset<-NULL
    orig.dataset<- plyr::laply(roi.cat[[tissue]], function (x){
      fisher.test(x,dat[[tissue]]$clinical[,3])$p
    })
    
    cl.roi<-rbind(cl.roi, orig.dataset)
  }
  
  if(tissue == "nblood"){ 
    cl.roi<-plyr::laply(dat[[tissue]]$clinical[,n.qual.var], function(y) {
      plyr::laply(data.frame(roi.cat[[tissue]]), function(x){
        #fisher.test(y,x, workspace=2e+07, hybrid=TRUE)$p
        chisq.test(table(y,x))$p.value
        })
    })
  
    rownames(cl.roi)<-names(dat[[tissue]]$clinical)[n.quant.var]
    colnames(cl.roi)<-names(roi.cat[[tissue]])
  }
  
  cl.roi <- cl.roi[,match(names(MEs[[tissue]]), paste("ME",colnames(cl.roi), sep=""))]
  
  cols = colnames(cl.roi) 
  cl.roi = cbind(rownames(cl.roi), cl.roi)
  colnames(cl.roi) = c("Clinical", cols) 
  
  
  
  return(cl.roi)
}

#' Get graph nodes for a TOM graph for a given tissue 
#' @export
getTOMGraphNodes <- function(tissue) {

  edges = getTOMGraphEdges(tissue)
  g = igraph::graph_from_data_frame(edges)
  
  layout = igraph::layout_with_fr(g, niter=15)
  nodes <- net[[tissue]]$nodeData
  nodes = nodes[nodes$nodeName %in% igraph::V(g)$name, ]
  nodes = cbind(nodes, layout)
  colnames(nodes) <- c("id", "altId", "color", "x", "y")
  return(nodes)
}

#' Get graph edges for a TOM graph for a given tissue 
#' @export
getTOMGraphEdges <- function(tissue){
  edges = net[[tissue]]$edgeData
  edges = dplyr::mutate(edges, id=paste0(fromNode,"-",toNode))
  nodes <- net[[tissue]]$nodeData
  #edges = mutate(edges, module=nodes[nodes$nodeName == fromNode,]$color)
  #edges = edges[edges$weight > 0.15, ]
  colnames(edges) <- c("source", "target", "weight", "direction", "sourceAltId", "targetAltId", "id")
  return(edges)
}


