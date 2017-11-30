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
#' @param cohort the patient cohort we select patients from, defaults to all patients
#' @import colorspace
#' @export
#' @examples
#' heatmap("blood", "green")
#' heatmap("biopsy", "blue")
#'@export
cohort_heatmap <- function(tissue, module, cohort.name="all", patient.ids=NULL,
													 gene.names=NULL,  orderByModule=NULL,
													 orderByTissue=NULL, cl.height=6, title=title) {
	return (mixtR::cohort_heatmap(mixt.dat=dat, mixt.ranksum = bresat,
															 tissue = tissue, module = module,
															 cohort.name=cohort.name, orderByModule=orderByModule,
															 orderByTissue = orderByTissue, cl.height=cl.height))
}

# Generate cohort scatterplot.
# Needs some refactoring re: variable names etc.
#' @export
cohort_scatterplot <-  function (x.tissue, x.module, y.tissue, y.module,
												cohort.name = "all") {
	return(mixtR::cohort_scatterplot(bresat, x.tissue, x.module,
																	 y.tissue, y.module, cohort.name))
}

#' Generate boxplot.
#' @export
cohort_boxplot<-function (tissue, module, orderByTissue, orderByModule,
													cohort.name="all", patient.ids=NULL){

	return(mixtR::cohort_boxplot(bresat, tissue, module, cohort.name,
															 orderByTissue, orderByModule))
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
	tissues = getAllTissues()
	genes = NULL
	for(tissue in tissues){
		genes = c(genes, names(moduleColors[[tissue]]))
	}
  return(unique(genes))
}

#' Get all modules a specific gene is found in.
#' @param gene the interesting gene
#' @return vector with module name from blood and module name from biopsy
#' @export
getAllModules <- function(gene) {
	tissues = getAllTissues()
  modules = NULL
  for(tissue in tissues) {
  	id = match(gene, names(moduleColors[[tissue]])) # same id blood or biopsy
  	module = as.character(moduleColors[[tissue]][id])
  	modules = c(modules, module)
  }
  return(modules)
}

#' Retrieves an overview of all genes and the modules they participate in.
#' @export
getAllGenesAndModules <- function(cohort="all") {
  res <- NULL
  tissues <- getAllTissues()
  for (tissue in tissues){
    for(module in names(bresat[[tissue]])) {
      if(module == "grey"){
        next
      }
      gs <- bresat[[tissue]][[module]][[cohort]]$gene.order
      for(gene in gs){
        gene = as.character(gene)
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
getGeneList <- function(tissue, module, cohort="all"){
  genes <- rownames(bresat[[tissue]][[module]][[cohort]]$dat)
  up.dn <- bresat[[tissue]][[module]][[cohort]]$up.dn

  res <- matrix(c(genes,up.dn), nrow=length(genes))
  colnames(res) <- c("Gene", "up.dn")
  res = data.frame(res)

  # get correlation and merge
  a = matrix(unlist(c(bresat[[tissue]][[module]][[cohort]]$up.cor, bresat[[tissue]][[module]][[cohort]]$dn.cor)))
  colnames(a) <- c("cor")
  df = data.frame(a)
  res = cbind(res, df)
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
  res$updn.pval = as.character(res$updn.pval)
  return(res)
}

#' Get gene set names available to the MIxT app
#' @export
getGeneSetNames <- function() {
  return (names(msigdb.enrichment[[1]][[1]]$updn.common))
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
userEnrichmentScores <- function(tissue, genelist, cohort="all") {
  modules = names(bresat[[tissue]])

  all_genes = names(moduleColors[[tissue]])
  genelist = genelist[genelist %in% all_genes]
  results <- lapply(modules, function(module){
    mod.genes <- names(moduleColors[[tissue]])[which(moduleColors[[tissue]]==module)]
    s <- length(intersect(mod.genes, all_genes))
    e <- length(intersect(genelist,all_genes))
    com <-length(intersect(mod.genes, genelist))
    intersections <- intersect(mod.genes, genelist)
    p_val<- sum(dhyper(com:e,s,length(all_genes)-s, e))
  return(list(p_values=p_val, common=intersections))
  })

  p_values <- unlist(lapply(results, function (x) x$p_values))
  common <- lapply(results, function (x) x$common)
  names(common) <- modules

  ret <- data.frame(p_values=p_values, module=modules, stringsAsFactors=F)
  ret$common <- common
  return(ret)
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
	tissue = names(moduleColors[1])
	module = getModules(tissue)[1]
  return(goterms[[tissue]][[module]]$GO.table$Term)
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



#' Compute gene overlap between genes from two tissues
#' @import WGCNA
#' @export
geneOverlapTest <- function(tissueA=NULL, tissueB=NULL){

		if(is.null(tissueA)){
			tissueA = names(moduleColors)[1]
		} else if(is.null(tissueB)){
			tissueB = names(moduleColors)[2]
		}

    tissue.overlap = NULL
    tissue.overlap = WGCNA::overlapTable(moduleColors[[tissueA]], moduleColors[[tissueB]])
    adjusted = NULL
    adjusted = matrix(p.adjust(tissue.overlap$pTable, method="BH"),
                         nrow = nrow(tissue.overlap$pTable),
                         ncol=ncol(tissue.overlap$pTable))
    adjusted = as.data.frame(adjusted)
    rownames(adjusted) = rownames(tissue.overlap$pTable)
    adjusted = cbind(rownames(adjusted), adjusted)
    colnames(adjusted) = c("module", colnames(tissue.overlap$pTable))
    return(adjusted)
}


#' Get available cohorts
#' @export
getCohorts <- function(tissue=NULL){
	if(is.null(tissue)){
		tissue = names(moduleColors)[1]
	}
  cohorts = names(dat[[tissue]]$cohorts)
  return(cohorts)
}

#' Calculate patient rank sum scores for given tissues.
#' @export
patientRankSum <- function(tissueA=NULL,tissueB=NULL,cohort="all") {

	if(is.null(tissueA)){
		tissueA = names(moduleColors)[1]
	} else if(is.null(tissueB)){
		tissueB = names(moduleColors)[2]
	}

  # The p-values have already been computed so we can just extract them
  # from the data frame.

  if(tissueA != tissueB & tissueA == names(moduleColors)[1]){
        object <-  paste(tissueA, tissueB, sep="_")
        correlation_p_value <- perm_cor_p[[object]][[cohort]]
    }
    if(tissueA != tissueB & tissueB == names(moduleColors)[1]){
        object = paste(tissueB, tissueA, sep="_")
        correlation_p_value <- t(perm_cor_p[[object]][[cohort]])
    }
    if(tissueA == tissueB){
	object = paste0(tissueA, "2")
	correlation_p_value = perm_cor_p[[object]][[cohort]]
    }

  correlation_p_value = as.data.frame(correlation_p_value)

  cols = colnames(correlation_p_value)
  correlation_p_value = cbind(rownames(correlation_p_value), correlation_p_value)
  colnames(correlation_p_value) = c("module", cols)
  rownames(correlation_p_value) = NULL
  return(correlation_p_value)
}

#' Compute all 4 different analyses for modules from two tissues.
#' @param tissueA is the first tissue
#' @param tissueB is the second tissue
#' @param moduleA is a module from the first tissue
#' @param moduleB is a module from the second tissue
#' @export
comparisonAnalyses <- function(tissueA, tissueB, moduleA, moduleB, cohort="all"){

  analyses = NULL
  overlap = geneOverlapTest(tissueA, tissueB)
  ranksum = patientRankSum(tissueA,tissueB, cohort)
  
  analyses$ranksum = as.numeric(ranksum[ranksum[,1] == moduleA, colnames(ranksum) == moduleB])
  analyses$overlap =  as.numeric(overlap[overlap[,1] == moduleA , colnames(overlap) == moduleB])
  analyses$common =  intersect(rownames(bresat[[tissA]][[modA]][[cohort]]$dat),rownames(bresat[[tissB]][[modB]][[cohort]]$dat))

  return(analyses)
}

#' Compute clinical ranksum for given tissue
#' @export
clinicalRanksum <- function(tissue, cohort="all") {
  ### Now results are pre-computed. Preseting now FDR stat adjusted for multiple testing
  clinicalVars = row.names(mod_clinical_fdr[[cohort]][[tissue]])
  cols = colnames(mod_clinical_fdr[[cohort]][[tissue]])
  mod_clinical_fdr[[cohort]][[tissue]] = cbind(clinicalVars, mod_clinical_fdr[[cohort]][[tissue]])
  colnames(mod_clinical_fdr[[cohort]][[tissue]]) = c("Clinical", cols)
  select.var<-c("lymph", "er", "MKS","pam50.parker", "hybrid", "cit",
                "lumC", "t.size", "claudin.low", "weight", "LUMS", "hrt",
                "her2", "HER2S", "age", "menopause", "medication")

  tmp = NULL
  tmp = mod_clinical_fdr[[cohort]][[tissue]]
  tmp = tmp[rownames(tmp) %in% select.var, ]

       return(tmp)
    }


#' Get graph nodes for a TOM graph for a given tissue
#' @export
getTOMGraphNodes <- function(tissue) {

  n<-network::network(net[[tissue]]$edgeData[,1:2], directed=F)

  #  n %e% "weight" <- net[[tissue]]$edgeData$weight
  network::set.edge.attribute(n, "weight", net[[tissue]]$edgeData$weight)

  x = network::network.vertex.names(n)
  #n %v% "module"<- moduleColors[[tissue]][match(x,names(moduleColors[[tissue]]))]
  network::set.vertex.attribute(n, "module", moduleColors[[tissue]][match(x,names(moduleColors[[tissue]]))])

  g = GGally::ggnet2(n, color="module", mode = "fruchtermanreingold", label=T, label.size=1,label.color="grey60", layout.par = list(cell.jitter = 0.2),edge.size="weight", alpha=0.75, size=3)

  nodes = NULL
  nodes = cbind(g$data$label, g$data$color, g$data$x, g$data$y)
  colnames(nodes) <- c("id", "color", "x", "y")

  nodes = as.data.frame(nodes)

  return(nodes)
}

#' Get graph edges for a TOM graph for a given tissue
#' @export
getTOMGraphEdges <- function(tissue){
  edges = net[[tissue]]$edgeData
  edges = dplyr::mutate(edges, id=paste0(fromNode,"-",toNode))
  nodes <- net[[tissue]]$nodeData
  #edges = mutate(edges, module=nodes[nodes$nodeName == fromNode,]$color)
  #edges = edges[edges$weight > 0.1, ]
  colnames(edges) <- c("source", "target", "weight", "direction", "sourceAltId", "targetAltId", "id")
  return(edges)
}

getEdgeColor <- function(edge,tissue){
  color = as.character(moduleColors[[tissue]][[edge[1]]])
  color2 = as.character(moduleColors[[tissue]][[edge[2]]])
  return (color)
}

getNodeColor <- function(node, tissue){
  color = as.character(moduleColors[[tissue]][[node[1]]])
}

edgeColors <- function(a){
  return(a)
}
