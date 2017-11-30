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
#' @import BCIutils
#' @import huc
#' @export
#' @examples
#' heatmap("blood", "green")
#' heatmap("biopsy", "blue") 
#'@export
cohort_heatmap <- function(tissue, module, cohort.name="all", patient.ids=NULL, gene.names=NULL,  orderByModule=NULL, orderByTissue=NULL, cl.height=6, title=title) 
  {
    plot.new()
    title = "" 
    if(is.null(orderByModule)){
      orderByModule = module
    }
    if(is.null(orderByTissue)){
      orderByTissue = tissue
    }

    if(orderByModule == module && orderByTissue == tissue) {
      title = paste(module," module from ",tissue, sep="") 
    } else {
      title = paste(cohort.name, module, tissue, "ordered by", orderByModule, orderByTissue)
    }
    
    
    #heatmap variables
    col.clust = FALSE
    layout.m = matrix(c("key","title","","","",
                        "","","","","",
                        "ranksum.text","ranksum.line","","","",
                        "row.labels.rjust","heatmap","","","" ,
                        "","","","","",
                        "ranks.text","ranks","","","",
                        "","","","","",
                        "clinical.labels.rjust","clinical","","","",
                        "","","","",""),nrow=9,ncol=5,byrow=TRUE)
    
    layout.m.sum = matrix(c("","","","",""  ,  "","","","",""  ,"","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "row.labels.rjust","heatmap","","",""  ,  "","","","",""  ,   "","","","",""  ,  "","","","",""),nrow=9,ncol=5,byrow=TRUE)
    layout.m.updn = matrix(c("","","","",""  ,  "","","","",""  ,"","","","",""  ,  "","","","heatmap",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""),nrow=9,ncol=5,byrow=TRUE)
    widths = c(2,5,0.25,0.25,0.25,0.25)
    heights = c(1,0.25,0.5,3,0.25,0.25,0.25,cl.height,0.25)
    
    scale = "none"
    min.val=-5
    max.val=5
    key.min=-5
    key.max=5
    
    ## define reordered variables
    bs.order.by <- bresat[[orderByTissue]][[orderByModule]][[cohort.name]]
    
    order.by<-bs.order.by$pat.order
    roi<-bs.order.by$roi
    roi.cat<-bs.order.by$roi.cat
    
    ## define patients to include
    if (is.null(patient.ids)) {
      patients <- pat.cohorts(dat$bnblood)[[cohort.name]]
    } else {
      patients <- patient.ids
    }

    # define clinical and exprs ------------------    
    ## define reordered clinical data    
    mclinical = NULL
    cl = NULL
    cl<-dat[[tissue]]$clinical[rownames(dat[[tissue]]$clinical) %in% patients,]
    mclinical = cl[order.by,]
    
    ## define reordered expression and select genes
    bs <- bresat[[tissue]][[module]][[cohort.name]]
    data = bs$dat[, match(rownames(mclinical), colnames(bs$dat))]
    if (!is.null(gene.names)){
      data<-data[rownames(data) %in% gene.names,]
    }
        
    bc.var<- c("er","her2" ,"pam50.parker", "hybrid", "cit", "intclust",
               "claudin.low", 
               "lymph" , "t.size",
               "menopause","hrt","medication",
               "age","weight", "MKS", "LUMS", "HER2S")
    if(tissue=="bnblood"){
      bc.var = c("cancer",bc.var)
    }
    
    ## plot heatmap
    ddrs = heatmap.simple(data,
                          clinical = huc.color.clinical(mclinical)[,bc.var], 
                          layout.mat = layout.m, widths = widths, heights = heights, col.clust = FALSE, 
                          row.clust = FALSE, title= title,
                          row.labels=rownames(data),
                          col.labels=rep("", ncol(data)))
    
    ## plot updn top left heatmap
    up.dn = as.vector(array(1,dim=c(1,length(bs$gene.order))))
    names(up.dn) = rownames(bs$dat)
    if(length(bs$dn)>0){
      up.dn[names(up.dn) %in% rownames(dat[[tissue]]$exprs)[bs$dn]] = -1}
    to.plot = (as.matrix(up.dn,ncol=1))
    color.scheme = heatmap.color.scheme(low.breaks=c(-1.5,0),high.breaks=c(0,1.5))
    heatmap.simple(to.plot, scale=scale, layout.mat = layout.m.updn, widths = widths, heights = heights, col.clust = FALSE, row.clust = FALSE, color.scheme = color.scheme)
    
    ## plot ranks for top left heatmap
    the.layout <- grid.layout(nrow(layout.m), ncol(layout.m), widths=widths, heights=heights)
    mid.vp <- viewport(layout=the.layout, name="heatmap.mid.vp")
    pushViewport(mid.vp)
    elem = 'ranks'
    idx <- which(layout.m == elem, arr.ind=TRUE)
    pushViewport(viewport(name=elem,
                          layout.pos.row=unique(idx[,1]),
                          layout.pos.col=unique(idx[,2])))
    
    rank.colors<-rev(diverge_hcl(n=ncol(bs$dat)))
    names(rank.colors)<-colnames(bs$dat)
    rank.colors<-rank.colors[match(rownames(mclinical), colnames(bs$dat))]
    ranksum = t(as.matrix(rank.colors))[,names(rank.colors) %in% patients,drop=FALSE]
    heatmap.clinical(ranksum)
    upViewport()
    elem = 'ranks.text'
    idx <- which(layout.m == elem, arr.ind=TRUE)
    pushViewport(viewport(name=elem,
                          layout.pos.row=unique(idx[,1]),
                          layout.pos.col=unique(idx[,2])))
    heatmap.labels('ranks', type="row.labels", just="right")
    upViewport()
    
    ## ranksum bottom left
    the.layout <- grid.layout(nrow(layout.m), ncol(layout.m), widths=widths, heights=heights)
    top.vp <- viewport(layout=the.layout, name="heatmap.top.vp")
    pushViewport(top.vp)
    ranksum.plot = bs$ranksum[order.by][colnames(data) %in% patients]
    elem = 'ranksum.line'
    idx <- which(layout.m == elem, arr.ind=TRUE)
    pushViewport(viewport(name=elem,
                          layout.pos.row=unique(idx[,1]),
                          layout.pos.col=unique(idx[,2]),
                          xscale=c(0.5, length(ranksum.plot) + 0.5),
                          yscale=range(ranksum.plot)))
    
    grid.rect(gp=gpar(lwd=0.1))
    grid.polyline(rep(c(0, 1), 4), rep(c(0.2, 0.4, 0.6, 0.8), each=2), id.lengths=rep(2, 4), gp=gpar(lwd=0.1, col="grey70"))
    xrange <- range(1:length(ranksum.plot))
    n <- length(ranksum.plot)
    grid.segments(unit(1:length(ranksum.plot),"native"), rep(0,n), unit(1:length(ranksum.plot),"native"),unit(ranksum.plot, "native"))
    upViewport()
    
    elem = 'ranksum.text'
    idx <- which(layout.m == elem, arr.ind=TRUE)
    pushViewport(viewport(name=elem,
                          layout.pos.row=unique(idx[,1]),
                          layout.pos.col=unique(idx[,2])))
    
    heatmap.labels('ranksum', type="row.labels", just="right")
    upViewport()
    
    ## plot roi lines top left heatmap
    first.ind = length(which(roi.cat[rownames(cl) %in% patients]==3))
    last.ind = first.ind + length(which(roi.cat[rownames(cl) %in% patients]==2))
    
    res.random.dist.begin = first.ind
    res.random.dist.end = last.ind
    
    elem = 'heatmap'
    idx <- which(layout.m == elem, arr.ind=TRUE)
    pushViewport(viewport(name=elem,
                          layout.pos.row=unique(idx[,1]),
                          layout.pos.col=unique(idx[,2])))
    
    par(new=TRUE, fig=gridFIG(), mar=c(0,0,0,0))
    grid.lines(x = unit(c(first.ind/length(which(rownames(cl) %in% patients)),first.ind/length(which(rownames(cl) %in% patients))), "npc"),gp=gpar(col='yellow',lty = 1, lwd = 2))
    grid.lines(x = unit(c(last.ind/length(which(rownames(cl) %in% patients)),last.ind/length(which(rownames(cl) %in% patients))), "npc"),gp=gpar(col='yellow',lty = 1, lwd = 2))
    upViewport()
}

# Generate cohort scatterplot. 
# Needs some refactoring re: variable names etc. 
#' @export
cohort_scatterplot <-
  function (x.tissue,
            x.module,
            y.tissue,
            y.module,
            cohort.name = "all",
            patient.ids = NULL) {
    ## define patients to include
    if (is.null(patient.ids)) {
      patients <- pat.cohorts(dat[[x.tissue]])[[cohort.name]]
    } else {
      patients <- patient.ids
    }
    # we need to swap around these to fix accessing the per.cor.p object
    if (x.tissue != "blood") {
      tmp = NULL
      xmodule = NULL
      ymodule = NULL
      tmp = x.module
      xmodule = y.module
      ymodule = tmp
      xtissue = y.tissue
      ytissue = x.tissue
    } else {
      xmodule = x.module
      ymodule = y.module
      xtissue = x.tissue
      ytissue = y.tissue
    }
    
    if(x.tissue==y.tissue){
      comp=paste0(xtissue, 2)
    } else {
      comp=paste0(xtissue,"_", ytissue)
    }
    
    
    ### data for scatterplot
    plot.data <-
      data.frame(
        x.ranksum = bresat[[x.tissue]][[x.module]][[cohort.name]]$ranksum,
        y.ranksum = bresat[[y.tissue]][[y.module]][[cohort.name]]$ranksum,
        cohort = rep(cohort.name, length(patients)),
        subtype = huc::huc.color.clinical(dat[[x.tissue]]$clinical)$hybrid[rownames(dat[[x.tissue]]$clinical) %in% patients]
      )
    
    sub.col <-
      c(
        normal = "white",
        all = "grey",
        erp = "green",
        ern = "firebrick2",
        her2p = "hotpink2",
        her2n = "#21B6A8",
        erp.her2p = "orange",
        ern.her2p = "hotpink2",
        erp.her2n = "blue",
        ern.her2n = "firebrick2",
        luma = "blue4",
        erp.luma = "blue4",
        lumb = "deepskyblue",
        erp.lumb = "deepskyblue",
        normL = "green4",
        erp.normL = "green4",
        basalL = "firebrick2",
        her2E = "hotpink2",
        erp.her2E = "orange",
        erp.basalL = "#7fffd4",
        cit.luma = "blue4",
        cit.lumb = "deepskyblue",
        cit.normL = "green4",
        cit.mApo = "hotpink2",
        cit.lumc = "#7fffd4",
        cit.basalL = "firebrick2",
        intclust1="deepskyblue",
        intclust3="blue4", 
        intclust7="blue4", 
        intclust8="blue4",
        intclust9="deepskyblue", 
        intclust4="green4",
        intclust4.ern="green4",
        intclust4.erp="green4", 
        intclust5="hotpink2",
        intclust10="firebrick2"
      )
    plot.data$sub.col <- sub.col[as.character(plot.data$cohort)]
    plot.data$sub.col <-
      ifelse(plot.data$sub.col == "grey",
             as.character(plot.data$subtype),
             plot.data$sub.col)
    
    p1 <-
      ggplot2::ggplot(plot.data, ggplot2::aes(x = x.ranksum, y = y.ranksum)) +
      ggplot2::geom_smooth(
        method = "lm",
        colour = "white",
        alpha = 0.2,
        size = 0.4
      ) +
      ggplot2::geom_point(ggplot2::aes(colour = sub.col), size = 2) +
      ggplot2::scale_colour_manual(values = levels(factor(plot.data$sub.col))) +
      ggplot2::labs(
        y = paste(y.module, y.tissue, "module ranksum"),
        x = paste(x.module, x.tissue, "module ranksum"),
        title = paste(
          cohort.name,
          " patients",
          " (cor=",
          as.character(signif(
            cor.test(plot.data$x.ranksum, plot.data$y.ranksum)$estimate,
            digits = 1
          )),
          ", p=",
          signif(perm.cor.p[[comp]][[cohort.name]][xmodule, ymodule], digits = 3),
          ")",
          sep = ""
        )
      ) +
      ggplot2::theme(
        legend.position = "none",
        panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        axis.line.x   = ggplot2::element_line(colour = "grey60"),
        axis.line.y   = ggplot2::element_line(colour = "grey60"),
        axis.title = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(
          hjust = 0,
          vjust = 1,
          size = 12
        )
      )
    
    print(p1)
  }

#' Generate boxplot. 
#' @export 
cohort_boxplot<-function (blood.module, orderByTissue, orderByModule, cohort.name="all", patient.ids=NULL){
  if(is.null(orderByModule)){
    orderByModule = module
  }
  if(is.null(orderByTissue)){
    orderByTissue = tissue
  }
  roi.cat<-bresat[[orderByTissue]][[orderByModule]][[cohort.name]]$roi.cat
  
  ## define patients to include
  if (is.null(patient.ids))
  {
    patients<-pat.cohorts(dat$bnblood)[[cohort.name]]
  }else
  {
    patients<-patient.ids
  }
  
  bs <- bresat$bnblood[[blood.module]][[cohort.name]]
  bnclinical = dat$bnblood$clinical[rownames(dat$bnblood$clinical) %in% patients, ]
  
 plot.data<-data.frame(bnbl.ranksum=bs$ranksum[c(which(bnclinical$cancer==TRUE), which(bnclinical$cancer==FALSE))],
                                                   cohort=c(rep(cohort.name, length(which(bnclinical$cancer==TRUE))), rep("normal", length(which(bnclinical$cancer==FALSE)))),
                                                   roi.cat=c(roi.cat,  
                                                             rep(NA, length(which(dat$bnblood$clinical$cancer==FALSE)))))

  plot.data$cancer<-1
  plot.data$cancer<-ifelse( plot.data$roi.cat==1 & !is.na(plot.data$roi.cat), 4, as.character(plot.data$cancer))
  plot.data$cancer<-ifelse(plot.data$roi.cat==2 & !is.na(plot.data$roi.cat), 2, as.character(plot.data$cancer))
  plot.data$cancer<-ifelse(plot.data$roi.cat==3 & !is.na(plot.data$roi.cat), 1, as.character(plot.data$cancer))
  plot.data$cancer<-as.numeric(plot.data$cancer)

  plot.data$tumor.cat<-"control"
  plot.data$tumor.cat<-ifelse(plot.data$roi.cat==1 & !is.na(plot.data$roi.cat), "low", as.character(plot.data$tumor.cat))
  plot.data$tumor.cat<-ifelse(plot.data$roi.cat==2 & !is.na(plot.data$roi.cat), "mid", as.character(plot.data$tumor.cat))
  plot.data$tumor.cat<-ifelse(plot.data$roi.cat==3 & !is.na(plot.data$roi.cat), "high", as.character(plot.data$tumor.cat))
  plot.data$tumor.cat<-factor(plot.data$tumor.cat, levels=c("low", "mid", "high", "control"))
  plot.data$tumor.cat.ordered<-factor(plot.data$tumor.cat, levels=c("high", "mid", "low", "control"), ordered=T)

  sub.col <- c(normal="white", all="grey", erp="green", ern="firebrick2", her2p="hotpink2", her2n="#21B6A8",
               erp.her2p="orange", ern.her2p="hotpink2", erp.her2n="blue", ern.her2n="firebrick2",
               luma="blue4", erp.luma="blue4", lumb="deepskyblue", erp.lumb="deepskyblue", normL="green4", erp.normL="green4", basalL="firebrick2", her2E="hotpink2", erp.her2E="orange", erp.basalL="#7fffd4", 
               cit.luma="blue4", cit.lumb="deepskyblue", cit.normL="green4", cit.mApo="hotpink2", cit.lumc="#7fffd4",  cit.basalL="firebrick2",
               intclust1="deepskyblue", intclust3="blue4", intclust7="blue4", intclust8="blue4", intclust9="deepskyblue", 
               intclust4="green4", intclust4.ern="green4", intclust4.erp="green4",intclust5="hotpink2", intclust10="firebrick2")
  plot.data$sub.col<-sub.col[as.character(plot.data$cohort)]
  
  p<-ggplot2::ggplot(data = plot.data, ggplot2::aes(x=tumor.cat.ordered, y=bnbl.ranksum)) + 
    ggplot2::geom_boxplot(ggplot2::aes(fill=sub.col, alpha=1/cancer))+
    ggplot2::geom_boxplot(colour="black", outlier.shape ="+", outlier.size = 1, fill=NA)+
    ggplot2::geom_point(shape="+")+
    ggplot2::scale_fill_manual(values=levels(factor(plot.data$sub.col)))+
    ggplot2::labs(x=paste(orderByModule, orderByTissue, "ROI module category"),
         y=paste(blood.module, "blood module ranksum"), 
         title=paste(cohort.name, " patients and controls\n(aov p=", as.character(signif(anova(lm(plot.data$bnbl.ranksum~plot.data$tumor.cat))$`Pr(>F)`[1], digits=1)), ")", sep=""))+
    ggplot2::theme(legend.position="none",
          panel.background = ggplot2::element_rect(fill = "transparent",colour = NA),
          axis.line.x = ggplot2::element_line(colour="grey60"),
          axis.line.y = ggplot2::element_line(colour="grey60"),
          axis.title.x=ggplot2::element_text(hjust=0.2),
          axis.title=ggplot2::element_text(size=10),
          plot.title = ggplot2::element_text(hjust=0,vjust=1, size=12)
    )
  print(p)
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
geneOverlapTest <- function(tissueA="blood", tissueB="biopsy"){
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
getCohorts <- function(tissue="blood"){
  tissue="blood"
  pat.dat <- NULL  
  pat.dat[[tissue]]$cohorts <- pat.cohorts(dat[[tissue]])
  cohorts<- names(pat.dat[[tissue]]$cohorts)
  return(cohorts)
}

#' Calculate patient rank sum scores for given tissues. 
#' @export 
patientRankSum <- function(tissueA="blood",tissueB="biopsy",cohort="all") { 
  
  # The p-values have already been computed so we can just extract them 
  # from the data frame.
  
  if(tissueA == "biopsy" && tissueB=="blood") {
    correlation_p_value = t(perm.cor.p$blood_biopsy[[cohort]])
  }
  if(tissueA=="blood" && tissueB == "biopsy") {
    correlation_p_value = perm.cor.p$blood_biopsy[[cohort]]
  }
  if(tissueA=="biopsy" && tissueB == "biopsy") {
    correlation_p_value = perm.cor.p$biopsy2[[cohort]]
  }
  if(tissueA== "blood" && tissueB=="blood") {
    correlation_p_value = perm.cor.p$blood2[[cohort]]
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
  ranksum = patientRankSum(tissueA,tissueB,"all")
  
  analyses$ranksum = as.numeric(ranksum[ranksum[,1] == moduleA, colnames(ranksum) == moduleB])
  analyses$overlap =  as.numeric(overlap[overlap[,1] == moduleA , colnames(overlap) == moduleB])
  analyses$common =  intersect(rownames(bresat[[tissueA]][[moduleA]][[cohort]]$dat),rownames(bresat[[tissueB]][[moduleB]][[cohort]]$dat))

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
