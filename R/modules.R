remove.outliers<- function (exprs){
  altered = NULL
  res = (sapply(1:nrow(exprs), function(i){
    ratio = 4.652788
    gene.exprs = exprs[i,]
    mad = mad(gene.exprs)
    diff = abs(median(gene.exprs) - gene.exprs)
    indices = which(diff>(ratio*abs(mad)))
    out = diff[indices]
    NAs = matrix(FALSE,nrow=1,ncol=length(gene.exprs))
    colnames(NAs) = colnames(exprs)
    if (length(out) > 0) {
      NAs[indices] = TRUE
    }
    return(list(exprs=gene.exprs,NAs=NAs,gene=rownames(exprs)[i]))
  }))
  altered = do.call(rbind,res['NAs',])
  rownames(altered) = do.call(rbind,res['gene',])
  NAs = sort(rowSums(altered>0),decreasing=TRUE)
  to.change = exprs[rownames(altered),colnames(altered)]
  to.change[altered] = NA
  exprs[rownames(altered),colnames(altered)] = to.change
}

match.datasets <- function(ds1, ds2, include.normals=TRUE, remove.outliers=TRUE) {
  if(is.null(ds1$exprs) || is.null(ds1$exprs))
    stop("Both datasets should have expression matrix stored as dataset$exprs in them")
  if(is.null(ds1$probe.info$gene.name) || is.null(ds2$probe.info$gene.name))
    stop("Both datasets should have gene names stored as dataset$probe.info$gene.name in them")
  
  
  if(remove.outliers) {
    ds1.exprs <-remove.outliers(exprs=ds1$exprs)
    ds2.exprs <-remove.outliers(exprs=ds2$exprs)
  }
  
  ### renoming objects
  ds1.exprs = ds1$exprs
  ds1.clinical = ds1$clinical
  #ds1.clinical.heatmap = ds1$clinical.heatmap
  
  ds2.exprs = ds2$exprs
  ds2.clinical = ds2$clinical
  #ds2.clinical.heatmap = ds2$clinical.heatmap
  
  if(remove.outliers) {
    ds1.exprs <-remove.outliers(exprs=ds1.exprs)
    ds2.exprs <-remove.outliers(exprs=ds2.exprs)
  }
  
  ds1.genes = ds1$probe.info$gene.name
  ds2.genes = ds2$probe.info$gene.name
  
  ## create matched.ds[[ds$name]]$matched.tumor including only MATCHED BC CASES AND GENES from both datasets
  ds1.samples = colnames(ds1.exprs)[!ds1.clinical$Normal][ds1.clinical$orig.id[!ds1.clinical$Normal] %in% ds2.clinical$orig.id[!ds2.clinical$Normal]]
  ds2.samples = colnames(ds2.exprs)[!ds2.clinical$Normal][ds2.clinical$orig.id[!ds2.clinical$Normal] %in% ds1.clinical$orig.id[!ds1.clinical$Normal]]
  
  ds1.matched.genes = ds1.genes %in% ds2.genes
  ds2.matched.genes = ds2.genes %in% ds1.genes
  
  matched.ds = NULL
  matched.ds[[ds1$name]]$matched.tumor = NULL
  matched.ds[[ds2$name]]$matched.tumor = NULL
  
  matched.ds[[ds1$name]]$matched.tumor$exprs = ds1.exprs[ds1.matched.genes,sort(ds1.samples)]
  matched.ds[[ds1$name]]$matched.tumor$genes = ds1$genes[ds1.matched.genes,]
  matched.ds[[ds1$name]]$matched.tumor$probe.info = ds1$probe.info[ds1.matched.genes,]
  matched.ds[[ds1$name]]$matched.tumor$clinical = ds1.clinical[sort(ds1.samples),,drop=F]
  #matched.ds[[ds1$name]]$matched.tumor$clinical.heatmap = ds1.clinical.heatmap[sort(ds1.samples),,drop=F]
  matched.ds[[ds1$name]]$matched.tumor$name = ds1$name
  matched.ds[[ds2$name]]$matched.tumor$exprs = ds2.exprs[ds2.matched.genes,sort(ds2.samples)]
  matched.ds[[ds2$name]]$matched.tumor$genes = ds2$genes[ds2.matched.genes,]
  matched.ds[[ds2$name]]$matched.tumor$probe.info = ds2$probe.info[ds2.matched.genes,]
  matched.ds[[ds2$name]]$matched.tumor$clinical = ds2.clinical[sort(ds2.samples),,drop=F]
  #matched.ds[[ds2$name]]$matched.tumor$clinical.heatmap = ds2.clinical.heatmap[sort(ds2.samples),,drop=F]
  matched.ds[[ds2$name]]$matched.tumor$name = ds2$name
  
  ## create matched.ds[[ds$name]]$normal including only controls and MATCHED GENES from both datasets
  if(include.normals){ 
    ds1.samples = colnames(ds1.exprs)[ds1.clinical$Normal]
    ds2.samples = colnames(ds2.exprs)[ds2.clinical$Normal]
    
    matched.ds[[ds1$name]]$normal = NULL
    matched.ds[[ds2$name]]$normal = NULL
    matched.ds[[ds1$name]]$normal$exprs = ds1.exprs[ds1.matched.genes,sort(ds1.samples)]
    matched.ds[[ds1$name]]$normal$genes = ds1$genes[ds1.matched.genes,]
    matched.ds[[ds1$name]]$matched.tumor$probe.info = ds1$probe.info[ds1.matched.genes,]
    matched.ds[[ds1$name]]$normal$clinical = ds1.clinical[sort(ds1.samples),,drop=F]
    #matched.ds[[ds1$name]]$normal$clinical.heatmap = ds1.clinical.heatmap[sort(ds1.samples),,drop=F]
    matched.ds[[ds1$name]]$normal$name = ds1$name
    matched.ds[[ds2$name]]$normal$exprs = ds2.exprs[ds2.matched.genes,sort(ds2.samples)]
    matched.ds[[ds2$name]]$normal$genes = ds2$genes[ds2.matched.genes,]
    matched.ds[[ds2$name]]$matched.tumor$probe.info = ds2$probe.info[ds2.matched.genes,]
    matched.ds[[ds2$name]]$normal$clinical = ds2.clinical[sort(ds2.samples),,drop=F]
    #matched.ds[[ds2$name]]$normal$clinical.heatmap = ds2.clinical.heatmap[sort(ds2.samples),,drop=F]
    matched.ds[[ds2$name]]$normal$name = ds2$name
  }
  return(matched.ds)
}  

################################
### Function to plot scale free topology criterion and mean connectivify given power (default=6 for unsigned network)
### choose soft-thesholding powers that (a) give approximate scale-free topology
### in each data set, and (b) give roughly comparable mean or median
### connectivities across the data sets. 

soft.threshold.exprs.data <- function (exprs, powers = c(c(1:10), seq(from = 12, to=20, by=2))) {
  sft = pickSoftThreshold(t(exprs), powerVector = powers, verbose = 1, networkType="unsigned")
  
  dev.new()
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  return(sft)
}

find.modules <- function (exprs, power = 6, deepSplit = 2, 
                          minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, pick.soft.threshold=FALSE){
  ### to interactively pick power for soft-threshold
  if (pick.soft.threshold) {
    soft.threshold.exprs.data(exprs)
    cat ("Please enter a valid integer as the soft threshold power to find modules (default=6) :")
    line <- readline()
    while(is.na(as.integer(line))) {
      cat ("Please enter a valid integer as the soft threshold power to find modules (default=6) :")
      line <- readline()
    }
    power = as.integer(line)
  }
  
  ### identify modules using the blockwise function
  gene.modules = NULL
  gene.modules = blockwiseModules(t(exprs), power = power, deepSplit = deepSplit, minModuleSize = minModuleSize, 
                                  reassignThreshold = reassignThreshold, mergeCutHeight = mergeCutHeight, 
                                  numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 3)
  
  moduleLabels = gene.modules$colors
  moduleColors = labels2colors(gene.modules$colors)
  MEs = gene.modules$MEs;
  geneTree = gene.modules$dendrograms;
  gene.modules = list(MEs=MEs, moduleLabels=moduleLabels, moduleColors=moduleColors, geneTree=geneTree, modules.GeneList=NULL, blockGenes=gene.modules$blockGenes)
  
  modules.gene.list = list()
  for (item in (unique(gene.modules$moduleColors))) {
    modules.gene.list[[item]] = rownames(exprs)[which(gene.modules$moduleColors == item)]
  }
  
  gene.modules$modules.GeneList = modules.gene.list 
  
  return(gene.modules)
}

load.modules <- function(dat, mod)
{
  mymodules$biopsy$exprs <- dat$biopsy$matched.tumor$exprs
  mymodules$blood$exprs <- dat$blood$matched.tumor$exprs
  mymodules$biopsy$clinical <- dat$biopsy$matched.tumor$clinical
  mymodules$blood$clinical <- dat$blood$matched.tumor$clinical
  mymodules$biopsy$heatmap.clinical <- huc.color.clinical(dat$biopsy$matched.tumor$clinical)
  mymodules$blood$heatmap.clinical <- huc.color.clinical(dat$blood$matched.tumor$clinical)
  mymodules$biopsy$modules <- mod$biopsy$modules.GeneList
  mymodules$blood$modules <- mod$blood$modules.GeneList
  ## make sure the grey module always comes first
  biopsy.grey.idx <- which(names(mymodules$biopsy$modules) == "grey")
  blood.grey.idx <- which(names(mymodules$blood$modules) == "grey")
  mymodules$biopsy$modules <- c(mymodules$biopsy$modules[biopsy.grey.idx], mymodules$biopsy$modules[-biopsy.grey.idx])
  mymodules$blood$modules <- c(mymodules$blood$modules[blood.grey.idx], mymodules$blood$modules[-blood.grey.idx])
  
  ### add normals
  mymodules$nblood<-mymodules$blood
  mymodules$nblood$exprs<-cbind(dat$blood$matched.tumor$exprs,dat$blood$normal$exprs)
  mymodules$nblood$clinical<-rbind(dat$blood$matched.tumor$clinical, dat$blood$normal$clinical)
  mymodules$nblood$clinical$Normal<-c(rep(FALSE, nrow(dat$blood$matched.tumor$clinical)), 
                                      rep(TRUE,nrow(dat$blood$normal$clinical)))
  mymodules$nblood$heatmap.clinical<-huc.color.clinical(mymodules$nblood$clinical)
  
  return(mymodules)
}

gene.overlap.test <- function(modules)
{
  all.genes <- intersect(unlist(modules$biopsy$modules[-1]), unlist(modules$blood$modules[-1]))
  
  pvals <- sapply(modules$biopsy$modules[-1], function(biopsy.mod) {
    sapply(modules$blood$modules[-1], function(blood.mod) {
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


create.modules.heatmap <- function(bs, exprs, clinical, re.order=TRUE, order.by, title=title) 
{
  stopifnot(ncol(bs$dat) == nrow(clinical))  
  
  #heatmap variables
  col.clust = FALSE
  layout.m = matrix(c("key","title","","",""  ,  "","","","",""  ,"ranksum.text","ranksum.line","","",""  ,  "row.labels.rjust","heatmap","","",""  ,  "","col.labels.rjust","","",""  ,  "","","","",""  ,  "","","","",""  ,  "clinical.labels.rjust","clinical","","",""  ,  "","","","",""),nrow=9,ncol=5,byrow=TRUE)
  layout.m.sum = matrix(c("","","","",""  ,  "","","","",""  ,"","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "row.labels.rjust","heatmap","","",""  ,  "","","","",""  ,   "","","","",""  ,  "","","","",""),nrow=9,ncol=5,byrow=TRUE)
  #layout.m.dist.cdf = matrix(c("","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","col.labels.ljust","","",""  ,  "row.labels.rjust","heatmap","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""),nrow=11,ncol=5,byrow=TRUE)
  layout.m.updn = matrix(c("","","","",""  ,  "","","","",""  ,"","","","",""  ,  "","","","heatmap",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""),nrow=9,ncol=5,byrow=TRUE)
  widths = c(2,5,0.25,0.25,0.25)
  heights = c(1,0.25,0.5,3,0.25,0.25,0.25,8,0.25)
  
  scale = "none"
  min.val=-5
  max.val=5
  key.min=-5
  key.max=5
  
  if (re.order == FALSE){
    order.by<-bs$pat.order
  }

  data = exprs[bs$gene.order, order.by, drop=FALSE]
  mclinical = clinical[order.by,]
  
  #plot.new()
  ddrs = heatmap.simple(data,layout.mat = layout.m, widths = widths, heights = heights, col.clust = FALSE, 
                        clinical = mclinical, row.clust = FALSE, title=title,row.labels=rownames(data), 
                        col.labels=rep("", ncol(data)))
  
  ranks = bs$pat.order
  names(ranks) = 1:length(bs$pat.order)
  ranks = as.integer(names(sort(ranks)))
  ranksum = t(as.matrix(ranks))[,order.by,drop=FALSE]
  rownames(ranksum) = "ranks"
  heatmap.simple(-ranksum, layout.mat = layout.m.sum, widths = widths, heights = heights, col.clust = FALSE, row.clust = FALSE)
  
  up.dn = as.vector(array(1,dim=c(1,length(bs$gene.order))))
  names(up.dn) = unique(c(bs$up,bs$dn))
  up.dn[names(up.dn) %in% bs$dn] = -1
  to.plot = (as.matrix(up.dn,ncol=1)[bs$gene.order,,drop=FALSE])
  color.scheme = heatmap.color.scheme(low.breaks=c(-1.5,0),high.breaks=c(0,1.5))
  
  heatmap.simple(to.plot, scale=scale, layout.mat = layout.m.updn, widths = widths, heights = heights, col.clust = FALSE, row.clust = FALSE, color.scheme = color.scheme)
  
  the.layout <- grid.layout(nrow(layout.m), ncol(layout.m), widths=widths, heights=heights)
  top.vp <- viewport(layout=the.layout, name="heatmap.top.vp")
  pushViewport(top.vp)
  elem = 'ranksum.line'
  idx <- which(layout.m == elem, arr.ind=TRUE)
  pushViewport(viewport(name=elem,
                        layout.pos.row=unique(idx[,1]),
                        layout.pos.col=unique(idx[,2])))
  ranksum.plot = bs$ranksum[order.by]
  par(new=TRUE, fig=gridFIG(), mar=c(0,0,0,0))
  plot(1:length(ranksum.plot), ranksum.plot, ann=FALSE, xaxs='i', yaxt='n', xaxt='n',bty='n',type='l')
  upViewport()
  elem = 'ranksum.text'
  idx <- which(layout.m == elem, arr.ind=TRUE)
  pushViewport(viewport(name=elem,
                        layout.pos.row=unique(idx[,1]),
                        layout.pos.col=unique(idx[,2])))
  
  heatmap.labels('ranksum', type="row.labels", just="right")
  upViewport()
  
  first.ind = length(which(bs$roi.cat==3))+1
  last.ind = first.ind + length(which(bs$roi.cat==2))
  
  res.random.dist.begin = first.ind
  res.random.dist.end = last.ind
  
  elem = 'heatmap'
  idx <- which(layout.m == elem, arr.ind=TRUE)
  pushViewport(viewport(name=elem,
                        layout.pos.row=unique(idx[,1]),
                        layout.pos.col=unique(idx[,2])))
  
  par(new=TRUE, fig=gridFIG(), mar=c(0,0,0,0))
  grid.lines(x = unit(c(first.ind/length(bs$roi),first.ind/length(bs$roi)), "npc"),gp=gpar(col='yellow',lty = 1, lwd = 2))
  grid.lines(x = unit(c(last.ind/length(bs$roi),last.ind/length(bs$roi)), "npc"),gp=gpar(col='yellow',lty = 1, lwd = 2))
  upViewport()
}


mod.p.heatmap<- function (mat.p, title, blood.biopsy=TRUE){
  col.scheme <- heatmap.color.scheme(low.breaks=c(), high.breaks=c(0, -log10(c(0.05, 0.01, 0.001, 0.0001))))
  mylayout<-matrix(c("key", "title", "","",
                     "", "", "", "",
                     "row.ddr.left", "heatmap", "row.labels.ljust", "rowSideLabels",
                     "", "col.labels.rjust", "", "",
                     "","colSideLabels", "", ""), ncol=4, byrow=TRUE)
  mywidths = c(1,5,1,1)
  myheights <- c(0.8,0.2,5,1)
  plot.new()
  vals <- -log10(mat.p)
  vals[vals > 10] <- 10
  heatmap.simple(vals, scale="none", row.clust=FALSE, col.clust=FALSE, 
                 layout.mat=mylayout, widths=mywidths, heights=myheights,
                 color.scheme=col.scheme, key.min=-2, key.max=10, title=title)
  
  if(blood.biopsy){
  idx <- which(mylayout == "rowSideLabels", arr.ind=TRUE) 
  pushViewport(viewport(name="rowSideLabels", layout.pos.row=unique(idx[,1]), 
  layout.pos.col=unique(idx[,2])))
  
  grid.text("blood", x=0.95, y=0.55,default.units="native", rot=-90,
            gp=gpar(cex=0.2 * heatmap.labels.cex("blood")))
  upViewport()
  
  idx <- which(mylayout == "colSideLabels", arr.ind=TRUE)
  pushViewport(viewport(name="colSideLabels",
                        layout.pos.row=unique(idx[,1]),
                        layout.pos.col=unique(idx[,2])))
  
  grid.text("biopsy", y=0.05, default.units="native", 
            gp=gpar(cex=0.2 * heatmap.labels.cex("biopsy")))
  upViewport()
  } 
}

plot.modules.correlation <- function(moduleColors, exprs, power = 6, max.size = 0, min.size = 0, dir="", filename="") {
	allowWGCNAThreads()
	
	if (max.size == 0 && min.size == 0) {
		colors = "all"
	} else if (max.size == 0) {
		max.size =sum(table(moduleColors))
		colors = paste("contains between",min.size,"and",max.size,"genes")
	} else {
		colors = paste("contains between",min.size,"and",max.size,"genes")
	}
	
	if (max.size > 0) {
		selected.colors.min = names(sort(table(moduleColors)))[which(sort(table(moduleColors)) > min.size)]
		selected.colors.max = names(sort(table(moduleColors)))[which(sort(table(moduleColors)) < max.size)]
		selected.colors = selected.colors.max[selected.colors.max %in% selected.colors.min]
	}
	
	ordered.indices = NULL
	ordered.expressions = NULL
	ordered.colors = NULL
	
	for (color in selected.colors) {
		print(color)
		if(color != 'grey') {
			select = sort(which(moduleColors == color))
			ordered.expressions = rbind(ordered.expressions,exprs[select,])
			
			dissTOM = 1-TOMsimilarityFromExpr(t(exprs[select,]), power = power);
			selectTree = flashClust(as.dist(dissTOM), method = "average")
			selectColor = moduleColors[select]
			
			ordered.indices = c(ordered.indices,(selectTree$order+length(ordered.indices)))
			ordered.colors = c(ordered.colors,selectColor)
		}
	}
	print('TOM disstance for the whole matrix')
	
	# calculated during module detection, but let us do it again here.
	dissTOM = 1-TOMsimilarityFromExpr(t(ordered.expressions), power = power);
	# set the diagonal of the dissimilarity to NA and raised it to the power of 4 to bring out the module structure
	#(these changes effectively amount to a change in the color scale of the plot). 
  plotDiss = dissTOM^6
	diag(plotDiss) = NA;
	png(filename=paste(dir,filename,' modules (', colors ,').png',sep=""), width=5000, height=5000, units='px')
	heatmap(plotDiss[rev(ordered.indices),ordered.indices], RowSideColors=rev(ordered.colors), ColSideColors=ordered.colors, main = paste(filename,"modules heatmap plot, (", colors ,")"), Rowv = NA, Colv = NA, scale='none', labRow = F, labCol = F)
	dev.off()
	png(filename=paste(dir,filename,' modules (', colors ,') - small.png',sep=""), width=500, height=500, units='px')
	heatmap(plotDiss[rev(ordered.indices),ordered.indices], RowSideColors=rev(ordered.colors), ColSideColors=ordered.colors, Rowv = NA, Colv = NA, scale='none', labRow = F, labCol = F)
	dev.off()
	
	gc()
}

