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

match.datasets <- function(ds1, ds2, include.normals=TRUE, remove.outliers=FALSE) {
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

  ds2.exprs = ds2$exprs
  ds2.clinical = ds2$clinical

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

find.modules <- function (exprs, power = 6, deepSplit = 2, saveTOMs=FALSE,
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
  gene.modules = blockwiseModules(t(exprs), randomSeed = 234567, power = power, deepSplit = deepSplit, minModuleSize = minModuleSize,
                                  reassignThreshold = reassignThreshold, mergeCutHeight = mergeCutHeight,
                                  numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = saveTOMs, verbose = 3)

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
  mymodules<-NULL
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
  all.genes <- intersect(unlist(modules$biopsy$modules), unlist(modules$blood$modules))

  pvals <- sapply(modules$biopsy$modules[-1], function(biopsy.mod) {
    sapply(modules$blood$modules[-1], function(blood.mod) {
      s <- length(intersect(blood.mod, all.genes))
      e <- length(intersect(biopsy.mod, all.genes))
      com <- length(intersect(biopsy.mod, blood.mod))
      #hyper.test(s, length(all.genes) - s, e, com)
      ## new computation of hypergeometri test so it matches fisher exact test
      ## ref http://rpackages.ianhowson.com/bioc/GeneOverlap/man/GeneOverlap.html
      sum(dhyper(com:e,s,length(all.genes)-s, e))
      })
  })
  ret <- p.adjust(pvals, method="BH")
  dim(ret) <- dim(pvals)
  dimnames(ret) <- dimnames(pvals)
  return(ret)
}

hyper.fisher <- function(pop1, pop2, bckg)
{
      s <- length(intersect(pop1, bckg))
      e <- length(intersect(pop2,bckg))
      com <- length(intersect(pop1, pop2))
      #hyper.test(s, length(all.genes) - s, e, com)
      ## new computation of hypergeometri test so it matches fisher exact test
      ## ref http://rpackages.ianhowson.com/bioc/GeneOverlap/man/GeneOverlap.html
      ret<- sum(dhyper(com:e,s,length(bckg)-s, e))
  return(ret)
}


create.modules.heatmap <- function(bs, dat, tissue, cohort.name="all", mod, cl.var, huc.clinical=FALSE, cl.height=8, title=title,
                                   patient.ids=NULL)
{

  #heatmap variables
  col.clust = FALSE
  layout.m = matrix(c("key","title","","",""  ,  "","","","",""  ,"ranksum.text","ranksum.line","","",""  ,  "row.labels.rjust","heatmap","","",""  ,  "","col.labels.rjust","","",""  ,  "","","","",""  ,  "","","","",""  ,  "clinical.labels.rjust","clinical","","",""  ,  "","","","",""),nrow=9,ncol=5,byrow=TRUE)
  layout.m.sum = matrix(c("","","","",""  ,  "","","","",""  ,"","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "row.labels.rjust","heatmap","","",""  ,  "","","","",""  ,   "","","","",""  ,  "","","","",""),nrow=9,ncol=5,byrow=TRUE)
  layout.m.updn = matrix(c("","","","",""  ,  "","","","",""  ,"","","","",""  ,  "","","","heatmap",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""),nrow=9,ncol=5,byrow=TRUE)
  widths = c(2,5,0.25,0.25,0.25)
  heights = c(1,0.25,0.5,3,0.25,0.25,0.25,cl.height,0.25)

  scale = "none"
  min.val=-5
  max.val=5
  key.min=-5
  key.max=5

  if (is.null(patient.ids)){
  patients <- pat.cohorts(dat[[tissue]])[[cohort.name]]
}
  else {patients <- patient.ids}

  sub.bs <- bs[[tissue]][[mod]][[cohort.name]]

  data = sub.bs$dat
  mclinical = dat[[tissue]]$clinical[patients,][sub.bs$pat.order,cl.var]

  if(!huc.clinical){
  ddrs = heatmap.simple(data,layout.mat = layout.m, widths = widths, heights = heights, col.clust = FALSE,
                        clinical = huc.color.clinical(mclinical), row.clust = FALSE, title=title,row.labels=rownames(data),
                        col.labels=rep("", ncol(data)))
  }

  if(huc.clinical){
    ddrs = heatmap.simple(data,layout.mat = layout.m, widths = widths, heights = heights, col.clust = FALSE,
                          clinical = mclinical, row.clust = FALSE, title=title,row.labels=rownames(data),
                          col.labels=rep("", ncol(data)))
  }
  ranks = sub.bs$pat.order
  names(ranks) = 1:length(sub.bs$pat.order)
  ranks = as.integer(names(sort(ranks)))
  ranksum = t(as.matrix(ranks))[,sub.bs$pat.order,drop=FALSE]
  rownames(ranksum) = "ranks"
  heatmap.simple(-ranksum, layout.mat = layout.m.sum, widths = widths, heights = heights, col.clust = FALSE, row.clust = FALSE)

  up.dn = as.vector(array(1,dim=c(1,length(sub.bs$gene.order))))
  names(up.dn) = rownames(data)
  if(length(sub.bs$dn)>0){
  up.dn[names(up.dn) %in% rownames(dat[[tissue]]$exprs)[sub.bs$dn]] = -1}
  to.plot = (as.matrix(up.dn,ncol=1))
  color.scheme = heatmap.color.scheme(low.breaks=c(-1.5,0),high.breaks=c(0,1.5))
  heatmap.simple(to.plot, scale=scale, layout.mat = layout.m.updn, widths = widths, heights = heights, col.clust = FALSE, row.clust = FALSE, color.scheme = color.scheme)

  the.layout <- grid.layout(nrow(layout.m), ncol(layout.m), widths=widths, heights=heights)
  top.vp <- viewport(layout=the.layout, name="heatmap.top.vp")
  pushViewport(top.vp)
  elem = 'ranksum.line'

  ranksum.plot = sub.bs$ranksum[sub.bs$pat.order]
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

  first.ind = length(which(sub.bs$roi.cat==3))+1
  last.ind = first.ind + length(which(sub.bs$roi.cat==2))

  res.random.dist.begin = first.ind
  res.random.dist.end = last.ind

  elem = 'heatmap'
  idx <- which(layout.m == elem, arr.ind=TRUE)
  pushViewport(viewport(name=elem,
                        layout.pos.row=unique(idx[,1]),
                        layout.pos.col=unique(idx[,2])))

  par(new=TRUE, fig=gridFIG(), mar=c(0,0,0,0))
  grid.lines(x = unit(c(first.ind/length(sub.bs$roi),first.ind/length(sub.bs$roi)), "npc"),gp=gpar(col='yellow',lty = 1, lwd = 2))
  grid.lines(x = unit(c(last.ind/length(sub.bs$roi),last.ind/length(sub.bs$roi)), "npc"),gp=gpar(col='yellow',lty = 1, lwd = 2))
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
  plotDiss = dissTOM^7
	diag(plotDiss) = NA;
	png(filename=paste(dir,filename,' modules (', colors ,').png',sep=""), width=5000, height=5000, units='px')
	heatmap(plotDiss[rev(ordered.indices),ordered.indices], RowSideColors=rev(ordered.colors), ColSideColors=ordered.colors, main = paste(filename,"modules heatmap plot, (", colors ,")"), Rowv = NA, Colv = NA, scale='none', labRow = F, labCol = F)
	dev.off()
	png(filename=paste(dir,filename,' modules (', colors ,') - small.png',sep=""), width=500, height=500, units='px')
	heatmap(plotDiss[rev(ordered.indices),ordered.indices], RowSideColors=rev(ordered.colors), ColSideColors=ordered.colors, Rowv = NA, Colv = NA, scale='none', labRow = F, labCol = F)
	dev.off()

	gc()
}

parallelset <- function(..., freq, col="gray", border=0, layer,
                        alpha=0.5, gap.width=0.05) {
  p <- data.frame(..., freq, col, border, alpha, stringsAsFactors=FALSE)
  n <- nrow(p)
  if(missing(layer)) { layer <- 1:n }
  p$layer <- layer
  np <- ncol(p) - 5
  d <- p[ , 1:np, drop=FALSE]
  p <- p[ , -c(1:np), drop=FALSE]
  p$freq <- with(p, freq/sum(freq))
  col <- col2rgb(p$col, alpha=TRUE)
  if(!identical(alpha, FALSE)) { col["alpha", ] <- p$alpha*256 }
  p$col <- apply(col, 2, function(x) do.call(rgb, c(as.list(x), maxColorValue = 256)))
  getp <- function(i, d, f, w=gap.width) {
    a <- c(i, (1:ncol(d))[-i])
    o <- do.call(order, d[a])
    x <- c(0, cumsum(f[o])) * (1-w)
    x <- cbind(x[-length(x)], x[-1])
    gap <- cumsum( c(0L, diff(as.numeric(d[o,i])) != 0) )
    gap <- gap / max(gap) * w
    (x + gap)[order(o),]
  }
  dd <- lapply(seq_along(d), getp, d=d, f=p$freq)
  par(mar = c(0, 0, 2, 0) + 0.1, xpd=TRUE )
  plot(NULL, type="n",xlim=c(0, 1), ylim=c(np, 1),
       xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlab='', ylab='', frame=FALSE)
  for(i in rev(order(p$layer)) ) {
    for(j in 1:(np-1) )
      polygon(c(dd[[j]][i,], rev(dd[[j+1]][i,])), c(j, j, j+1, j+1),
              col=p$col[i], border=p$border[i])
  }
  text(0, seq_along(dd), labels=names(d), adj=c(0,-2), font=2)
  for(j in seq_along(dd)) {
    ax <- lapply(split(dd[[j]], d[,j]), range)
    for(k in seq_along(ax)) {
      lines(ax[[k]], c(j, j))
      text(ax[[k]][1], j, labels=names(ax)[k], adj=c(0, -0.25))
    }
  }
}

### pat.cohorts()
###
### returns a list containing information about each cohort for the
### given dataset.
###
### Arguments:
###
###     dat:
###
###         rearranged dataset object
###
### Returns:
###
###     a list with each element corresponding to a cohort.
###
###         a vector of (integer) indices pointing out the patients that belong to the cohort
###

pat.cohorts <- function(dat)
{
  ret <- list()
  cl <- dat$clinical

  ## define all the patient indices for all our cohorts
  ret$all <- rownames(cl)
  ret$erp <- rownames(cl)[which(cl$er | is.na(cl$er))]
  ret$ern <- rownames(cl)[which(!cl$er| is.na(cl$er))]
  ret$her2p <- rownames(cl)[which(cl$her2| is.na(cl$er))]
  ret$her2n <- rownames(cl)[which(!cl$her2| is.na(cl$er))]
  ret$erp.her2p <- rownames(cl)[which(cl$er & cl$her2 | is.na(cl$er))]
  ret$ern.her2p <- rownames(cl)[which(!cl$er & cl$her2 | is.na(cl$er))]
  ret$erp.her2n <- rownames(cl)[which(cl$er & !cl$her2 | is.na(cl$er))]
  ret$ern.her2n <- rownames(cl)[which(!cl$er & !cl$her2| is.na(cl$er))]
  ret$luma <- rownames(cl)[which(cl$pam50.parker == "LumA" | is.na(cl$er))]
  ret$erp.luma <- rownames(cl)[which(cl$pam50.parker == "LumA" & cl$er | is.na(cl$er))]
  ret$lumb <- rownames(cl)[which(cl$pam50.parker == "LumB"| is.na(cl$er))]
  ret$erp.lumb <- rownames(cl)[which(cl$pam50.parker == "LumB" & cl$er | is.na(cl$er))]
  ret$normL <- rownames(cl)[which(cl$pam50.parker == "Normal" | is.na(cl$er))]
  ret$erp.normL <- rownames(cl)[which(cl$pam50.parker == "Normal" & cl$er | is.na(cl$er))]
  ret$basalL <- rownames(cl)[which(cl$pam50.parker == "Basal"| is.na(cl$er))]
  ret$erp.basalL <- rownames(cl)[which(cl$pam50.parker == "Basal" & cl$er| is.na(cl$er))]
  ret$her2E <- rownames(cl)[which(cl$pam50.parker == "Her2"| is.na(cl$er))]
  ret$erp.her2E <- rownames(cl)[which(cl$pam50.parker == "Her2" & cl$er| is.na(cl$er))]
  ret$cit.luma <- rownames(cl)[which(cl$cit == "lumA"| is.na(cl$er))]
  ret$cit.lumb <- rownames(cl)[which(cl$cit == "lumB"| is.na(cl$er))]
  ret$cit.lumc <- rownames(cl)[which(cl$cit == "lumC"| is.na(cl$er))]
  ret$cit.normL <- rownames(cl)[which(cl$cit == "normL"| is.na(cl$er))]
  ret$cit.mApo <- rownames(cl)[which(cl$cit == "mApo"| is.na(cl$er))]
  ret$cit.basalL <- rownames(cl)[which(cl$cit == "basL"| is.na(cl$er))]
  ret$intclust1 <- rownames(cl)[which(cl$intclust == "intclust1"| is.na(cl$er))]
  ret$intclust3 <- rownames(cl)[which(cl$intclust == "intclust3"| is.na(cl$er))]
  ret$intclust4 <- rownames(cl)[which(cl$intclust %in% c("intclust4.ern", "intclust4.erp")| is.na(cl$er))]
  ret$intclust4.erp <- rownames(cl)[which(cl$intclust == "intclust4.erp"| is.na(cl$er))]
  ret$intclust4.ern <- rownames(cl)[which(cl$intclust == "intclust4.ern"| is.na(cl$er))]
  ret$intclust5 <- rownames(cl)[which(cl$intclust == "intclust5"| is.na(cl$er))]
  ret$intclust7 <- rownames(cl)[which(cl$intclust == "intclust7"| is.na(cl$er))]
  ret$intclust8 <- rownames(cl)[which(cl$intclust == "intclust8"| is.na(cl$er))]
  ret$intclust9 <- rownames(cl)[which(cl$intclust == "intclust9"| is.na(cl$er))]
  ret$intclust10 <- rownames(cl)[which(cl$intclust == "intclust10"| is.na(cl$er))]

return(ret)
}

### plot.pat.bs()
###
### plot bresat heatmap for the specified data set and
### cohort. Genes can also be specified.
###
### Arguments:
###     bs: bresat object (list) returned by sig.ranksum for data set 'dat'
###
###     cl: data.frame 'clinical'
###
###     pat.cohort: list of patient cohorts
###
###     cohort.name
###
###         name of cohort
###
###     gene.names
###
###         character vector with gene names. If this is missing or is
###         NULL, then all genes will be used. See 'rm.missing.genes'
###         for handling of gene names that are not found in the data
###         set.
###
###    title
###



plot.pat.bs <- function(bs, dat, cohorts, cohort.name="all", patient.ids=NULL, gene.names=NULL, blood.mod, biopsy.mod, tissue.order.by="biopsy", title, reorder.genes=F)
{
# Plot layouts -----------------
  ## layout heatmap top left corner
  layout.m = matrix(c("key","title","","","","",
                      "","","","","","",
                      "ranksum.text","ranksum.line","","","","",
                      "row.labels.rjust","heatmap","","","" ,"boxplot",
                      "","","","","","",
                      "ranks.text","ranks","","","","",
                      "","","","","","",
                      "key.2","title.2","","","","",
                      "test","","","","","",
                      "ranksum.text.2","ranksum.line.2","","","","",
                      "row.labels.rjust2","heatmap.2","","","" ,"scatterplot",
                      "","col.labels.rjust","","",""  ,"",
                      "ranks.text.2","ranks.2","","","","",
                      "","","","","","",
                      "clinical.labels.rjust","clinical","","","","",
                      "","","","","", ""),nrow=16,ncol=6,byrow=TRUE)
  ## layout heatmap bottom left corner
  layout.m.2 = matrix(c("","","","","","",
                        "","","","","","",
                        "","","","","","",
                        "","","","","","",
                        "","","","","","",
                        "","","","","","",
                        "","","","","","",
                      "key","title","","","","",
                      "","","","","","",
                      "ranksum.text","ranksum.line","","","","",
                      "row.labels.rjust","heatmap","","","" ,"",
                      "","col.labels.rjust","","",""  ,"",
                      "ranks.text","ranks","","","","",
                      "","","","","","",
                      "clinical.labels.rjust","clinical","","","","",
                      "","","","","", ""),nrow=16,ncol=6,byrow=TRUE)
  ## layout updn heatmap for top left corner
  layout.m.updn = matrix(c("","","","","","",
                           "","","","","" ,"",
                           "","","","","","",
                           "","","","heatmap","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","",""),nrow=16,ncol=6,byrow=TRUE)
  ## layout updn heatmap for bottom left corner
  layout.m.updn.2 = matrix(c("","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                           "","","","","","",
                         "","","","heatmap","","",
                         "","","","","","",
                         "","","","","","",
                         "","","","","","",
                         "","","","","","",
                         "","","","","",""),nrow=16,ncol=6,byrow=TRUE)

  ## layout dimensions
  widths = c(2,5,0.25,0.25,0.25,8)
  heights = c(1,0.25,0.5,3,0.25,0.25,0.25,1,0.25,0.5,3,0.25,0.25,0.25,5,0.25)

  ## key
  scale = "none"
  min.val=-5
  max.val=5
  key.min=-5
  key.max=5

# define var ------------------
  ## define reordered variables
  if(tissue.order.by=="blood"){
  bs.order.by <- bs$blood[[blood.mod]][[cohort.name]]}
  else {bs.order.by <- bs$biopsy[[biopsy.mod]][[cohort.name]]}

  order.by<-bs.order.by$pat.order
  roi<-bs.order.by$roi
  roi.cat<-bs.order.by$roi.cat

  ## define patients to include
  if (is.null(patient.ids))
  {
    patients<-cohorts[[cohort.name]]
  }
  else
  {
    patients<-patient.ids
  }

# define clinical and exprs ------------------
  ## define reordered clinical data
  cl<-dat$blood$clinical[rownames(dat$blood$clinical) %in% patients,]
  mclinical = cl[order.by,]
  bnclinical=dat$bnblood$clinical[rownames(dat$bnblood$clinical) %in% patients, ]


  ## define blood reordered expression and select genes
  ss.bs.blood <- bs$blood[[blood.mod]][[cohort.name]]
  blood.data = ss.bs.blood$dat[, match(rownames(mclinical), colnames(ss.bs.blood$dat))]
  if (!is.null(gene.names))
  {blood.data<-blood.data[rownames(blood.data) %in% gene.names,]}
  if(reorder.genes == T ){
    genes <- rownames(dat$blood$exprs)[bs.order.by$gene.order]
    blood.data<-blood.data[match(genes, rownames(blood.data)),]
  }


  ## define biopsy reordered expression data and select genes
  ss.bs.biopsy <- bs$biopsy[[biopsy.mod]][[cohort.name]]
  biopsy.data = ss.bs.biopsy$dat[, match(rownames(mclinical), colnames(ss.bs.biopsy$dat))]
  if (!is.null(gene.names))
  {biopsy.data<-biopsy.data[rownames(biopsy.data) %in% gene.names,]}
  if(reorder.genes == T ){
    genes <- rownames(dat$biopsy$exprs)[bs.order.by$gene.order]
    biopsy.data<-biopsy.data[match(genes, rownames(biopsy.data)),]
  }


# plot biopsy heatmap ----------------------
  ### Biopsy heatmap top left (without clinical)
  ddrs = heatmap.simple(biopsy.data,
                        layout.mat = layout.m, widths = widths, heights = heights, col.clust = FALSE,
                        row.clust = FALSE, title=paste(cohort.name, biopsy.mod, "biopsy ordered by", biopsy.mod, "biopsy"),
                        row.labels=rownames(biopsy.data),
                        col.labels=rep("", length(colnames(biopsy.data))))

  up.dn = as.vector(array(1,dim=c(1,length(ss.bs.biopsy$gene.order))))
  names(up.dn) = rownames(ss.bs.biopsy$dat)
  if(length(ss.bs.biopsy$dn)>0){
    up.dn[names(up.dn) %in% rownames(dat$biopsy$exprs)[ss.bs.biopsy$dn]] = -1}
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

  rank.colors<-rev(diverge_hcl(n=ncol(ss.bs.biopsy$dat)))
  names(rank.colors)<-colnames(ss.bs.biopsy$dat)
  rank.colors<-rank.colors[match(rownames(cl)[order.by], colnames(ss.bs.biopsy$dat))]
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
  ranksum.plot = ss.bs.biopsy$ranksum[order.by][colnames(biopsy.data) %in% patients]
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

  ## plot roi lines bottom left
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

  # plot blood heatmap ---------------------
  ## plot blood heatmap bottom left
  bc.var<- c("er","her2" ,"pam50.parker","hybrid","cit", "claudin.low", "intclust",
             "lymph" ,"t.size",
             "menopause","hrt","medication","hospital",
             "age","weight", "MKS", "LUMS", "HER2S")

  ddrs = heatmap.simple(blood.data,
                        clinical = huc.color.clinical(mclinical)[rownames(mclinical) %in% patients,][,bc.var],
                        layout.mat = layout.m.2, widths = widths, heights = heights, col.clust = FALSE,
                        row.clust = FALSE, title=paste(cohort.name, blood.mod, "blood ordered by", biopsy.mod, "biopsy"),
                        row.labels=rownames(blood.data),
                        col.labels=rep("", length(colnames(blood.data))))

  ## plot updn top left heatmap
  up.dn = as.vector(array(1,dim=c(1,length(ss.bs.blood$gene.order))))
  names(up.dn) = rownames(ss.bs.blood$dat)
  if(length(ss.bs.blood$dn)>0){
    up.dn[names(up.dn) %in% rownames(dat$blood$exprs)[ss.bs.blood$dn]] = -1}
  to.plot = (as.matrix(up.dn,ncol=1))
  color.scheme = heatmap.color.scheme(low.breaks=c(-1.5,0),high.breaks=c(0,1.5))
  heatmap.simple(to.plot, scale=scale, layout.mat = layout.m.updn.2, widths = widths, heights = heights, col.clust = FALSE, row.clust = FALSE, color.scheme = color.scheme)

  ## plot ranks for top left heatmap
  the.layout <- grid.layout(nrow(layout.m), ncol(layout.m), widths=widths, heights=heights)
  mid.vp <- viewport(layout=the.layout, name="heatmap.mid.vp")
  pushViewport(mid.vp)
  elem = 'ranks.2'
  idx <- which(layout.m == elem, arr.ind=TRUE)
  pushViewport(viewport(name=elem,
                        layout.pos.row=unique(idx[,1]),
                        layout.pos.col=unique(idx[,2])))

  rank.colors<-rev(diverge_hcl(n=ncol(ss.bs.blood$dat)))
  names(rank.colors)<-colnames(ss.bs.blood$dat)
  rank.colors<-rank.colors[match(rownames(cl)[order.by], colnames(ss.bs.blood$dat))]
  ranksum = t(as.matrix(rank.colors))[,names(rank.colors) %in% patients,drop=FALSE]
  heatmap.clinical(ranksum)
  upViewport()

  elem = 'ranks.text.2'
  idx <- which(layout.m == elem, arr.ind=TRUE)
  pushViewport(viewport(name=elem,
                        layout.pos.row=unique(idx[,1]),
                        layout.pos.col=unique(idx[,2])))
  heatmap.labels('ranks', type="row.labels", just="right")
  upViewport()

  ## plot ranksum line top left heatmap
  the.layout <- grid.layout(nrow(layout.m), ncol(layout.m), widths=widths, heights=heights)
  top.vp <- viewport(layout=the.layout, name="heatmap.top.vp")
  pushViewport(top.vp)

  ranksum.plot = ss.bs.blood$ranksum[order.by][colnames(blood.data) %in% patients]
  elem = 'ranksum.line.2'
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

  elem = 'ranksum.text.2'
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

  elem = 'heatmap.2'
  idx <- which(layout.m == elem, arr.ind=TRUE)
  pushViewport(viewport(name=elem,
                        layout.pos.row=unique(idx[,1]),
                        layout.pos.col=unique(idx[,2])))

  par(new=TRUE, fig=gridFIG(), mar=c(0,0,0,0))
  grid.lines(x = unit(c(first.ind/length(which(rownames(cl) %in% patients)),first.ind/length(which(rownames(cl) %in% patients))), "npc"),gp=gpar(col='yellow',lty = 1, lwd = 2))
  grid.lines(x = unit(c(last.ind/length(which(rownames(cl) %in% patients)),last.ind/length(which(rownames(cl) %in% patients))), "npc"),gp=gpar(col='yellow',lty = 1, lwd = 2))
  upViewport()


# plot blood boxplot in cases and controls ------------------------
  ### Bnblood boxplot data
  ss.bs.bnblood <- bs$bnblood[[blood.mod]][[cohort.name]]
  bnblood<-data.frame(bnbl.ranksum=ss.bs.bnblood$ranksum[c(which(bnclinical$cancer==TRUE), which(bnclinical$cancer==FALSE))],
                      bl.ranksum=c(ss.bs.blood$ranksum, rep(NA, length(which(bnclinical$cancer==FALSE)))),
                      t.ranksum=c(ss.bs.biopsy$ranksum, rep(NA, length(which(bnclinical$cancer==FALSE)))),
                      subtype=c(rep(cohort.name, length(which(bnclinical$cancer==TRUE))),
                               rep("normal", length(which(bnclinical$cancer==FALSE)))))
  bnblood$roi.cat<-c(roi.cat, rep(NA, length(which(as.character(bnblood$subtype)=="normal"))))
  bnblood$cancer<-1
  bnblood$cancer<-ifelse(bnblood$roi.cat==1 & !is.na(bnblood$roi.cat), 4, as.character(bnblood$cancer))
  bnblood$cancer<-ifelse(bnblood$roi.cat==2 & !is.na(bnblood$roi.cat), 2, as.character(bnblood$cancer))
  bnblood$cancer<-ifelse(bnblood$roi.cat==3 & !is.na(bnblood$roi.cat), 1, as.character(bnblood$cancer))
  bnblood$cancer<-as.numeric(bnblood$cancer)

  bnblood$tumor.cat<-"control"
  bnblood$tumor.cat<-ifelse(bnblood$roi.cat==1 & !is.na(bnblood$roi.cat), "low", as.character(bnblood$tumor.cat))
  bnblood$tumor.cat<-ifelse(bnblood$roi.cat==2 & !is.na(bnblood$roi.cat), "mid", as.character(bnblood$tumor.cat))
  bnblood$tumor.cat<-ifelse(bnblood$roi.cat==3 & !is.na(bnblood$roi.cat), "high", as.character(bnblood$tumor.cat))
  bnblood$tumor.cat<-factor(bnblood$tumor.cat, levels=c("low", "mid", "high", "control"))
  bnblood$tumor.cat.ordered<-factor(bnblood$tumor.cat, levels=c("high", "mid", "low", "control"), ordered=T)

  sub.col <- c(normal="white", all="grey", erp="green", ern="firebrick2", her2p="hotpink2", her2n="#21B6A8",
               erp.her2p="orange", ern.her2p="hotpink2", erp.her2n="blue", ern.her2n="firebrick2",intclust10="firebrick2",
              luma="blue4", erp.luma="blue4", intclust3= "blue4", intclust7= "blue4", intclust8= "blue4",
              lumb="deepskyblue", erp.lumb="deepskyblue", intclust1= "deepskyblue", intclust9= "deepskyblue",
              normL="green4", erp.normL="green4", intclust4.erp="green4", basalL="firebrick2", 
              intclus5= "hotpink2",her2E="hotpink2", erp.her2E="orange", erp.basalL="#7fffd4",
              cit.luma="blue4", cit.lumb="deepskyblue", cit.normL="green4", cit.mApo="hotpink2", cit.lumc="#7fffd4",  cit.basalL="firebrick2")
  bnblood$sub.col<-sub.col[as.character(bnblood$subtype)]

  ### plot boxplot top right corner
  elem = 'boxplot'
  idx <- which(layout.m == elem, arr.ind=TRUE)
  vp<-viewport(name=elem,
               layout.pos.row=2:10,
               layout.pos.col=unique(idx[,2]))
  pushViewport(vp)

  p<-ggplot(data = bnblood, aes(x=tumor.cat.ordered, y=bnbl.ranksum)) +
    geom_boxplot(aes(fill=sub.col, alpha=1/cancer))+
    geom_boxplot(colour="black", outlier.shape ="+", outlier.size = 1, fill=NA)+
    geom_point(shape="+")+
    scale_fill_manual(values=levels(factor(bnblood$sub.col)))+
    labs(x="tumor module category",
         y="blood module ranksum",
         title=paste(cohort.name, " patients and controls\n(aov p=", as.character(signif(anova(lm(bnblood$bnbl.ranksum~bnblood$tumor.cat))$`Pr(>F)`[1], digits=1)), ")", sep=""))+
    theme(legend.position="none",
          panel.background = element_rect(fill = "transparent",colour = NA),
          axis.line.x = element_line(colour="grey60"),
          axis.line.y = element_line(colour="grey60"),
          axis.title.x=element_text(hjust=0.2),
          axis.title=element_text(size=10),
          plot.title = element_text(hjust=0,vjust=1, size=12)
          )
  print(p, vp=vp)
  upViewport()

  ### plot ranksum scatterplot for BC patients in blood vs tumor
  elem = 'scatterplot'
  idx <- which(layout.m == elem, arr.ind=TRUE)
  vp1<-viewport(name=elem,
               layout.pos.row=11:15,
               layout.pos.col=unique(idx[,2]))
  pushViewport(vp1)

  perm.cor.p <- get(load("../../data/mixt/perm_cor_p.RData"))
  p1<-ggplot(bnblood[bnblood$tumor.cat != "control",], aes(x=bl.ranksum,y=t.ranksum))+
    geom_smooth(method="lm", colour="white", alpha=0.2, size=0.4)+
    geom_point(aes(colour=sub.col, alpha=1/cancer), size=2)+
    geom_point(shape = 1,size = 2,colour = "black")+
    scale_colour_manual(values=levels(factor(bnblood$sub.col)))+
    labs(y="Tumor module ranksum",
         x="Blood module ranksum",
         title=paste(cohort.name, " patients", " (cor=",as.character(signif(cor.test(bnblood[bnblood$cancer != "control",]$bl.ranksum, bnblood[bnblood$cancer != "control",]$t.ranksum)$estimate, digits=1)),
                     ", p=",signif(perm.cor.p$blood.biopsy[[cohort.name]][blood.mod, biopsy.mod], digits = 3),")",sep="")) +
    theme(legend.position="none",
          panel.background = element_rect(fill = "transparent",colour = NA),
          axis.line.x   = element_line(colour="grey60"),
          axis.line.y   = element_line(colour="grey60"),
          axis.title = element_text(size=10),
          plot.title = element_text(hjust=0, vjust=1, size=12))

  print(p1, vp=vp1)
  upViewport()
  invisible(NULL)

}


sampledCorModules = function(
  tissue1,
  tissue2,
  nRuns,
  randomSeed = 12345,
  skipUnsampledCalculation = FALSE,
  mc.cores=80,
  datRank,
  corType = "p",
  ...,
  verbose = 2, indent = 0)

{
  spaces = indentSpaces(indent);

  result = list();
  nSamples = length(datRank[[tissue1]][[1]]);

  seedSaved = FALSE;
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
      seedSaved = TRUE;
      savedSeed = .Random.seed
    }
    set.seed(randomSeed);
  }

  mods <- mclapply(1:nRuns, function(i) {
    set.seed(randomSeed + 2*i + 1);
    if (verbose > 0) printFlush(paste(spaces, "...working on run", i, ".."));

    if (i > 1 || skipUnsampledCalculation)
    {
      useTissue1Samples = sample(nSamples)
      useTissue2Samples = sample(nSamples)

    } else {
      useTissue1Samples = c(1:nSamples)
      useTissue2Samples = c(1:nSamples)}

    mods <- laply(datRank[[tissue1]], function (x){
      laply(datRank[[tissue2]], function (y){
        cor(x[useTissue1Samples],y[useTissue2Samples], use=corType)
      })
    })

    rownames(mods) <- names(datRank[[tissue1]])
    colnames(mods) <- names(datRank[[tissue2]])

    return(mods)
  }, mc.cores = mc.cores)
}

