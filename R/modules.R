load.modules <- function(dat, mod)
{
  mymodules$biopsy$exprs <- dat$biopsy$matched.tumor$exprs
  mymodules$blood$exprs <- dat$blood$matched.tumor$exprs
  mymodules$biopsy$clinical <- dat$biopsy$matched.tumor$clinical.heatmap
  mymodules$blood$clinical <- dat$blood$matched.tumor$clinical.heatmap
  mymodules$biopsy$modules <- mod$biopsy$modules.GeneList
  mymodules$blood$modules <- mod$blood$modules.GeneList
  ## make sure the grey module always comes first
  biopsy.grey.idx <- which(names(mymodules$biopsy$modules) == "grey")
  blood.grey.idx <- which(names(mymodules$blood$modules) == "grey")
  mymodules$biopsy$modules <- c(mymodules$biopsy$modules[biopsy.grey.idx], mymodules$biopsy$modules[-biopsy.grey.idx])
  mymodules$blood$modules <- c(mymodules$blood$modules[blood.grey.idx], mymodules$blood$modules[-blood.grey.idx])
  ## remove periods and '_' from clinical names
  names(mymodules$biopsy$clinical) <- sub("_", " ", sub("\\.", " ", names(mymodules$biopsy$clinical)))
  names(mymodules$blood$clinical) <- sub("_", " ", sub("\\.", " ", names(mymodules$blood$clinical)))
  
  return(mymodules)
}

gene.overlap.test <- function(modules)
{
  all.genes <- intersect(unlist(modules$biopsy$modules[-1]), unlist(modules$blood$modules[-1]))
  
  pvals <- sapply(modules$blood$modules[-1], function(blood.mod) {
    sapply(modules$biopsy$modules[-1], function(biopsy.mod) {
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


create.modules.heatmap <- function(bs, clinical, title=title) 
{
  stopifnot(ncol(bs$dat) == nrow(clinical))  
  
  #heatmap variables
  col.clust = FALSE
  layout.m = matrix(c("key","title","","",""  ,  "","","","",""  ,"ranksum.text","ranksum.line","","",""  ,  "row.labels.rjust","heatmap","","",""  ,  "","col.labels.rjust","","",""  ,  "","","","",""  ,  "","","","",""  ,  "clinical.labels.rjust","clinical","","",""  ,  "","","","",""),nrow=9,ncol=5,byrow=TRUE)
  layout.m.sum = matrix(c("","","","",""  ,  "","","","",""  ,"","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "row.labels.rjust","heatmap","","",""  ,  "","","","",""  ,   "","","","",""  ,  "","","","",""),nrow=9,ncol=5,byrow=TRUE)
  #layout.m.dist.cdf = matrix(c("","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","col.labels.ljust","","",""  ,  "row.labels.rjust","heatmap","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""),nrow=11,ncol=5,byrow=TRUE)
  layout.m.updn = matrix(c("","","","",""  ,  "","","","",""  ,"","","","",""  ,  "","","","heatmap",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""  ,  "","","","",""),nrow=9,ncol=5,byrow=TRUE)
  widths = c(2,5,0.25,0.25,0.25)
  heights = c(1,0.25,0.5,5,0.5,0.5,0.25,2,0.25)
  
  scale = "none"
  min.val=-5
  max.val=5
  key.min=-5
  key.max=5
  
  module.info = bs$gene.order
  data = bs$dat
  mclinical = clinical[bs$pat.order,]
  
  #plot.new()
  ddrs = heatmap.simple(data,layout.mat = layout.m, widths = widths, heights = heights, col.clust = FALSE, 
                        clinical = mclinical, row.clust = FALSE, title=title,row.labels=rownames(data))
  
  ranks = bs$pat.order
  names(ranks) = 1:length(bs$pat.order)
  ranks = as.integer(names(sort(ranks)))
  ranksum = t(as.matrix(ranks))[,bs$pat.order,drop=FALSE]
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
  ranksum.plot = bs$ranksum[bs$pat.order]
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
