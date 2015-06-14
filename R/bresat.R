if(R.version$major >= 3){
  require(parallel)
} else {
  require(multicore)
}

###
### file: bresat.R
### author: Ali Tofigh
###
### Contains functions related to Robert Lesurf's BreSAT ideas and
### signature database sigdb.


### sig.ranksum()
###
### Robert Lesurf's main BreSAT function that ranks the samples given a
### signature. In this implementation, patients with the same ranksum
### are assigned their "average" rank, i.e., they all receive the same
### rank value, which might not be an integer. One benefit of this
### strategy is that exchanging the up- and down-genes reverses the
### ranks, and ordering the patients in decreasing order (with a stable
### sort algorithm) will be the same as the order of the patients in
### increasing order before the exchange. This consistent behavior in
### turn makes sig.correlation behave consistently when computing
### correlations between patient ranks with respect to different gene
### signatures.
###
### Note that 'up' and 'dn' are used to index into exprdata. Hence, these can
### be either numerical indices or row names. NAs in exprdata are not
### supported.
###
### Arguments:
###     exprdata        the gene expression data, must be a matrix
###     up              indices of up-regulated genes
###     dn              indices of down-regulated genes
###     ns              indices of genes whose directions have not been specified.
###                     These genes will be clustered into two parts based on
###                     their correlation to each other and up and dn indices
###                     will be assigned by the function. It is currently not
###                     possible to mix 'up'/'dn' with 'ns'.
###     full.return     returns a struct when true, otherwise only the rank
###
### Returns:
###     if 'full.return' is FALSE, then the computed rank is
###     returned. Otherwise, a list is returned with the following members.
###
###     rank            rank of each patient, higher = better
###     dat             exprdata for the signature genes only and ordered both by
###                     patient and gene.
###     up.dn           vector with nrow(dat) values. -1, indicates down-gene
###                     and 1 indicates up-gene
###     pat.order       patient ordering, with best patient first.
###                     exprdata[, pat.order] orders the columns according to the
###                     strength of the signature.
###     gene.order      exprdata[gene.order, ] sorts rows of exprdata so that all
###                     up-genes are first and all down genes last.  Also, the
###                     most correlated up-gene is at top and most correlated
###                     down-gene is at the bottom.
sig.ranksum <- function(exprdata, up=NULL, dn=NULL, ns=NULL, full.return = FALSE, ranks.matrix=FALSE)
{
  if(ranks.matrix && full.return)
    stop("Cannot return full results of bresat ranksum with ranks as input")
  if(ranks.matrix && length(ns) > 0)
    stop("Cannot compute bresat ranksum with ranks as input on a signature without identified up and dn genes")
  
  if (length(dim(exprdata)) != 2)
        stop("'exprdata' must be a matrix")
    if (length(up) == 0 && length(dn) == 0 && length(ns) == 0)
        stop("no indices were specified")
    if (length(ns) > 0 && (length(up) > 0 || length(dn) > 0))
        stop("both directional and non-directional indices were specified")
    if (is.logical(up))
        up <- which(up)
    if (is.logical(dn))
        dn <- which(dn)
    if (is.logical(ns))
        ns <- which(ns)

    if (ncol(exprdata) < 2) {
      ret <- list()
      up <- c(up,ns)
      ret$rank <- 1
      ret$pat.order <- 1
      ret$gene.order <- c(up,dn)
      ret$dat <- exprdata[ret$gene.order,, drop=FALSE]
      ret$up.dn <- c(rep(1, length(up)), rep(-1, length(dn)))
      return(ret)
    }

    if (length(ns) > 0)
    {
        ## identify rows in exprdata that have zero variance
        library(matrixStats)
        tmp <- rowSds(exprdata[ns,,drop=F],na.rm=TRUE) == 0 ## faster
        #tmp <- sapply(ns, function(idx) {sd(exprdata[idx, ], na.rm=TRUE)}) == 0
        zero.sd.idx <- ns[tmp]
        ns <- ns[!tmp]

        if (length(ns) == 1)
        {
            up <- ns
        }
        else if (length(ns) == 2)
        {
            if (cor(exprdata[ns[1], ], exprdata[ns[2], ], use="pairwise") < 0)
            {
                up <- ns[1]
                dn <- ns[2]
            }
            else
            {
                up <- ns
            }
        }
        else if (length(ns) > 2)
        {
            library(cluster)
            diss <- 1-cor(t(exprdata[ns, , drop=FALSE]), use="pairwise")
            diss[which(is.na(diss))] <- 1
            diss[diss >= 1] <- diss[diss >= 1] + 1

            clustering <- pam(diss, k=2, diss=TRUE, cluster.only=TRUE)
            up.cluster <- which.max(table(clustering))
            up <- ns[which(clustering == up.cluster)]
            dn <- ns[which(clustering != up.cluster)]
            up.dn.cor <- cor(t(exprdata[up, , drop=FALSE]), t(exprdata[dn, , drop=FALSE]), use="pairwise")
            if (sum(up.dn.cor < 0,na.rm=T) < length(up) * length(dn) / 2)
            {
                up <- ns
                dn <- NULL
            }
        }

        if (length(zero.sd.idx) > 0)
        {
            up <- c(up, zero.sd.idx)
        }
    }

    ranksum <- double(ncol(exprdata))
    col.counts <- rep(0, ncol(exprdata))
    if (length(up) != 0)
    {
        dat <- exprdata[up, , drop = FALSE]
        if(ranks.matrix)
          ranksum <- colSums(dat, na.rm=TRUE)
        else
        ranksum <- rowSums(apply(dat, 1, function(x) {.Internal(rank(x, length(x), "average"))}))
        col.counts <- colSums(!is.na(dat))
    }
    if (length(dn) != 0)
    {
        dat <- exprdata[dn, , drop = FALSE]
        if(ranks.matrix)
          ranksum <- ranksum + colSums(ncol(exprdata) - dat + 1, na.rm=TRUE)
        else
        ranksum <- ranksum + rowSums(ncol(exprdata) - apply(dat, 1, function(x) {.Internal(rank(x, length(x), "average"))}) + 1)
        col.counts <- col.counts + colSums(!is.na(dat))
    }

  ranksum <- ranksum / col.counts  
  rank <- .Internal(rank(ranksum, length(ranksum),"average"))

    if (full.return == FALSE)
        return(rank)

    if (length(up) == 0)
        up <- NULL
    if (length(dn) == 0)
        dn <- NULL

    ## computes correlation of a single row of the expression data with the
    ## patient ordering. Used to obtain ret$gene.order.
    gene.cor <- function(gene.idx, is.up, exprdata, pat.order)
    {
        gene.expr <- exprdata[gene.idx, pat.order]
        if (is.up == TRUE)
            gene.expr <- -gene.expr
        cor(rank(gene.expr), 1:ncol(exprdata))
    }

    ret <- list()
    ret$rank <- rank
    ## best first, patient nr pat.order[1] is the best patient for this sig so
    ## that dat[ , pat.order] sorts the patients correctly with nr 1 = best.
    ret$pat.order <- order(-rank)
    ## Also, gene.order orders the genes so that the up-genes come first,
    ## followed by the down-genes. The "best" up-gene first and the "best"
    ## down-gene last.
    up.cor <- sapply(up, gene.cor, TRUE, exprdata, ret$pat.order)
    dn.cor <- sapply(dn, gene.cor, FALSE, exprdata, ret$pat.order)
    ret$gene.order <- c(up[rev(order(up.cor))], dn[order(dn.cor)])
    ret$up <- up
    ret$dn <- dn
    ret$up.cor<-up.cor
    ret$dn.cor<-dn.cor
    ret$dat <- exprdata[ret$gene.order, ret$pat.order, drop=FALSE]
    ret$up.dn <- c(rep(1, length(up)), rep(-1, length(dn)))
    ret$ranksum <- ranksum

    return(ret)
}


### sig.correlation()
###
### Uses BreSAT patient ordering to obtain Pearson correlation between two
### signatures in epi and stroma.
###
### Arguments:
###     epi.exprdata:   the epithelial exprdata
###     str.exprdata:   the stromal exprdata
###     epi.up:         indices of up-regulated genes in epi-signature
###     epi.dn:         indices of down-regulated genes in epi-signature
###     str.up:         indices of up-regulated genes in stroma-signature
###     str.dn:         indices of down-regulated genes in stroma-signature
###     full.return:    returns a struct when TRUE, otherwise only the
###                     correlation
###
###     Note that epi.up, epi.dn, str.up, str.dn are used to index into
###     epi.exprdata and str.exprdata. Hence, these can be either numerical
###     indices or row names. NAs are not supported.
###
### Returns:
###     if 'full.return' is FALSE, then the correlation between signatures in
###     epi and stroma. Otherwise, a list is returned with the following
###     members:
###
###     cor                 correlation of the signatures
###     epi$rank            see below for the rest
###     epi$dat
###     epi$up.dn
###     epi$pat.order
###     epi$gene.order
###     str$rank
###     str$dat
###     str$up.dn
###     str$pat.order
###     str$gene.order
###
###     the returned 'epi' and 'str' lists are exactly the structs returned by
###     sig.ranksum(), so see the comments there for more info.
sig.correlation <- function(epi.exprdata, str.exprdata,
                            epi.up, epi.dn, str.up, str.dn,
                            full.return = FALSE)
{
    epi.rank <- sig.ranksum(epi.exprdata, epi.up, epi.dn, full.return)
    str.rank <- sig.ranksum(str.exprdata, str.up, str.dn, full.return)

    if (full.return == FALSE)
    {
        ret <- cor(epi.rank, str.rank)
    }
    else
    {
        ret <- list()
        ret$cor <- cor(epi.rank$rank, str.rank$rank)
        ret$epi <- epi.rank
        ret$str <- str.rank
    }

    return(ret)
}

### Ali's original function
# random.ranks <- function(bs, n=1000)
# {
#   datrank.up <- t(apply(bs$dat[bs$up.dn > 0, , drop=FALSE], 1, function(x) {rank(x, "average", na.last="keep")}))
#   datrank.dn <- ncol(bs$dat) - t(apply(bs$dat[bs$up.dn < 0, , drop=FALSE], 1, function(x) {rank(x, "average", na.last="keep")})) + 1
#   datrank <- rbind(datrank.up, datrank.dn)
#   nvals.cols <- nrow(datrank) - c(colSums(is.na(datrank)), 0)
#   nvals.rows <- ncol(datrank) - rowSums(is.na(datrank))
#   random.cols <- t(sapply(1:nrow(datrank), function(i) {runif(n, 0, nvals.rows[i] + 1)}))
#   rand.dist <- sapply(1:n, function(i) {
#     datrank.rand.col <- cbind(datrank, random.cols[, i])
#     datrank.rand.col <- t(apply(datrank.rand.col, 1, function(x) {rank(x, "average", na.last="keep")}))
#     ranksums <- colSums(datrank.rand.col, na.rm=TRUE) / nvals.cols
#     tail(rank(ranksums, "average", na.last=TRUE), 1)
#   }) - 1
#   ret <- hist(rand.dist, breaks=0:(ncol(datrank)+1), plot=FALSE)$counts
#   return(ret)
# }

random.ranks <- function(bs, n=1000, seed=123456, mc.cores=2)
{
  datrank.up <- datrank.dn <- NULL
  if(any(bs$up.dn > 0)) datrank.up <- t(apply(bs$dat[bs$up.dn > 0, , drop=FALSE], 1, function(x) {rank(x, "average", na.last="keep")}))
  if(any(bs$up.dn < 0)) datrank.dn <- ncol(bs$dat) - t(apply(bs$dat[bs$up.dn < 0, , drop=FALSE], 1, function(x) {rank(x, "average", na.last="keep")})) + 1
  datrank <- rbind(datrank.up, datrank.dn)
  nvals.cols <- nrow(datrank) - c(colSums(is.na(datrank)), 0)
  nvals.rows <- ncol(datrank) - rowSums(is.na(datrank))
  set.seed(seed)
  random.cols <- t(sapply(1:nrow(datrank), function(i) {runif(n, 0, nvals.rows[i] + 1)}))
  rand.dist <- unlist(mclapply(1:n, function(i) {
    datrank.rand.col <- cbind(datrank, random.cols[, i])
    datrank.rand.col <- t(apply(datrank.rand.col, 1, function(x) {rank(x, "average", na.last="keep")}))
    ranksums <- colSums(datrank.rand.col, na.rm=TRUE) / nvals.cols
    tail(rank(ranksums, "average", na.last=TRUE), 1)
  }, mc.cores=mc.cores)) - 1
  ret <- hist(rand.dist, breaks=0:(ncol(datrank)+1), plot=FALSE)$counts
  return(ret)
}


define.roi.regions <- function(bresat, random.rank.counts, middle.range=0.95){
  random.ranks.cdf <- cumsum(random.rank.counts) / sum(random.rank.counts)
  left <- max(c(0, which(random.ranks.cdf < (1 - middle.range) / 2)))
  right <- min(which(random.ranks.cdf > 1 - ((1 - middle.range) / 2))) - 1
  ret <- rep(2, length(bresat$rank))
  ret[bresat$rank < left] <- 1
  ret[bresat$rank > right] <- 3
  ret
}



# bresat.plot <- function(bs, clinical, random.rank.counts, middle.range=0.95,
#                                    title, row.labels, col.labels, clinical.labels,
#                                    add.ranksums=TRUE, add.rank.colors=TRUE, add.random.ranks=TRUE,
#                                    heights=c(title=0.15, heatmap=1, col.labels=0.07, clinical=0.3, ranksums=0.07, rank.colors=0.05, random.ranks=0.07),
#                                    widths=c(row.labels=0.15, heatmap=1, space=0.02, ns=0.02))
# {
#   stopifnot(ncol(bs$dat) == nrow(clinical))
#   user.heights <- heights
#   user.widths <- widths
#   
#   clinical <- clinical[bs$pat.order, ]
#   cols <- ncol(bs$dat)
#   ## reorder for consistency
#   bs$dat <- bs$dat[, cols:1]
#   clinical<-clinical[cols:1,]
#   bs$rank <- rev(bs$rank)
#   bs$ranksum <- rev(bs$ranksum)
#   bs$pat.order <- order(-bs$rank)
# 
#   
#   layout <- rbind(c("key", "title", "", ""),
#                   c("row.labels.rjust", "heatmap", "", "ns"),
#                   c("", "col.labels.rjust", "", ""))
#   widths <- c(row.labels=0.2, heatmap=1, space=0.02, ns=0.02)
#   heights=c(title=0.15, heatmap=1, col.labels=0.07)
#   
#   if (!missing(clinical))
#   {
#     layout <- rbind(layout, c("clinical.labels.rjust", "clinical", "", ""))
#     heights <- c(heights, clinical=0.3)
#   }
#   if (isTRUE(add.ranksums))
#   {
#     layout <- rbind(layout, c("ranksums.label", "ranksums", "", ""))
#     heights <- c(heights, ranksums=0.07)
#     
#   }
#   if (isTRUE(add.rank.colors))
#   {
#     layout <- rbind(layout, c("rank.colors.label", "rank.colors", "", ""))
#     heights <- c(heights, rank.colors=0.05)
#   }
#   if (isTRUE(add.random.ranks))
#   {
#     layout <- rbind(layout, c("random.ranks.label", "random.ranks", "", ""))
#     heights <- c(heights, random.ranks=0.07)
#   }
#   
#   user.heights <- user.heights[names(user.heights) %in% names(heights)]
#   user.widths <- user.widths[names(user.widths) %in% names(widths)]
#   heights[names(user.heights)] <- user.heights
#   widths[names(user.widths)] <- user.widths
#   
#   heatmap.simple(bs$dat, clinical=clinical,
#                  layout.mat=layout, widths=widths, heights=heights,
#                  row.clust=FALSE, col.clust=FALSE,
#                  title=title)
#   
#   if (isTRUE(add.ranksums))
#   {
#     ranksums<-viewport(y=0.45, height=0.9, xscale=c(0.5, cols + 0.5), yscale=c(0, max(bs$ranksum)), name="ranksums")
#     ndown <- downViewport("ranksums")
#     grid.polygon(c(1, 1:cols, cols), c(0, bs$ranksum[bs$pat.order], 0), gp=gpar(fill="red"), default.units="native")
#     popViewport()
#     upViewport(ndown)
#   }
#   if (isTRUE(add.rank.colors))
#   {
#     ndown <- downViewport("rank.colors")
#     heatmap.map(t(scale(1:cols)))
#     upViewport(ndown)
#     
#     ndown <- downViewport("rank.colors.label")
#     heatmap.labels("Rank", just="right")
#     upViewport(ndown)
#   }
#   if (isTRUE(add.random.ranks))
#   {
#     stopifnot(!missing(random.rank.counts) && length(random.rank.counts) == cols + 1
#               && middle.range < 1 && middle.range > 0)
#     tot <- sum(random.rank.counts)
#     random.ranks.cdf <- cumsum(random.rank.counts) / tot
#     left <- max(c(0, which(random.ranks.cdf < (1 - middle.range) / 2)))
#     right <- min(which(random.ranks.cdf > 1 - ((1 - middle.range) / 2))) - 1
#     
#     ## stopifnot(abs(cor(col.order, 1:cols) == 1))
#     ndown <- downViewport("heatmap.top.vp")
#     pushViewport(viewport(xscale=c(0, cols), layout.pos.row=2:nrow(layout), layout.pos.col=2))
#     grid.lines(left, c(0, 1), default.units="native", gp=gpar(col="black", lwd=2, lty=2))
#     grid.lines(right, c(0, 1), default.units="native", gp=gpar(col="black", lwd=2, lty=2))
#     popViewport()
#     upViewport(ndown)
#     
#     ndown <- downViewport("random.ranks")
#     pushViewport(viewport(xscale=c(0, cols)))
#     grid.lines(c(rep(0:cols, each=2), cols), c(0, rep(random.ranks.cdf, each=2)), default.units="native")
#     grid.rect(gp=gpar(lwd=0.5))
#     popViewport()
#     upViewport(ndown)
#     
#     ndown <- downViewport("random.ranks.label")
#     heatmap.labels("CDF random ranks", just="right")
#     upViewport(ndown)
#   }
#   
#   ndown <- downViewport("ns")
#   ## heatmap.map(matrix(c(bs$up.dn*3)))
#   pushViewport(viewport(yscale=c(length(bs$up.dn), 0)))
#   col.scheme <- heatmap.color.scheme()$col
#   grid.rect(0.5, diff(range(which(bs$up.dn > 0))) / 2 + 0.5, height=diff(range(which(bs$up.dn > 0))) + 1, gp=gpar(col=NA, fill=tail(col.scheme, 1)), default.units="native")
#   grid.rect(0.5, max(which(bs$up.dn < 0)) - (diff(range(which(bs$up.dn < 0))) / 2 + 0.5), height=diff(range(which(bs$up.dn < 0))) + 1, gp=gpar(col=NA, fill=head(col.scheme, 1)), default.units="native")
#   popViewport()
#   upViewport(ndown)
#   invisible(NULL)
# }
