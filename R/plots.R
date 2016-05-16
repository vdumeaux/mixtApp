### Contains functions for frequently used plot types

#source("../../src/heatmap.R")

## plot.patient.heatmap()
##
## plots heatmaps showing correct/incorrect predictions. The rows of the
## heatmap correspond to predictors while the columns correspond to
## patients. The columns will be ordered such that all good outcome
## patients will be to the left of all bad outcome patients. Within the
## good outcomes, the patients will be ordered by increasing median
## posteriors, while within the bad outcomes, the patients will be
## ordered by decreasing median posteriors. The order of the predictors
## will be preserved (hence the user should reorder 'posteriors'
## appropriately before calling this function).
##
## this function uses only functions from the grid package, so the user
## might have to call plot.new() prior to calling this function.
##
## After a call to the function, viewports are left in place with
## meaningful names so that the caller can easily access the different
## plotting areas and add/touch up the figure. The names of the
## viewports can be deduced by looking at the code.
##
## Argument:
##
##      posteriors
##
##          either a data.frame with columns representing predictors or
##          a matrix with rows representing predictors. The values
##          should be scores between 0 and 1 for predicting bad outcome
##          (higher scores indicating higher probability of bad
##          outcome). The ordering of the predictors will be preserved
##          while the patients will be reordered using the median of the
##          posteriors.
##
##      event
##
##          logical vector with TRUE indicating bad outcome
##
##      threshold
##
##          threshold to use for predicting bad outcomes
##
##      clinical
##
##          a clinical data frame matching the size of 'posteriors'
##
##      clinical.labels
##
##          labels describing the clinical variables. Defaults to
##          names(clinical).
##
##      predictor.labels
##
##          labels to be printed for the rows of the heatmap. Defaults
##          to the names given to the predictors in 'posteriors'.
##
##      predictor.stat
##
##          a numeric vector of same length as the number of
##          predictors. Will be used to sort the rows/predictors of the
##          patient heatmap. The predictors will be sorted in decreasing
##          order, the assumption being that higher values indicates
##          better predictor performance. If not given, the predictors
##          will not be reordered.
##
##      patient.stat
##
##          a numeric vector of same length as the number of
##          samples/patients. Will be used to order the columns/patients
##          of the patient heatmap. The patients will be sorted in
##          increasing order (but separately within good- and
##          bad-outcome patients). If not given, the median posteriors
##          will be used by default to order the patients.
##
##      plot.predictor.stat
##
##          logical. If TRUE and 'predictor.stat' has been provided, a
##          plot of the predictor stat will be drawn to the right of the
##          main patient heatmap.
##
##      predictor.stat.label
##
##          a label to be printed next to the plot of 'predictor.stat'
##
##      predictor.stat.range
##
##          the range of the y-axis of the plot of 'predictor.stat'
##
##      predictor.order.increasing
##
##          if TRUE, then order the predictors in increasing order of predictor.stat
##          otherwise, in decreasing order.
##
##      plot.patient.stat
##
##          logical. If TRUE, then a plot of 'sample.stat' will be drawn
##          below the main patient heatmap (and also below the heatmap
##          of the additional posteriors, if included).
##
##
##      patient.stat.label
##
##          a label to be printed next to the plot of 'patient.stat'
##
##      patient.stat.range
##
##          the range of y-axis of the plot of 'patient.stat'
##
##      patient.order.increasing
##
##          if TRUE, then order the patients in increasing order of
##          'patient.stat'.  otherwise, in decreasing order.
##
##      plot.patient.labels
##
##          logical. If TRUE, then patient labels will be plotted at the
##          bottom of the main heatmap, and below the patient stats if
##          any.
##
##      patient.labels
##
##          labels for the patients to be plotted if
##          'plot.patient.labels' is TRUE.
##
##
##      patient.labels.height
##
##          the height used in the heatmap layout for the patient
##          labels, given relative to the height of the main heatmap.
##
##      title
##
##          a title to be printed at the top of the plot
##
##      label.width
##
##          the width to use for the area where labels are printed (the
##          left side) given relative to the size of the main heatmap
##          width.
##
##      title.height
##
##          the hight used for the title area given relative to the size
##          of the heatmap height.
##
##      clinical.height.multiplier
##
##          The height of each row of the clinical variables is set to
##          the same height as for the individual signatures multiplied
##          by 'clinical.height.multiplier'
##
##      predictor.stat.width
##
##          the width to use for the area where 'predictor.stat' is
##          plotted given relative the size of the main heatmap width.
##
##      additional.posteriors
##
##          optional. same kind of data as 'posteriors', except that these will
##          not affect the ordering of the patients and will be plotted
##          as a separate heatmap below the heatmap of 'posteriors'.
##
##      additional.predictor.labels
##
##          labels to be printed for the rows of the heatmap
##          corresponding to 'additional.posteriors'.
##
##      just.labels
##
##          do not plot the heatmap or clinical variables
##          (necessary to avoid rasterized text in downstream photoshop editing)
##
##      ...
##
##          additional arguments passed to heatmap.simple() when
##          plotting the additional posteriors. Used mainly to set
##          'na.color'.
##
##  Returns:
##
##      list returned by heatmap.simple
plot.patient.heatmap <- function(posteriors, event, threshold=0.5,
                                 clinical=NULL,
                                 clinical.labels=NULL,
                                 predictor.labels=NULL,
                                 predictor.stat=NULL,
                                 patient.stat=NULL,
                                 plot.predictor.stat=FALSE,
                                 predictor.stat.label="",
                                 predictor.stat.range=c(0, 1),
                                 predictor.order.increasing=TRUE,
                                 plot.patient.stat=FALSE,
                                 patient.stat.label="Median posteriors",
                                 patient.stat.range=c(0, 1),
                                 patient.order.increasing=TRUE,
                                 plot.patient.labels=FALSE,
                                 patient.labels=NULL,
                                 patient.labels.height=0.1,
                                 title=NULL,
                                 label.width=0.3,
                                 title.height=0.1,
                                 clinical.height.multiplier=1,
                                 predictor.stat.width=0.1,
                                 additional.posteriors=NULL,
                                 additional.predictor.labels=NULL,
                                 just.labels=FALSE,
                                 ...)
{
    library(matrixStats)

    ## check input and assign default values where needed
    stopifnot(class(posteriors) == "data.frame" || class(posteriors) == "matrix")
    if (class(posteriors) == "data.frame")
    {
        posteriors <- t(as.matrix(posteriors))
    }
    stopifnot(is.numeric(posteriors))


    stopifnot(is.vector(event, "logical") && length(event) == ncol(posteriors))

    if (!is.null(clinical))
    {
        stopifnot(class(clinical) == "data.frame" && nrow(clinical) == ncol(posteriors))
        if (is.null(clinical.labels))
            clinical.labels <- names(clinical)
        else
            stopifnot(length(clinical.labels) == length(clinical))
    }

    if (is.null(predictor.labels))
        predictor.labels <- rownames(posteriors)
    else
        stopifnot(length(predictor.labels) == nrow(posteriors))

    stopifnot(is.null(predictor.stat) || (is.numeric(predictor.stat) && length(predictor.stat) == nrow(posteriors)))

    stopifnot(is.null(patient.stat) || (is.numeric(patient.stat) && length(patient.stat) == ncol(posteriors)))
    if (is.null(patient.stat))
        patient.stat <- colMedians(x=posteriors, na.rm=TRUE)

    if (!is.null(additional.posteriors))
    {
        stopifnot(class(additional.posteriors) %in% c("data.frame", "matrix"))
        if (class(additional.posteriors) == "data.frame")
            additional.posteriors <- t(as.matrix(additional.posteriors))
        stopifnot(is.numeric(additional.posteriors))
        stopifnot(ncol(additional.posteriors) == ncol(posteriors))
        if (is.null(additional.predictor.labels))
            additional.predictor.labels <- rownames(additional.posteriors)
        else
            stopifnot(length(additional.predictor.labels) == nrow(additional.posteriors))
    }

    ## order rows (i.e., predictors) of 'posteriors'
    if (!is.null(predictor.stat))
    {
        row.order <- order(predictor.stat, decreasing=!predictor.order.increasing)
        posteriors <- posteriors[row.order, , drop=FALSE]
        predictor.stat <- predictor.stat[row.order]
        predictor.labels <- predictor.labels[row.order]
    }

    ## order columns (i.e., patients) of 'posteriors'
    col.order <- order(patient.stat, decreasing=!patient.order.increasing)
    posteriors <- posteriors[, col.order, drop=FALSE]
    event <- event[col.order]
    patient.stat <- patient.stat[col.order]
    if (!is.null(clinical))
        clinical <- clinical[col.order, , drop=FALSE]
    if (!is.null(additional.posteriors))
        additional.posteriors <- additional.posteriors[, col.order, drop=FALSE]
    if (!is.null(patient.labels))
        patient.labels <- patient.labels[col.order]

    ## and separate the good and bad outcomes
    col.order <- order(event)
    posteriors <- posteriors[, col.order, drop=FALSE]
    event <- event[col.order]
    patient.stat <- patient.stat[col.order]
    if (!is.null(clinical))
        clinical <- clinical[col.order, , drop=FALSE]
    if (!is.null(additional.posteriors))
        additional.posteriors <- additional.posteriors[, col.order, drop=FALSE]
    if (!is.null(patient.labels))
        patient.labels <- patient.labels[col.order]

    ## compute values of posteriors for plotting heatmaps. values should be:
    ## correctly predicted good: -2,
    ## incorrectly predicted good: -1,
    ## correctly predicted bad: 2,
    ## incorrectly predicted bad: 1
    if (ncol(posteriors) < 2)
        posteriors.rounded <- as.matrix(apply(posteriors, 1, function(p) {event * 3 - 2 + (p > threshold)}))
    else
        posteriors.rounded <- t(apply(posteriors, 1, function(p) {event * 3 - 2 + (p > threshold)}))

    ## create the heatmap layout, widths, and heights depending on what
    ## is supposed to be included in the plot
    layout <- matrix(c("", "", "",
                       "row.labels.rjust", "heatmap", "",
                       "", "", ""), nrow=3, byrow=TRUE)
    widths <- c(label.width, 1, 0.02)
    heights <- c(0.02, 1, 0.02)

    if (!is.null(additional.posteriors))
    {
        layout <- rbind(layout, c("additional.labels", "additional", ""))
        layout <- rbind(layout, c("", "", ""))
        heights <- c(heights[-length(heights)], 1/nrow(posteriors), nrow(additional.posteriors)/nrow(posteriors), tail(heights, 1))
    }

    if (!is.null(patient.stat) && isTRUE(plot.patient.stat))
    {
        layout <- rbind(layout, c("patient.stat.label", "patient.stat", ""))
        layout <- rbind(layout, c("", "", ""))
        heights <- c(heights[-length(heights)], 1/nrow(posteriors), 5/nrow(posteriors), tail(heights, 1))
    }

    if (!is.null(patient.labels) && isTRUE(plot.patient.labels))
    {
        layout <- rbind(layout, c("", "col.labels.ljust", ""))
        layout <- rbind(layout, c("", "", ""))
        heights <- c(heights[-length(heights)], 1/nrow(posteriors), patient.labels.height, tail(heights, 1))
    }

    if (!is.null(clinical))
    {
        layout <- rbind(layout, c("clinical.labels.rjust", "clinical", ""))
        layout <- rbind(layout, c("", "", ""))
        heights <- c(heights[-length(heights)], 1/nrow(posteriors), (length(clinical) / nrow(posteriors))*clinical.height.multiplier, tail(heights, 1))
    }

    if (!is.null(title))
    {
        layout[1, 2] <- "title"
        heights[1] <- title.height
    }

    if (!is.null(predictor.stat) && isTRUE(plot.predictor.stat))
    {
        idx <- which(layout == "heatmap", arr.ind=TRUE)[1]
        layout <- cbind(layout, matrix("", nrow=nrow(layout), ncol=2))
        layout[idx, ncol(layout) - 1] <- "predictor.stat"
        layout[idx, ncol(layout)] <- "predictor.stat.label"
        widths <- c(widths, predictor.stat.width, predictor.stat.width/3)
    }

    if (just.labels) {
      layout[which(layout == "heatmap")] <- ""
      layout[which(layout == "clinical")] <- ""
    }

    ## draw the main patient heatmap possibly together with the clinical
    ## variables and title.
    ret <- heatmap.simple(posteriors.rounded,
                          layout.mat=layout,
                          widths=widths,
                          heights=heights,
                          scale="none",
                          row.labels=predictor.labels,
                          row.clust=FALSE,
                          col.clust=FALSE,
                          clinical=clinical,
                          clinical.labels=clinical.labels,
                          col.labels=patient.labels,
                          title=title)

    ## plot 'patient.stat' below the main heatmap
    if (!is.null(patient.stat) && isTRUE(plot.patient.stat))
    {
        idx <- which(layout == "patient.stat", arr.ind=TRUE)
        downViewport("heatmap.top.vp")
        pushViewport(viewport(layout.pos.row=idx[1], layout.pos.col=idx[2],
                              xscale=c(0.5, ncol(posteriors) + 0.5),
                              yscale=patient.stat.range))
        grid.rect()
        grid.polyline(rep(c(0, 1), 4), rep(c(0.2, 0.4, 0.6, 0.8), each=2), id.lengths=rep(2, 4), gp=gpar(lwd=0.2))
        grid.lines(1:ncol(posteriors), patient.stat, gp=gpar(col="blue"), default.units="native")
        popViewport()

        idx <- which(layout == "patient.stat.label", arr.ind=TRUE)
        pushViewport(viewport(layout.pos.row=idx[1], layout.pos.col=idx[2]))
        heatmap.labels(labels=patient.stat.label, type="row.labels", just="right")
        popViewport()

        upViewport()
    }

    ## plot 'predictor.stat' to the right of the main heatmap
    if (!is.null(predictor.stat) && isTRUE(plot.predictor.stat))
    {
        idx <- which(layout == "predictor.stat", arr.ind=TRUE)
        downViewport("heatmap.top.vp")
        pushViewport(viewport(layout.pos.row=idx[1], layout.pos.col=idx[2],
                              xscale=predictor.stat.range,
                              yscale=c(nrow(posteriors) + 0.5, 0.5)))
        grid.rect()
        grid.polyline(rep(c(0.2, 0.4, 0.6, 0.8), each=2), rep(c(0, 1), 4), id.lengths=rep(2, 4), gp=gpar(lwd=0.2))
        grid.polygon(c(0, predictor.stat), c(1, 1:nrow(posteriors)),
                     gp=gpar(fill="blue", col="blue", alpha=0.5), default.units="native")
        popViewport()

        idx <- which(layout == "predictor.stat.label", arr.ind=TRUE)
        pushViewport(viewport(layout.pos.row=idx[1], layout.pos.col=idx[2]))
        heatmap.labels(labels=predictor.stat.label, type="col.labels", just="centre", rot=180)
        popViewport()

        upViewport()
    }

    ## plot heatmap of any additional posteriors below the main heatmap
    if (!is.null(additional.posteriors))
    {
        labels.idx <- which(layout == "additional.labels")
        heatmap.idx <- which(layout == "additional")
        layout.additional <- layout
        layout.additional[1:length(layout)] <- ""
        layout.additional[labels.idx] <- "row.labels.rjust"
        layout.additional[heatmap.idx] <- "heatmap"

        if (ncol(additional.posteriors) < 2)
            additional.posteriors.rounded <- as.matrix(apply(additional.posteriors, 1, function(p) {event * 3 - 2 + (p > threshold)}))
        else
            additional.posteriors.rounded <- t(apply(additional.posteriors, 1, function(p) {event * 3 - 2 + (p > threshold)}))

        if (just.labels)
          layout.additional[which(layout.additional == "heatmap")] <- ""

        heatmap.simple(additional.posteriors.rounded,
                       layout.mat=layout.additional,
                       widths=widths,
                       heights=heights,
                       scale="none",
                       row.labels=additional.predictor.labels,
                       row.clust=FALSE,
                       col.clust=FALSE,
                       title="", ...)
    }

    ## leave behind a viewport tree that the caller can use to access
    ## the different plotting areas.
    pushViewport(viewport(layout=grid.layout(nrow(layout), ncol(layout), widths=widths, heights=heights),
                          name="heatmap.top.vp"))
    apply(which(layout != "", arr.ind=TRUE), 1, function(idx) {
        pushViewport(viewport(layout.pos.row=idx[1], layout.pos.col=idx[2], name=layout[idx[1], idx[2]]))
        upViewport()
    })
    upViewport()

    return(invisible(ret))
}

## plot.rank.heatmap()
##
## plots rank matrices, usually with each row corresponding to ranks of
## patients obtained from a gene set using bresat. 'col.ref' and
## 'flip.ranks' will determine the overall results.
##
## Arguments:
##
##      rank.matrix
##
##          a matrix where each row consists of ranks. It is assumed,
##          and is very important to obtain correct plots, that the
##          ranks are in the range 1..ncol(rank.matrix).
##
##      col.ref, flip.ranks
##
##          These two parameters together determine how the rows and
##          columns are ordered and how the ranks are displayed.
##
##          If given, 'col.ref' must be either a numeric or logical
##          vector with length equal to the number of columns of
##          'rank.matrix'. It is used to determine the column ordering.
##
##          If not specified, column medians will be used to order the
##          columns. First, the ranks of rows are reversed if their
##          correlation to the median ranks are negative. After this
##          step, the final median rank is obtained and used to order
##          the columns. Rows will be ordered based on correlation to
##          the final median ranks. If 'flip.ranks' is TRUE, then the
##          heatmap will display the ranks after the reversal of some
##          ranks, otherwise, the original ranks are displayed.
##
##          If 'col.ref' is a numeric vector, then it will be used to
##          order the columns (in increasing order). The rows of the
##          rank matrix will be ordered based on correlation with this
##          vector. If 'flip.ranks' is TRUE, then ranks of rows with
##          negative correlation to 'col.ref' will be reversed.
##
##          If 'col.ref' is a logical vector, it is assumed to be the
##          associated event vector. The columns will be grouped by
##          event and an attempt will be made to create a nice column
##          and row ordering based on median ranks and correlation to
##          the "best" row, respectively. Ranks may be reversed if that
##          induces a better association with 'col.ref'. If 'flip.ranks'
##          is TRUE, the heatmap will display the ranks after reversing
##          of the ranks. Otherwise, the original ranks are displayed.
##
##      row.labels
##
##          lables to be display for the rows of 'rank.matrix'
##
##      clinical
##
##          data.frame with colors to be displayed under the heatmap
##
##      clinical.labels
##
##          names of the clinical variables. If not specified, the names
##          of 'clinical' are used.
##
##      title
##
##          optional title to be plotted at the top
##
##      label.width, title.height
##
##          space to be reserved for displaying row and clinical labels
##          to the left of the heatmap and the for the title at the
##          top. The value is interpreted as units of one heatmap width/height.
##
##  Returns:
##
##      invisible(NULL)
plot.rank.heatmap <- function(rank.matrix, col.ref=NULL, flip.ranks=TRUE, row.labels=NULL,
                              clinical=NULL, clinical.labels=NULL, title=NULL,
                              label.width=0.3, title.height=0.05)
{
    if (class(rank.matrix) != "matrix")
        stop("'rank.matrix' must be a a matrix")
    if (!is.numeric(rank.matrix))
        stop("'rank.matrix' must be a numeric matrix")

    if (!is.null(col.ref))
    {
        if (length(col.ref) != ncol(rank.matrix))
            stop("size of 'col.ref' does not match 'rank.matrix'")
        if (!is.numeric(col.ref) && !is.logical(col.ref))
            stop("bad 'col.ref'")
    }

    if (!is.null(clinical))
    {
        if (class(clinical) != "data.frame")
            stop("'clinical' must be a data.frame")
        if (nrow(clinical) != ncol(rank.matrix))
            stop("size of 'clinical' does not match size of 'rank.matrix'")
    }

    if (is.null(row.labels))
    {
        row.labels <- rownames(rank.matrix)
    }
    else
    {
        if (length(row.labels) != nrow(rank.matrix))
            stop("size of 'row.labels' does not match size of 'rank.matrix'")
    }

    if (!is.null(clinical))
    {
        if (is.null(clinical.labels))
            clinical.labels <- names(clinical)
        else if (length(clinical.labels) != length(clinical))
            stop("size of 'clinical.labels' does not match size of 'clinical'")
    }

    if (!is.null(title) && !is.character(title))
        stop("'title' must be of type character")

    ## create a layout for the heatmaps based on whether or not
    ## 'clinical' and 'title' have been given.
    layout <- matrix(c("", "", "",
                       "row.labels.rjust", "heatmap", "",
                       "", "", ""), nrow=3, byrow=TRUE)
    widths <- c(label.width, 1, 0.02)
    heights <- c(0.02, 1, 0.02)

    if (!is.null(clinical))
    {
        layout <- rbind(layout, c("clinical.labels.rjust", "clinical", ""))
        layout <- rbind(layout, c("", "", ""))
        heights <- c(heights[-length(heights)], 1/nrow(rank.matrix), length(clinical) / nrow(rank.matrix), tail(heights, 1))
    }

    if (!is.null(title))
    {
        layout[1, 2] <- "title"
        heights[1] <- title.height
    }

    if (is.null(col.ref))
    {
        ## use median ranks for columns to get the majority direction.
        meta.rank <- rank(apply(rank.matrix, 2, median))
        if (length(unique(meta.rank)) == 1)
            meta.rank <- rank.matrix[1, ]

        ## flip rows that are opposite to the majority
        flipped.matrix <- t(apply(rank.matrix, 1, function(x) {
            if (var(meta.rank) == 0 || var(x) == 0 || cor(x, meta.rank) > 0)
                x
            else
                length(x) - x + 1
        }))

        ## get a final rank for the columns based on the flipped matrix
        meta.rank <- rank(apply(flipped.matrix, 2, median))
        if (length(unique(meta.rank)) == 1)
            meta.rank <- rank.matrix[1, ]

        ## should the flipped matrix be displayed?
        if (isTRUE(flip.ranks))
            rank.matrix <- flipped.matrix

        ## order the rows based on correlation to the final column ranks
        cors <- apply(rank.matrix, 1, cor, meta.rank)
        row.order <- order(cors, decreasing=TRUE)
        rank.matrix <- rank.matrix[row.order, ]
        row.labels <- row.labels[row.order]
        cors <- cors[row.order]

        ## order the columns based on the final column ranks
        col.order <- order(meta.rank)
        rank.matrix <- rank.matrix[, col.order, drop=FALSE]
        if (!is.null(clinical))
            clinical <- clinical[col.order, ]
    }
    else if (is.logical(col.ref))
    {
        ## use the median of the ranks of the bad outcomes to decide which
        ## rows to flip.
        flipped.matrix <- t(apply(rank.matrix, 1, function(x) {
            flipped <- length(x) - x + 1
            s1 <- median(x[col.ref], na.rm=TRUE)
            s2 <- median(flipped[col.ref], na.rm=TRUE)
            if (s2 > s1)
                flipped
            else
                x
        }))

        ## should we display the flipped matrix?
        if (isTRUE(flip.ranks))
            rank.matrix <- flipped.matrix

        ## order rows twice, first by the medians of the bad outcomes,
        ## then by their correlation to the first (best) row.
        row.order <- order(apply(rank.matrix, 1, function(x) {median(x[col.ref], na.rm=TRUE)}), decreasing=TRUE)
        rank.matrix <- rank.matrix[row.order, , drop=FALSE]
        flipped.matrix <- flipped.matrix[row.order, , drop=FALSE]
        row.labels <- row.labels[row.order]

        cors <- apply(rank.matrix, 1, cor, rank.matrix[1, ])
        row.order <- c(which(cors > 0), which(cors <= 0))
        rank.matrix <- rank.matrix[row.order, , drop=FALSE]
        flipped.matrix <- flipped.matrix[row.order, , drop=FALSE]
        row.labels <- row.labels[row.order]

        ## order the columns twice, first by the median of the flipped
        ## matrix, then by the event vector
        col.order <- order(apply(flipped.matrix, 2, median))
        rank.matrix <- rank.matrix[, col.order, drop=FALSE]
        flipped.matrix <- flipped.matrix[, col.order, drop=FALSE]
        col.ref <- col.ref[col.order]
        if (!is.null(clinical))
            clinical <- clinical[col.order, ]

        col.order <- c(which(!col.ref), which(col.ref), which(is.na(col.ref)))
        rank.matrix <- rank.matrix[, col.order, drop=FALSE]
        flipped.matrix <- flipped.matrix[, col.order, drop=FALSE]
        if (!is.null(clinical))
            clinical <- clinical[col.order, ]
    }
    else ## col.ref is numeric
    {
        ## use the provided 'col.ref' as reference for ordering columns
        meta.rank <- rank(col.ref)

        ## flip rows that have negative correlation to 'col.ref'
        flipped.matrix <- t(apply(rank.matrix, 1, function(x) {
            if (var(x) == 0 || cor(x, meta.rank) > 0)
                x
            else
                length(x) - x + 1
        }))

        ## should the flipped matrix be displayed?
        if (isTRUE(flip.ranks))
            rank.matrix <- flipped.matrix

        ## order rows based on correlation to provided reference.
        cors <- apply(rank.matrix, 1, cor, meta.rank)
        row.order <- order(cors, decreasing=TRUE)
        rank.matrix <- rank.matrix[row.order, ]
        row.labels <- row.labels[row.order]
        cors <- cors[row.order]

        ## order columns based on provided reference.
        col.order <- order(meta.rank)
        rank.matrix <- rank.matrix[, col.order, drop=FALSE]
        if (!is.null(clinical))
            clinical <- clinical[col.order, ]
    }

    ## finally, plot the damn thing
    ret <- heatmap.simple(rank.matrix,
                          layout.mat=layout,
                          widths=widths,
                          heights=heights,
                          scale=function(x) {scale(x) * 1.2},
                          row.labels=row.labels,
                          row.clust=FALSE,
                          col.clust=FALSE,
                          clinical=clinical,
                          clinical.labels=clinical.labels,
                          title=title)
    return(invisible(ret))
}

## plot.coxph.by.cohort()
##
## plots a heatmap of the coxph p-values with rows representing genes
## and columns representing data sets. The p-values are -log10
## transformed.
##
## Arguments:
##
##      datasets
##
##          list of data sets with standard layouts
##
##      cohort
##
##          character with the name of the cohort
##
##      gene.set
##
##          character vector with names of genes
##
##      title
##
##          optional title for the plot
##
##      breaks
##
##          the breaks used for the color scheme for the
##          -log10-transformed p-values.
##
##      row.label.width, col.label.width, title.height
##
##          width and hieght parameters of the plot
##
## Returns:
##
##      invisibly returns the -log10-transformed p-value matrix
plot.coxph.by.cohort <- function(datasets, cohort, gene.set, title = NULL,
                                 breaks = -log10(c(1, 0.5, 0.1, 0.05, 0.01, 0.001, 0.0001)),
                                 row.label.width=0.1, col.label.width=0.1, title.height=0.1)
{
    pval.matrix <- sapply(datasets, function(dat) {
        coxph <- huc.get.coxph(dat, cohort, gene.names=gene.set, rm.missing.genes=FALSE)
        pval <- coxph$pval
        pval[pval == 0] <- .Machine$double.eps
        ret <- -log10(pval)
        dn.idx <- which(coxph$coef < 0)
        if (length(dn.idx) > 0)
            ret[dn.idx] <- ret[dn.idx] * -1
        return(ret)
    })

    color.scheme <- heatmap.color.scheme(low.breaks = rev(-breaks), high.breaks = breaks)
    layout.mat <- rbind(c("key", "col.labels.ljust"),
                        c("", ""),
                        c("row.labels.rjust", "heatmap"))
    widths <- c(row.label.width, 1)
    heights <- c(col.label.width, 1/nrow(pval.matrix), 1)
    if (!is.null(title))
    {
        heights <- c(title.height, heights)
        layout.mat <- rbind(c("", "title"), layout.mat)
    }

    ret <- heatmap.simple(pval.matrix,
                          color.scheme = color.scheme,
                          layout.mat=layout.mat,
                          widths=widths,
                          heights=heights,
                          row.labels = gene.set,
                          col.labels = bresect.name.map[names(datasets)],
                          row.clust=FALSE,
                          col.clust=FALSE,
                          scale="none",
                          title=title)
    ret$pval.matrix <- pval.matrix
    return(invisible(ret))
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


### plot palette

pal <- function(col, border = "light gray", ...)
{
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}
