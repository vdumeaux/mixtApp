### this file contains functions used to load, clean, and otherwise work
### on the data sets that we use in this project.

huc.variable.categorical <- function(variable) {
  variable %in% c("lymph","grade","er","her2","stage","event","event.5",
                  "pam50.genefu","hybrid.genefu","aims","intclust","claudin.low","chemo","tamoxifen","herceptin",
                  "lum", "lumN", "prolif", "basalL",
                  "MKS", "ERS", "LUMS", "HER2S",
                  "consistently.classified","type","dataset", 
                  "Normal", "cit", "orig.dataset", "TissueType")
}

huc.clinical.colors <- function() {
  library(RColorBrewer)

  heatmap.clinical <- NULL
  heatmap.clinical[["consistently.classified"]] <- c(ccgc="blue", ccgi="pink", ccbc="red", ccbi="lightblue", ncc="white")
  heatmap.clinical[["event.5"]] <- c(Good="white",Bad="red")
  heatmap.clinical[["event"]] <- c(Good="white",Bad="red")
  heatmap.clinical[["er"]] <- c(Negative="white",Positive="green")
  heatmap.clinical[["her2"]] <- c(Negative="white", Positive="orange")
  heatmap.clinical[["pam50.genefu"]] <- c(LumA="blue4", LumB="deepskyblue", Basal="firebrick2", Normal="green4", Her2="hotpink2")
  heatmap.clinical[["pam50.parker"]] <- c(LumA="blue4", LumB="deepskyblue", Basal="firebrick2", Normal="green4", Her2="hotpink2")
  heatmap.clinical[["aims"]] <- c(LumA="blue4", LumB="deepskyblue", Basal="firebrick2", Normal="green4", Her2="hotpink2")
  heatmap.clinical[["hybrid.genefu"]] <- c(luma.erp="blue4", lumb.erp="deepskyblue", basal.erp="firebrick2", normal.erp="green4", her2pam50p.erp="hotpink2", her2clinicaln.ern="tan", her2clinicalp.ern="orange")
  #heatmap.clinical[["intclust"]] <- c("tomato2", "palegreen3", "violetred3", "turquoise3", "firebrick4", "yellow2", "royalblue4", "darkgoldenrod1", "plum3", "purple3")
  heatmap.clinical[["lum"]]<-c(No="white", Yes="blue")
  heatmap.clinical[["lumN"]]<-c( No="white", Yes="darkturquoise")
  heatmap.clinical[["prolif"]]<-c( No="white", Yes="magenta")
  heatmap.clinical[["basalL"]]<-c(No="white", Yes="firebrick2")
  #names(heatmap.clinical[["intclust"]]) <- paste("intclust", 1:10, sep="")
  heatmap.clinical[["claudin.low"]] <- c(No="white", Yes="yellow3")
  heatmap.clinical[["grade"]] <- c(Low="white",Intermediate="pink",High="red")
  heatmap.clinical[["lymph"]] <- c(Negative="white",Positive="red")
  heatmap.clinical[["type"]] <- c(IDC="red","IDC/ILC"="purple",ILC="blue",Other="green")
  heatmap.clinical[["orig.dataset"]] <- c(cc1.blood="red", cc2.blood="purple",cc3.blood="blue",cc4.blood="green", biopsy="orange", cc.blood="firebrick2")
  heatmap.clinical[["TissueType"]] <- c(tumor="yellow", blood="red", blood.normal="pink")
  #heatmap.clinical[["chemo"]] <- c(No="white",Yes="red")
  #heatmap.clinical[["tamoxifen"]] <- c(No="white",Yes="red")
  #heatmap.clinical[["herceptin"]] <- c(No="white",Yes="red")
  heatmap.clinical[["cit"]] <- c(basL="firebrick2",lumA="blue4",lumB="deepskyblue", normL="green4",lumC="hotpink2",mApo="tan")
  heatmap.clinical[["lehmann"]] <- c(UNS="#F11A14", BL1="#F698C2", BL2="#3651A5", IM="#F99900", M="#BABBBE", MSL="#9B3067", LAR="#66BF3C")
  heatmap.clinical[["menopause"]]<-c(No="white", Yes="red")
  heatmap.clinical[["hrt"]]<-c(No="white", Yes="red")
  heatmap.clinical[["weight.70.plus"]]<-c(No="white", Yes="red")
  heatmap.clinical[["age.55.plus"]]<-c(No="white", Yes="red")
  heatmap.clinical[["smoking"]]<-c(No="white", Yes="red")
  heatmap.clinical[["medication"]]<-c(No="white", Yes="red")
  heatmap.clinical[["Normal"]]<-c(No="red", Yes="pink")
  
  ## age? size? dcis? inherently difficult? claudin low/high?
  hospitals <- c("Tromso", "Tonsberg", "Ulleval","Bodo", "Fredrikstad",
                 "Haukeland", "Molde", "Radiumhospitalet","Tronsdheim", 
                 "Stavanger")
  hospital.colors <- hcl(rep(seq(0, 350, 30), each=2), c(80, 50), c(50, 65, 50, 40))[1:length(hospitals)]
  names(hospital.colors) <- hospitals
  heatmap.clinical[["hospital"]] <- hospital.colors
  
  datasets <- huc.ls()
  dataset.colors <- hcl(rep(seq(0, 350, 30), each=2), c(80, 50), c(50, 65, 50, 40))[1:length(datasets)]
  names(dataset.colors) <- datasets
  heatmap.clinical[["dataset"]] <- dataset.colors
  
  ## if you add something here, make sure you give it a 'proper' name in bresect.name.map in common.R
  v <- setdiff(names(heatmap.clinical), names(bresect.name.map))
  if (length(v) > 0) {
    warning(paste("clinical variable(s) not given proper names in bresect.name.map:", paste(v, collapse=", ")))
    map.names <- names(bresect.name.map)
    bresect.name.map <- c(bresect.name.map, v)
    names(bresect.name.map) <- c(map.names, v)
  }

  heatmap.clinical
}

## huc.color.clinical()
##
## Given a standard dataset clinical data frame, returns a corresponding
## data frame with colors corresponding to the values in the clinical
## data frame.
huc.color.clinical <- function(clinical) {
  colors <- huc.clinical.colors()

  clinical1 <- clinical[,which(colnames(clinical) %in% names(colors)), drop=FALSE]
  if (ncol(clinical1) == 0) return(NULL)

  colors <- colors[which(names(colors) %in% colnames(clinical1))] ## remove variables in colors not in clinical
  clinical1 <- clinical1[,names(colors), drop=FALSE] ## order variables by their order in colors

  for (variable in names(colors)) {
    if (is.logical(clinical1[,variable]))
      clinical1[,variable] <- colors[[variable]][sign(clinical1[,variable])+1]
    else if (is.factor(clinical[,variable]))
      clinical1[,variable] <- colors[[variable]][as.character(clinical1[,variable])]
    else
      clinical1[,variable] <- colors[[variable]][clinical1[,variable]]
  }

  mks.colors <-sequential_hcl(n=length(clinical$MKS[!is.na(clinical$MKS)]), h=c(0,30))
  names(mks.colors)<-rownames(clinical)[order(clinical$MKS, decreasing=T, na.last=NA)]
  clinical1$MKS <- mks.colors[rownames(clinical)]
  
    ers.colors <-sequential_hcl(n=length(clinical$ERS[!is.na(clinical$ERS)]))
    names(ers.colors)<-rownames(clinical)[order(clinical$ERS, decreasing=T, na.last=NA)]
    clinical1$ERS <- ers.colors[rownames(clinical)]
    
    lums.colors <-sequential_hcl(n=length(clinical$LUMS[!is.na(clinical$LUMS)]))
    names(lums.colors)<-rownames(clinical)[order(clinical$LUMS, decreasing=T, na.last=NA)]
    clinical1$LUMS <- lums.colors[rownames(clinical)]
    
    her2s.colors <-sequential_hcl(n=length(clinical$HER2S[!is.na(clinical$HER2S)]), h=c(70,90))
    names(her2s.colors)<-rownames(clinical)[order(clinical$HER2S, decreasing=T, na.last=NA)]
    clinical1$HER2S <- her2s.colors[rownames(clinical)]
  
  return(clinical1)
}


### huc.ls()
###
### Returns:
###
###     a character vector with the names of all datasets that
###     are available for loading using the huc.load() function.
huc.ls <- function()
{    return(c("anders",
             "bild",
	     "boersma.epi",
	     "boersma.stroma",
             "curtis",
             "curtis.discovery",
             "curtis.discovery.1",
             "curtis.discovery.2",
             "curtis.validation",
             "curtis.validation.1",
             "curtis.validation.2",
             "guedj",
	     "li",
             "mcgill.epi",
             "mcgill.gq",
             "mcgill.str",
             "muggerud",
             "nki",
             "nowac.luma",
              "cc1.blood",
	     "cc2.blood",
	     "cc3.blood",
	     "cc4.blood",
	     "biopsy",
       "sample_matched_biopsy",
             "parker",
	     "servant",
             "sotiriou",
             "tcga",
             "vanvliet",
             "wang"))
}

### huc.load()
###
### loads data sets and rearranges them to ensure a uniform layout
### w.r.t. type and name of variables. Allows multiple ways of cleaning
### up the data sets based on missin clinical values and censoring.
###
### Arguments:
###
###     dat.names:
###
###         character vector of data set names as obtained from huc.ls()
###
###     data.dir:
###
###         the top data directory of the project
###
###     rm.na.fields:
###
###         a character vector with names of clinical variables that
###         must not be missing. Patients missing these variables will
###         be removed. Note that patients are removed *after* the coxph
###         computation so that as many patients as possible are used
###         for determining prognostic genes.
###
###     rm.censored.5:
###
###         logical. If TRUE, then all patients whose events are unknown
###         at 5 years (time < 60 and event == FALSE) are removed. The
###         outcome of such patients are essentially unkown to us and
###         they should not be used during training of NBCs. Note,
###         however, that the coxph coefficients and p-values are
###         computed before these patients are removed in order to use
###         all the available data.
###
### Returns:
###
###     A list of data set objects. The names of the elements will be
###     the same as 'dat.names'
huc.load <- function(dat.names, data.dir, rm.na.fields=c("time", "event", "er", "her2", "pam50.genefu"),
                     rm.censored.5=TRUE)
{
    stopifnot(typeof(dat.names) == "character")
    stopifnot(all(dat.names %in% huc.ls()))

    ## load the data set objects, but don't load anything for the
    ## curtis.discovery* and curtis.validation* sets. We will handle those
    ## later by splitting curtis up.
    dats.to.load <- dat.names
    curtis.subsets <- regexpr("curtis", dat.names) > 0
  
    if (any(curtis.subsets))
    {
        dats.to.load <- dat.names[!curtis.subsets]
        if (!"curtis" %in% dats.to.load)
            dats.to.load <- c(dats.to.load, "curtis")
    }
    paths <- sapply(dats.to.load, function(s) {file.path(data.dir, "huc", paste(s, ".rda", sep=""))})
    dataset.objects <- lapply(paths, function(p) {get(load(p))})
    names(dataset.objects) <- dats.to.load

    ## rearrange the objects by passing each data set to its own special rearrange function
    rearrange.funs <- lapply(paste("huc.rearrange.", dats.to.load, sep=""), get)
    datasets <- mapply(function(dat, fun) {fun(dat, data.dir)}, dataset.objects, rearrange.funs, SIMPLIFY=FALSE)

    ## add authors' pam50 assignments if any
    for (ii in 1:length(datasets))
    {
        datasets[[ii]]$clinical$pam50.authors <- NA_character_
        if (datasets[[ii]]$name == "curtis")
        {
            idx <- match(datasets[[ii]]$clinical$id.orig, dataset.objects[[ii]]$clinical$METABRIC_ID)
            datasets[[ii]]$clinical$pam50.authors <- dataset.objects[[ii]]$clinical$Pam50Subtype[idx]
            datasets[[ii]]$clinical$pam50.authors[datasets[[ii]]$clinical$pam50.authors == "NC"] <- NA
        }
    }


    ## add the in-house pam50 subtypes
    for (ii in 1:length(datasets))
    {
        datasets[[ii]]$clinical$pam50.inhouse <- huc.pam50.inhouse(datasets[[ii]], data.dir)
    }

    ## compute the pam50 subtypes using the genefu package
    for (ii in 1:length(datasets))
    {
        datasets[[ii]]$clinical$pam50.genefu <- huc.pam50.genefu(datasets[[ii]])
    }

    ## compute claudin-low calls and add them to the clinical data frames
    for (ii in 1:length(datasets))
    {
        datasets[[ii]]$clinical$claudin.low <- huc.claudin.low(datasets[[ii]], data.dir)$claudin.low
    }

    ## add the lehmann subtypes if they are defined in the lehmann paper
    for (ii in 1:length(datasets))
    {
        datasets[[ii]]$clinical$lehmann <- huc.lehmann(datasets[[ii]], data.dir)
    }

    ## add the guedj CIT subtypes
    for (ii in 1:length(datasets))
    {
      datasets[[ii]]$clinical$cit <- huc.cit(datasets[[ii]], data.dir)
    }

    ## split up curtis into curtis.discovery* and curtis.validation*
    if (any(curtis.subsets))
    {
        stopifnot(identical(dataset.objects$curtis$clinical$METABRIC_ID, datasets$curtis$clinical$id.orig))
        stopifnot(all(dataset.objects$curtis$clinical$SET %in% c("DISCOVERY", "VALIDATION")))

        discovery.idx <- which(dataset.objects$curtis$clinical$SET == "DISCOVERY")
        validation.idx <- which(dataset.objects$curtis$clinical$SET == "VALIDATION")

        load(file.path(data.dir, "huc", "curtis-subset-ids.rda")) ## introduces 'curtis.subset.ids'
        stopifnot(setequal(c(curtis.subset.ids$curtis.discovery.1, curtis.subset.ids$curtis.discovery.2),
                           dataset.objects$curtis$clinical$METABRIC_ID[discovery.idx]))
        stopifnot(setequal(c(curtis.subset.ids$curtis.validation.1, curtis.subset.ids$curtis.validation.2),
                           dataset.objects$curtis$clinical$METABRIC_ID[validation.idx]))

        if ("curtis.validation" %in% dat.names)
        {
            datasets$curtis.validation <- huc.rm.patients(datasets$curtis, discovery.idx)
            datasets$curtis.validation$name <- "curtis.validation"
        }
        if ("curtis.discovery" %in% dat.names)
        {
            datasets$curtis.discovery <- huc.rm.patients(datasets$curtis, validation.idx)
            datasets$curtis.discovery$name <- "curtis.discovery"
        }
        if ("curtis.discovery.1" %in% dat.names)
        {
            rm.idx <- which(!datasets$curtis$clinical$id.orig %in% curtis.subset.ids$curtis.discovery.1)
            datasets$curtis.discovery.1 <- huc.rm.patients(datasets$curtis, rm.idx)
            datasets$curtis.discovery.1$name <- "curtis.discovery.1"
        }
        if ("curtis.discovery.2" %in% dat.names)
        {
            rm.idx <- which(!datasets$curtis$clinical$id.orig %in% curtis.subset.ids$curtis.discovery.2)
            datasets$curtis.discovery.2 <- huc.rm.patients(datasets$curtis, rm.idx)
            datasets$curtis.discovery.2$name <- "curtis.discovery.2"
        }
        if ("curtis.validation.1" %in% dat.names)
        {
            rm.idx <- which(!datasets$curtis$clinical$id.orig %in% curtis.subset.ids$curtis.validation.1)
            datasets$curtis.validation.1 <- huc.rm.patients(datasets$curtis, rm.idx)
            datasets$curtis.validation.1$name <- "curtis.validation.1"
        }
        if ("curtis.validation.2" %in% dat.names)
        {
            rm.idx <- which(!datasets$curtis$clinical$id.orig %in% curtis.subset.ids$curtis.validation.2)
            datasets$curtis.validation.2 <- huc.rm.patients(datasets$curtis, rm.idx)
            datasets$curtis.validation.2$name <- "curtis.validation.2"
        }
    }

    ## reorder the edatasets according to input
    datasets <- datasets[dat.names]

    ## make sure all data sets have the correct layout
    layout.check <- sapply(datasets, huc.check.dataset.layout)
    for (ii in which(!sapply(layout.check, isTRUE)))
    {
        warning(dat.names[[ii]], " does not conform to the standard data set layout.")
    }

    ## add the cohort variables
    for (ii in 1:length(datasets))
    {
        datasets[[ii]]$cohorts <- huc.cohorts(datasets[[ii]], data.dir)
    }

    ## clean up the data sets by removing patients that lack required clinical variables
    datasets <- lapply(datasets, huc.clean.dataset, rm.na.fields=rm.na.fields, rm.censored.5=rm.censored.5)

    ## add some additional notes to the data sets about what has changed
    ## since last code base.
    if ("anders" %in% names(datasets))
    {
        datasets$anders$notes <- c(datasets$anders$notes,
                                   "DIFF: removed the 4 samples lacking er status (were previously kept)")
    }
    if ("bild" %in% names(datasets))
    {
        datasets$bild$notes <- c(datasets$bild$notes, "DIFF: no changes")
    }
    if ("nki" %in% names(datasets))
    {
        datasets$nki$notes <- c(datasets$nki$notes,
                                "DIFF: old time variable referred to time to recurrence including local recurrence. Now, time refers to distant metastasis.")
    }
    if ("parker" %in% names(datasets))
    {
        datasets$parker$notes <- c(datasets$parker$notes,
                                "DIFF: 36 samples were removed in order to ensure only one array type is used.")
    }
    if ("sotiriou" %in% names(datasets))
    {
        datasets$sotiriou$notes <- c(datasets$sotiriou$notes,
                                "DIFF: time now refers to DMFS instead of RFS (which includes local recurrence) and because of this 5 patients were removed who lacked DMFS follow-up time. Another 5 patients lacking ER status were also removed.")
    }
    if ("vanvliet" %in% names(datasets))
    {
        datasets$vanvliet$notes <- c(datasets$vanvliet$notes,
                                     "DIFF: time and event have changed to better reflect our preference for distant metastasis. At 5 years, 41 patients have switched outcome to good. No patient switched to bad. 10 patients were also removed who did not have DMFS event information (their RFS event was previously being used).")
    }
    if ("wang" %in% names(datasets))
    {
        datasets$wang$notes <- c(datasets$wang$notes,
                                 "DIFF: no changes")
    }

    invisible(datasets)
}

### huc.get.exprs()
###
### returns the expression matrix for the specified data set and
### cohort. Genes and patients can also be specified.
###
### Arguments:
###
###     dat
###
###         data set with standard layout
###
###     cohort.name
###
###         name of a cohort
###
###     gene.names
###
###         character vector with gene names. If this is missing or is
###         NULL, then all genes will be used. See 'rm.missing.genes'
###         for handling of gene names that are not found in the data
###         set.
###
###     patient.ids
###
###         character vector of patient IDs (as specified by the column
###         names of the expression matrix, which is identical to the
###         rownames of the clinical variable). If not specified, all
###         patients in the cohorts will be used.
###
###     rm.missing.genes
###
###         logical indicating what to do with gene names that do not
###         appear in the data set. If TRUE, those genes will be ignored
###         and an expression matrix with possibly fewer rows than the
###         number of specified genes will be returned. If FALSE, the
###         expression matrix will contain rows with all NA values for
###         the genes that were missing in the data set.
###
###     show.warnings
###
###         logical indicating whether or not to show warning if all the
###         patients in 'patient.ids' do not belong to the specified
###         cohort.
###
### Returns:
###
###     An expression matrix, which depending on the input and
###     'rm.missing.genes' may have zero rows or columns. The rownames
###     of the expression matrix will correspond to genes and the
###     columns to patient IDs.
huc.get.exprs <- function(dat, cohort.name, gene.names=NULL, patient.ids=NULL, rm.missing.genes=TRUE, show.warnings=TRUE)
{
    stopifnot(is.character(cohort.name) && length(cohort.name) == 1 && cohort.name %in% names(dat$cohorts))

    probes <- dat$cohorts[[cohort.name]]$probes
    genes <- dat$probe.info$gene.name[probes]
    patients <- dat$cohorts[[cohort.name]]$patients

    stopifnot (is.null(patient.ids) || all(patient.ids %in% colnames(dat$exprs)))
    if (isTRUE(show.warnings) && !all(patient.ids %in% colnames(dat$exprs)[patients]))
        warning("requested patients are not all in the specified cohort")

    exprs <- dat$exprs[probes, ]
    rownames(exprs) <- genes
    if (!is.null(gene.names))
    {
        row.idx <- match(gene.names, genes)
        if (isTRUE(rm.missing.genes))
            row.idx <- na.omit(row.idx)

        exprs <- exprs[row.idx, , drop=FALSE]
        if (!isTRUE(rm.missing.genes))
          rownames(exprs) <- gene.names
    }

    if (is.null(patient.ids))
    {
        exprs <- exprs[, patients, drop=FALSE]
    }
    else
    {
        exprs <- exprs[, patient.ids, drop=FALSE]
    }

    return(exprs)
}

### huc.get.coxph.coef()
###
### returns the coxph coefficients for the specified data set and
### cohort. If specified, only the coefficients for the genes in
### 'gene.names' will be returned. This is a wrapper around
### huc.get.coxph(); see that function for description of the
### arguments.
###
### Returns:
###
###     A numeric vector with coxph coefficients. The names of the
###     vector correspond to gene names (see huc.get.coxph() for how
###     missing genes are handled).
huc.get.coxph.coef <- function(dat, cohort.name, gene.names=NULL, rm.missing.genes=TRUE)
{
    coxph <- huc.get.coxph(dat=dat, cohort.name=cohort.name, gene.names=gene.names, rm.missing.genes=rm.missing.genes)
    ret <- coxph$coef
    names(ret) <- rownames(coxph)
    return(ret)
}

### huc.get.coxph.pval()
###
### returns the coxph p-values for the specified data set and cohort. If
### specified, only the coefficients for the genes in 'gene.names' will
### be returned. This is a wrapper around huc.get.coxph(); see that
### function for description of the arguments.
###
### Returns:
###
###     A numeric vector with coxph p-values. The names of the
###     vector correspond to gene names (see huc.get.coxph() for how
###     missing genes are handled).
huc.get.coxph.pval <- function(dat, cohort.name, gene.names=NULL, rm.missing.genes=TRUE)
{
    coxph <- huc.get.coxph(dat=dat, cohort.name=cohort.name, gene.names=gene.names, rm.missing.genes=rm.missing.genes)
    ret <- coxph$pval
    names(ret) <- rownames(coxph)
    return(ret)
}

### huc.get.coxph()
###
### returns the coxph coefficients and p-values for the given data set
### and cohort. If specified, data only for the genes in 'gene.names'
### will be returned.
###
### Arguments:
###
###     dat
###
###         data set with standard layout
###
###     cohort.name
###
###         name of a cohort
###
###     gene.names
###
###         character vector with gene names. If this is missing or is
###         NULL, then all genes will be used. See 'rm.missing.genes'
###         for handling of gene names that are not found in the data
###         set.
###
###     rm.missing.genes
###
###         logical indicating what to do with gene names that do not
###         appear in the data set. If TRUE, gene names not found in the
###         data set will be ignored and a data frame with possibly
###         fewer rows than the number of specified genes will be
###         returned. If FALSE, the returned vector will contain NAs for
###         the genes that were missing in the data set.
###
### Returns:
###
###     A data frame with three variables: 'coef', 'pval', and
###     'row.idx'. If 'gene.names' is NULL, all genes in the cohort will
###     be used, otherwise, the selection will be restricted to the
###     given gene names in 'gene.names'. The rownames of the returned
###     data frame will be the gene names in the data set if
###     'gene.names' is NULL and 'gene.names' otherwise (if
###     'rm.missing.genes' is FALSE, rows corresponding to missing genes
###     will be given names from 'gene.names').
huc.get.coxph <- function(dat, cohort.name, gene.names=NULL, rm.missing.genes=TRUE)
{
    stopifnot(is.character(cohort.name) && length(cohort.name) == 1 && cohort.name %in% names(dat$cohorts))

    probes <- dat$cohorts[[cohort.name]]$probes
    genes <- dat$probe.info$gene.name[probes]
    patients <- dat$cohorts[[cohort.name]]$patients

    ret <- data.frame(coef=dat$cohorts[[cohort.name]]$coxph.coef[probes],
                      pval=dat$cohorts[[cohort.name]]$coxph.pval[probes],
                      row.idx=probes,
                      row.names=genes)

    if (!is.null(gene.names))
    {
        gene.names <- unique(gene.names)
        idx <- match(gene.names, rownames(ret))

        if (isTRUE(rm.missing.genes))
        {
            idx <- na.omit(idx)
            ret <- ret[idx, , drop=FALSE]
        }
        else
        {
            ret <- ret[idx, , drop=FALSE]
            rownames(ret)[is.na(idx)] <- gene.names[is.na(idx)]
        }
    }

    return(ret)
}

### huc.get.top.genes
###
### returns the top prognostic genes in the specified data set and
### cohort.
###
### Arguments:
###
###     dat
###
###         a data set with standard layout
###
###     cohort.name
###
###         name of a cohort
###
###     n
###
###         number of top genes to return. if 'n' is negative, then all
###         genes are returned, sorted by their coxph-pvalues.
###
### Returns:
###
###     a character verctor of gene names
huc.get.top.genes <- function(dat, cohort.name, n=100)
{
    stopifnot(n != 0)

    coxph <- huc.get.coxph(dat, cohort.name)
    ordered.genes <- rownames(coxph)[order(coxph$pval)]
    if (n < 0)
        return(ordered.genes)
    else
        return(ordered.genes[1:n])
}

### huc.get.clinical()
###
### returns the clinical data frame for specified data set and cohort.
###
### Arguments:
###
###     dat
###
###         a data set with standard layout
###
###     cohort.name
###
###         name of cohort
###
###     patient.ids
###
###         If NULL, then all patients in 'cohort.name' are
###         chosen. otherwise, patient ids are assumed to be a subset of
###         the patients in the specified cohort. So if you want to
###         freely decide which patients you want the clinical variables
###         for, set 'cohort.name' to "all".
###
### Returns:
###
###     a subset of the clinical data.frame in data set 'dat'.
huc.get.clinical <- function(dat, cohort.name="all", patient.ids=NULL)
{
    stopifnot(cohort.name %in% names(dat$cohorts))
    if (!is.null(patient.ids))
        stopifnot(all(patient.ids %in% dat$clinical$id[dat$cohorts[[cohort.name]]$patients]))

    if (is.null(patient.ids))
        patient.ids <- dat$cohorts[[cohort.name]]$patients

    dat$clinical[patient.ids, ]
}

### huc.get.cohort.ids()
###
### returns the patient ids in the specified data set and cohort.
###
### Arguments:
###
###     dat
###
###         a data set with standard layout
###
###     cohort.name
###
###         name of cohort
###
### Returns:
###
###     a character vector with patient ids.
huc.get.cohort.ids <- function(dat, cohort.name="all")
{
    stopifnot(cohort.name %in% names(dat$cohorts))

    dat$clinical$id[dat$cohorts[[cohort.name]]$patients]
}

### huc.rearrange.<data set name>()
###
### These functions convert the original data set objects to a standard
### layout. These functions are used internally by huc.load().
###
### Arguments:
###
###     dat:
###
###         the original dataset object loaded from file
###
###     data.dir:
###
###         path to the project data directory
###
### Returns:
###
###     A new dataset object with standard layout. The layout can and
###     should be checked via the huc.check.dataset.layout() function.
huc.rearrange.anders <- function(dat, data.dir)
{
    ## a few consistency checks on the original data object
    stopifnot(!any(duplicated(dat$clinical$PatientID)))

    ret <- list()

    ret$name <- "anders"

    ## exprs
    ret$exprs <- dat$exprdata
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$bioconductor.probes.info$AffyProbeID)
    probe.info$gene.name <- as.character(dat$bioconductor.probes.info$GeneName)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- dat$clinical$PatientID
    clinical$geo.sample <- as.character(dat$clinical$GEOSampleID)
    clinical$er <- dat$clinical$ESR == "positive"
    clinical$her2 <- dat$clinical$Her2byAmplicon == "positive"
    clinical$size <- as.double(dat$clinical$Tumor_Size_.cm.)
    clinical$stage <- dat$clinical$Path_Stg
    levels(clinical$stage) <- c(1, 2, 2, 3, 4, NA)
    clinical$stage <- as.integer(clinical$stage)
    clinical$grade <- as.integer(match(dat$clinical$GRADE,
                                       c("Well Diff", "Moderately Diff", "Poorly Diff", "Undifferentiated")))
    ## Note that lymph node negative is encoded as
    ## NA. 'Regional_Nodes_Examined' has no NAs indicating that every
    ## patient was examined for lymph node status.
    clinical$lymph <- ifelse(!is.na(dat$clinical$Regional_Nodes_Positive), TRUE, FALSE)
    clinical$age <- as.double(dat$clinical$Age_at_Dx)
    clinical$time <- as.double(dat$clinical$DFS)
    clinical$event <- dat$clinical$TYPE_1st_Recur %in% c("Distant", "DIstant")
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- dat$clinical$Chemotherapy == "Yes"
    clinical$tamoxifen <- dat$clinical$Hormone_Tx == "Yes"
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- as.character(c("Other", "Other", "IDC", "IDC/ILC","ILC", "Other", "IDC", "Other")[match(dat$clinical$Histo_Desc, c("Atypical medullary carcinoma","Carcinoma; NOS","Ductal Carcinoma; NOS", "Infiltrating Duct & Lobualar Carcinoma", "Lobular carcinoma; NOS", "Mucinous Adenocarcinoma","Paget Disease & Infiltrating Duct Carcinoma", "Papillary Carcinoma; NOS"))])

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "18167534"
    ret$geo.series <- "GSE7849"
    ret$platform.ensembl <- "affy_hg_u95av2"
    ret$platform.bioconductor <- "hgu95av2.db"
    ret$notes <- c("60 samples are ductal carcinomas and 14 samples are not (as indicated by clinical$Histo_Desc). should the non-DCs be removed?",
                   "her2 status was determined using the probes of the her2 amplicon genes",
                   "time and event correspond to distant metastasis (DFS)",
                   "4 samples lack er status")

    return(ret)
}

huc.rearrange.bild <- function(dat, data.dir)
{
    ## a few consistency checks on the original data object
    stopifnot(!any(duplicated(dat$clinical$PatientID)))

    ret <- list()

    ret$name <- "bild"

    ## exprs
    ret$exprs <- dat$exprdata
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$bioconductor.probes.info$AffyProbeID)
    probe.info$gene.name <- as.character(dat$bioconductor.probes.info$GeneName)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$clinical$PatientID)
    clinical$geo.sample <- as.character(dat$clinical$GEOSampleID)
    clinical$er <- dat$clinical$ERlevel >= 1
    clinical$her2 <- dat$clinical$Her2byAmplicon == "positive"
    clinical$size <- NA_real_
    clinical$stage <- NA_integer_
    clinical$grade <- NA_integer_
    clinical$lymph <- NA
    clinical$age <- NA_real_
    clinical$time <- as.double(dat$clinical$SurvivalTime.months.)
    clinical$event <- dat$clinical$DeadAliveStatus == "Dead"
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- NA
    clinical$tamoxifen <- NA
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- NA_character_

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "16273092"
    ret$geo.series <- "GSE3143"
    ret$platform.ensembl <- "affy_hg_u95av2"
    ret$platform.bioconductor <- "hgu95av2.db"
    ret$notes <- c("some expression values have been imputed",
                   "her2 status was determined using the probes of the her2 amplicon genes",
                   "time and event correspond to survival")

    return(ret)
}

huc.rearrange.boersma.epi <- function(dat, data.dir)
{
#    stopifnot(!any(duplicated(dat$clinical$Patient.ID)))

    ret <- list()

    ret$name <- "boersma.epi"

    ## exprs
    ret$exprs <- dat$exprs
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$bioconductor.probes.info$AffyProbeID)
    probe.info$gene.name <- as.character(dat$bioconductor.probes.info$GeneName)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$clinical$Patient.ID)
    clinical$geo.sample <- as.character(dat$clinical$GSMID)
    clinical$er <- dat$clinical$ESR == "Postiive"
    clinical$her2 <- dat$clinical$Her2 == "Positive"
    clinical$size <- NA_real_
    clinical$stage <- as.integer(c(1,2,2,2,3,3,4)[match(dat$clinical$Stage, c("I", "IIA", "IIB", "IIB ", "IIIA", "IIIB", "IV"))])
    clinical$grade <- NA_integer_
    clinical$lymph <- NA
    clinical$age <- NA_real_
    clinical$time <- NA_real_
    clinical$event <- dat$clinical$Death.Status == "Alive"
    clinical$time.5 <- NA_real_
    clinical$event.5 <- clinical$event

    clinical$chemo <- dat$clinical$Chemo == "Yes"
    clinical$tamoxifen <- NA
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- as.character(dat$clinical$IBC.Status)

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "17999412"
    ret$geo.series <- "GSE5847"
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- NA_character_
    ret$notes <- ""

    return(ret)
}

huc.rearrange.boersma.stroma <- function(dat, data.dir)
{
#    stopifnot(!any(duplicated(dat$clinical$Patient.ID)))

    ret <- list()

    ret$name <- "boersma.stroma"

    ## exprs
    ret$exprs <- dat$exprs
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$bioconductor.probes.info$AffyProbeID)
    probe.info$gene.name <- as.character(dat$bioconductor.probes.info$GeneName)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$clinical$Patient.ID)
    clinical$geo.sample <- as.character(dat$clinical$GSMID)
    clinical$er <- dat$clinical$ESR == "Postiive"
    clinical$her2 <- dat$clinical$Her2 == "Positive"
    clinical$size <- NA_real_
    clinical$stage <- as.integer(c(1,2,2,2,3,3,4)[match(dat$clinical$Stage, c("I", "IIA", "IIB", "IIB ", "IIIA", "IIIB", "IV"))])
    clinical$grade <- NA_integer_
    clinical$lymph <- NA
    clinical$age <- NA_real_
    clinical$time <- NA_real_
    clinical$event <- dat$clinical$Death.Status == "Alive"
    clinical$time.5 <- NA_real_
    clinical$event.5 <- clinical$event

    clinical$chemo <- dat$clinical$Chemo == "Yes"
    clinical$tamoxifen <- NA
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- as.character(dat$clinical$IBC.Status)

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "17999412"
    ret$geo.series <- "GSE5847"
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- NA_character_
    ret$notes <- ""

    return(ret)
}

huc.rearrange.curtis <- function(dat, data.dir)
{
    ## a few consistency checks on the original data object
    stopifnot(!any(duplicated(dat$clinical$METABRIC_ID)))

    ret <- list()

    ret$name <- "curtis"

    ## exprs
    ret$exprs <- dat$exprs
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## set NAs in the expression matrix to the mean of the rows
    na.idx <- which(is.na(ret$exprs), arr.ind=TRUE)
    if (nrow(na.idx) > 0)
    {
        row.means <- rowMeans(ret$exprs[na.idx[, 1], ], na.rm=TRUE)
        ret$exprs[na.idx] <- row.means
    }

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$probe.info$ID)
    probe.info$gene.name <- as.character(dat$probe.info$Symbol)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$clinical$METABRIC_ID)
    clinical$geo.sample <-  NA_character_
    clinical$er <- dat$clinical$ER.Expr == "+"
    clinical$her2 <- dat$clinical$Her2.Expr == "+"
    clinical$size <- as.double(ifelse(dat$clinical$size == "null", NA, dat$clinical$size)) / 10
    clinical$stage <- as.integer(ifelse(dat$clinical$stage == "null", NA, dat$clinical$stage))
    clinical$grade <- as.integer(ifelse(dat$clinical$grade == "null", NA, dat$clinical$grade))
    clinical$lymph <- as.integer(ifelse(dat$clinical$lymph_nodes_positive == "null", NA, dat$clinical$lymph_nodes_positive)) > 0
    clinical$age <- as.double(dat$clinical$age_at_diagnosis)
    clinical$time <- as.double(dat$clinical$T / 365 * 12)
    clinical$event <- dat$clinical$last_follow_up_status == "d-d.s."
    clinical$event[dat$clinical$last_follow_up_status == "d"] <- NA
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- rep(NA, length(clinical$event))
    clinical$event.5[dat$clinical$last_follow_up_status == "a"] <- FALSE
    clinical$event.5[dat$clinical$last_follow_up_status == "d" & clinical$time >= 60] <- FALSE
    clinical$event.5[dat$clinical$last_follow_up_status == "d" & clinical$time < 60] <- NA
    clinical$event.5[dat$clinical$last_follow_up_status == "d-d.s." & clinical$time >= 60] <- FALSE
    clinical$event.5[dat$clinical$last_follow_up_status == "d-d.s." & clinical$time < 60] <- TRUE
    clinical$event.5[dat$clinical$last_follow_up_status == "d-o.c."] <- FALSE

    n <- ncol(ret$exprs)
    clinical$chemo <- rep(FALSE, n)
    clinical$chemo[grep("CT", dat$clinical$Treatment)] <- TRUE
    clinical$tamoxifen <- rep(FALSE, n)
    clinical$tamoxifen[grep("HT", dat$clinical$Treatment)] <- TRUE
    clinical$herceptin <- NA

    clinical$intclust <- paste("intclust", dat$clinical$IntClustMemb, sep="")
    clinical$intclust[is.na(dat$clinical$IntClustMemb)] <- NA
    clinical$type <- c("Other", "DCIS", "IDC", "IDC", "IDC", "IDC", "IDC/ILC", "ILC", "Other", "Other", NA, "Other", "Other", "Other")[match(dat$clinical$histological_type, c("BENIGN", "DCIS", "IDC",
	"IDC-MED", "IDC-MUC", "IDC-TUB", "IDC+ILC", "ILC", "INVASIVE TUMOUR", "MIXED NST AND A SPECIAL TYPE", "null", "OTHER", "OTHER INVASIVE", "PHYL"))]

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "22522925"
    ret$geo.series <- NA_character_
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- NA_character_
    ret$notes <- c("Some patients also received RT (radiotherapy?) so untreated patients cannot currently be detected",
                   "ER status is by expression",
                   "Her2 status is by expression",
                   "pam50 is by their calls, not by inhouse method although that is indicated in the clinical variables")

    return(ret)
}

huc.rearrange.li <- function(dat, data.dir)
{
    stopifnot(!any(duplicated(dat$clinical$id)))

    ret <- list()

    ret$name <- "li"

    ## exprs
    ret$exprs <- dat$exprs
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$bioconductor.probes.info$AffyProbeID)
    probe.info$gene.name <- as.character(dat$bioconductor.probes.info$GeneName)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$clinical$id)
    clinical$geo.sample <- as.character(dat$clinical$geo)
    clinical$er <- c("pos", "pos", "neg")[match(dat$clinical$er, c("pos","pos-l", "neg"))] == "pos"
    clinical$her2 <- c("pos", "pos", "neg")[match(dat$clinical$her2, c("pos (3+)", "low pos (2+)", "neg"))] == "pos"
    clinical$size <- as.double(dat$clinical$size)
    clinical$stage <- NA_integer_
    clinical$grade <- as.integer(c(1,2,3)[match(dat$clinical$grade, c("I", "II", "III"))])
    clinical$lymph <- c(NA, "pos", "neg")[match(dat$clinical$ln, c("pos.micromet", "positive", "negative"))] == "pos"
    clinical$age <- as.numeric(dat$clinical$age)
    clinical$time <- as.double(dat$clinical$distant.recurrence.free.survival.months)
    clinical$event <- dat$clinical$distant.rec == "Y"
    clinical$time[!clinical$event] <- dat$clinical$time.of.followup.months[!clinical$event]
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & !is.na(clinical$time) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- c("chemo", "chemo", NA, "none", "chemo")[match(dat$clinical$chemo.class,  c("Anthracycline-based", "Trastuzumab", "Unknown", "none", "Other"))] == "chemo"
    clinical$tamoxifen <- c("none", "horm", "horm", "horm")[match(dat$clinical$hormonal.rx, c("none", "Tam", "unknown", "Arimidex"))] == "horm"
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- as.character(c("IDC", "ILC", "IDC/ILC")[match(dat$clinical$histology,c("Ductal", "Lobular", "Mixed"))])

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "20098429"
    ret$geo.series <- "GSE19615"
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- NA_character_
    ret$notes <- ""

    return(ret)
}

huc.rearrange.guedj <- function(dat, data.dir)
{
    ## a few consistency checks on the original data object
    stopifnot(!any(duplicated(dat$clinical$SampleID)))

    ret <- list()

    ret$name <- "guedj"

    ## exprs
    ret$exprs <- dat$exprs
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$bioconductor.probes.info$probes)
    probe.info$gene.name <- as.character(dat$bioconductor.probes.info$GeneName)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- dat$clinical$Source.Name
    clinical$geo.sample <- NA_character_
    clinical$er <- dat$clinical$Characteristics.ESR1..Protein.expression. == "Y"
    clinical$her2 <- dat$clinical$Characteristics.ERBB2..Protein.expression. == "Y"
    clinical$size <- NA_real_
    clinical$stage <- NA_integer_
    clinical$grade <- as.integer(dat$clinical$Characteristics.Grade..Scarff.Bloom.Richardson.)
    clinical$lymph <- dat$clinical$Characteristics.TNM..N. == "N1"
    clinical$age <- as.double(ifelse(dat$clinical$Characteristics.Age. == "ND", NA, dat$clinical$Characteristics.Age.))

    clinical$time <- as.double(dat$clinical$MFS.delay)
    clinical$time <- ifelse(clinical$time < 0, 0, clinical$time)
    clinical$event <- dat$clinical$MFS.event == 1
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & !is.na(clinical$time) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- NA
    clinical$tamoxifen <- NA
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- as.character(c("IDC", "ILC", "Other", "IDC/ILC", "Other")[match(dat$clinical$Characteristics.Histology., c("ductal", "lobular", "other", "ductal-locular", "micropapillary"))])

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "21785460"
    ret$geo.series <- NA_character_
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- "hgu133plus2"
    ret$notes <- c("")
    return(ret)
}

huc.rearrange.mcgill.epi <- function(dat, data, data.dir)
{
    stopifnot(!any(duplicated(dat$clinical$Patient)))

    ret <- list()

    ret$name <- "mcgill.epi"

    ## exprs
    ret$exprs <- dat$exprdata
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$probes.info$ID)
    probe.info$gene.name <- as.character(dat$probes.info$GeneName)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$clinical$Patient)
    clinical$geo.sample <- NA_character_
    clinical$er <- dat$clinical$ESR == "Positive"
    clinical$her2 <- dat$clinical$Her2 == "Positive"
    clinical$size <- as.numeric(dat$clinical$TumorSize) / 10
    clinical$stage <- c(1L, 2L, 2L, 3L, 3L, 3L)[match(dat$clinical$pathoStage, c("I", "IIA", "IIB", "IIIA", "IIIB", "IIIC"))]
    clinical$grade <- c(1L, 2L, 3L)[match(dat$clinical$GRADE,c("I","II","III"))]
    clinical$lymph <- dat$clinical$LN == "Positive"
    clinical$age <- as.double(dat$clinical$AgeOfOperation)

    event <- rep(FALSE, nrow(dat$clinical))
    event[dat$clinical$Recurrence.Type == "distant"] <- TRUE

    get.months <- function(ym) {
        sapply(strsplit(as.character(ym), "-"), function(x) {as.numeric(x[1]) * 12 + as.numeric(x[2])})
    }
    clinical$time <- with(dat$clinical, get.months(DateOfLastFollowUp) - get.months(DateOfOperation))
    rec.time <- with(dat$clinical, get.months(DateOfRecurrenceDiagnosis) - get.months(DateOfOperation))
    clinical$time[event] <- rec.time[event]

    clinical$event <- event

    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- NA
    clinical$tamoxifen <- NA
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- as.character(ifelse(dat$clinical$TT == "IDC-Epithelial", "IDC", NA))

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- NA_character_
    ret$geo.series <- NA_character_
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- NA_character_
    ret$notes <- ""

    return(ret)
}

huc.rearrange.mcgill.gq <- function(dat, data.dir)
{
  ## a few consistency checks on the original data object
  stopifnot(!any(duplicated(dat$clinical$serialNumber)))
  
  ret <- list()
  
  ret$name <- "mcgill.gq"
  
  dat$clinical <- dat$clinicals
  
  ## exprs
  ret$exprs <- dat$exprs
  rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
  # names making them unsuitable as
  # rownames of the expression data
  colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")
  
  ## probe.info
  probe.info <- list()
  probe.info$probe.id <- as.character(dat$probe.info$ID)
  probe.info$gene.name <- as.character(dat$probe.info$Symbol)
  probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
  probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
  ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)
  
  ## clinical
  clinical <- list()
  clinical$id <- colnames(ret$exprs)
  clinical$id.orig <- as.character(dat$clinical$serialNumber)
  clinical$geo.sample <- NA_character_
  clinical$er <- dat$clinical$ER == "Positive"
  clinical$her2 <- dat$clinical$Her2 == "Positive"
  clinical$size <- as.numeric(dat$clinical$T..Size) / 10
  clinical$stage <- c(1L, 2L, 2L, 3L, 3L, 4L)[match(dat$clinical$pathoStage, c("I", "IIA", "IIB", "IIIA", "IIIC", "IV"))]
  clinical$grade <- c(1L, 2L, 3L)[match(dat$clinical$GRADE,c("I","II","III"))]
  clinical$lymph <- dat$clinical$LN == "Positive"
  clinical$age <- as.double(dat$clinical$AgeSurgery)
  clinical$time <- suppressWarnings(as.double(dat$clinical$rfs.time.months))
  
  dat$clinical$rfs.event[is.na(clinical$time)] = NA
  
  clinical$time[clinical$time < 0] <- 0 ## patient MHA101, V04-510 
  clinical$event <- dat$clinical$rfs.event
  clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
  clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)
  
  clinical$chemo <- dat$clinical$chemotherapy == "Yes"
  clinical$tamoxifen <- dat$clinical$Tamoxifen.AI == "Yes"
  clinical$herceptin <- dat$clinical$Herceptin == "Yes"
  
  clinical$intclust <- NA_character_
  clinical$type <- NA_character_
  
  ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)
  
  ## misc
  ret$pmid <- NA_character_
  ret$geo.series <- NA_character_
  ret$platform.ensembl <- NA_character_
  ret$platform.bioconductor <- NA_character_
  ret$notes <- c("clinical information constantly being updated by Nick.  Patient 085 is actually HER2-negative.  There are non-tumor samples and technical replicates which are marked in the clinical data.  PAM50 inhouse was run on tumor samples (including technical replicates)")
  
  ## extra clinical
  extra.clinical <- list()
  extra.clinical$tissue.type <- dat$clinical$Type
  extra.clinical$technical.replicates <- dat$clinicals$technical.replicates
  extra.clinical$patient.id <- dat$clinical$patient.identifier
  extra.clinical$tumor.block <- dat$clinical$tumorBlock
  
  ret$extra.clinical <- as.data.frame(extra.clinical, row.names=clinical$id, stringsAsFactors=FALSE)
  
  return(ret)
}

huc.rearrange.mcgill.str <- function(dat, data, data.dir)
{
    stopifnot(!any(duplicated(dat$clinical$Patient)))

    ret <- list()

    ret$name <- "mcgill.str"

    ## exprs
    ret$exprs <- dat$exprdata
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$probes.info$ID)
    probe.info$gene.name <- as.character(dat$probes.info$GeneName)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$clinical$Patient)
    clinical$geo.sample <- NA_character_
    clinical$er <- dat$clinical$ESR == "Positive"
    clinical$her2 <- dat$clinical$Her2 == "Positive"
    clinical$size <- as.numeric(dat$clinical$TumorSize) / 10
    clinical$stage <- c(1L, 2L, 2L, 3L, 3L, 3L)[match(dat$clinical$pathoStage, c("I", "IIA", "IIB", "IIIA", "IIIB", "IIIC"))]
    clinical$grade <- c(1L, 2L, 3L)[match(dat$clinical$GRADE,c("I","II","III"))]
    clinical$lymph <- dat$clinical$LN == "Positive"
    clinical$age <- as.double(dat$clinical$AgeOfOperation)

    event <- rep(FALSE, nrow(dat$clinical))
    event[dat$clinical$Recurrence.Type == "distant"] <- TRUE

    get.months <- function(ym) {
        sapply(strsplit(as.character(ym), "-"), function(x) {as.numeric(x[1]) * 12 + as.numeric(x[2])})
    }
    clinical$time <- with(dat$clinical, get.months(DateOfLastFollowUp) - get.months(DateOfOperation))
    rec.time <- with(dat$clinical, get.months(DateOfRecurrenceDiagnosis) - get.months(DateOfOperation))
    clinical$time[event] <- rec.time[event]

    clinical$event <- event

    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- NA
    clinical$tamoxifen <- NA
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- as.character(ifelse(dat$clinical$TT == "IDC-Connective", "IDC", NA))

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- NA_character_
    ret$geo.series <- NA_character_
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- NA_character_
    ret$notes <- ""

    return(ret)
}

huc.rearrange.muggerud <- function(dat, data.dir)
{
    stopifnot(!any(duplicated(dat$clinical$Patient)))

    ret <- list()

    ret$name <- "muggerud"

    ## exprs
    ret$exprs <- dat$coeff
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$genes$PROBE_NAME)
    probe.info$gene.name <- as.character(dat$genes$GeneName)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$clinical$Patient)
    clinical$geo.sample <- NA_character_
    clinical$er <- dat$clinical$ER == 1
    clinical$her2 <- dat$clinical$Her2Amp == 1
    dat$clinical$Size[which(dat$clinical$Size == 999)] <- NA
    clinical$size <- as.double(dat$clinical$Size / 10)
    clinical$stage <- NA_integer_
    clinical$grade <- as.integer(as.vector(dat$clinical$Grade))
    clinical$lymph <- dat$clinical$Lymph == 1
    clinical$age <- as.double(dat$clinical$Age)
    clinical$time <- as.double(dat$clinical$DeathMonths)
    clinical$event <- dat$clinical$Death == "2"
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- dat$clinical$Chemo == 2
    clinical$tamoxifen <- dat$clinical$HormoneTreatment == 2
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- as.character(c("DCIS", "IDC", "DCIS/IDC")[match(dat$clinical$Diagnosis, c(1,2,3))])

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "20663721"
    ret$geo.series <- "GSE26304"
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- NA_character_
    ret$notes <- "time and event correspond to survival"

    return(ret)
}

huc.rearrange.nki <- function(dat, data.dir)
{
    ## a few consistency checks on the original data object
    stopifnot(!any(duplicated(dat$clinical$SampleID)))

    ret <- list()

    ret$name <- "nki"

    ## exprs
    ret$exprs <- dat$exprdata
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$rosetta.probes.info$RosettaProbeID)
    probe.info$gene.name <- as.character(dat$rosetta.probes.info$GeneName)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- dat$clinical$SampleID
    clinical$geo.sample <- NA_character_
    clinical$er <- dat$clinical$ESR == "positive"
    clinical$her2 <- dat$clinical$Her2byAmplicon == "positive"
    clinical$size <- as.double(dat$clinical$diameter.mm. / 10)
    clinical$stage <- NA_integer_
    clinical$grade <- as.integer(match(as.character(dat$clinical$Grade_3_classes),
                                       c("Well diff", "Intermediate", "Poorly diff")))
    clinical$lymph <- dat$clinical$LN == "positive"
    clinical$age <- as.double(dat$clinical$Age.years.)
    ## the time variable is the same as
    ## dat$clinical$Follow_up_time_or_metastasis, but by computing it
    ## this way, we get more digits (Follow_up_time_or_metastasis was
    ## rounded in the original data
    clinical$time <- dat$clinical$TIMEsurvival
    clinical$time[dat$clinical$EVENTmeta == "yes"] <- dat$clinical$TIMEmeta[dat$clinical$EVENTmeta == "yes"]
    clinical$time <- clinical$time * 12
    clinical$event <- dat$clinical$EVENTmeta == "yes"
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- dat$clinical$Chemo == "Yes"
    clinical$tamoxifen <- dat$clinical$Hormonal == "Yes"
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- "IDC"

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "12490681"
    ret$geo.series <- NA_character_
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- NA_character_
    ret$notes <- c("some expression values have been imputed",
                   "her2 status was determined using the probes of the her2 amplicon genes",
                   "time and event correspond to distant metastasis",
                   "these are all custom arrays")

    return(ret)
}

huc.rearrange.nowac.luma <- function(dat, data.dir)
{
    ## a few consistency checks on the original data object
    stopifnot(inherits(dat, "EList"))
    stopifnot(!any(duplicated(dat$clinical$SampleID)))

    ret <- list()

    ret$name <- "nowac.luma"

    ## exprs
    ret$exprs <- dat$E
    rownames(ret$exprs) <- NULL
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$genes$Ill_ID)
    probe.info$gene.name <- as.character(dat$genes$Symbol)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$targets$SampleNames)
    clinical$geo.sample <- NA_character_
    clinical$er <- dat$targets$ER == 12391 # factor level equals 12391
    clinical$her2 <- NA
    clinical$size <- NA_real_
    clinical$stage <- as.integer(dat$targets$Stage2)
    clinical$grade <- NA_integer_
    clinical$lymph <- dat$targets$LN == 1
    clinical$age <- NA_real_
    clinical$time <- dat$targets$FollowUp / 30
    clinical$event <- dat$targets$Event == 1
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- NA
    clinical$tamoxifen <- NA
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- "IDC"

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- NA_character_
    ret$geo.series <- NA_character_
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- NA_character_
    ret$notes <- c("")

    return(ret)
}

huc.rearrange.parker <- function(dat, data.dir)
{
    ## a few consistency checks on the original data object
    stopifnot(!any(duplicated(dat$clinical$Patient.ID)))

    ret <- list()

    ret$name <- "parker"

    ## exprs
    ret$exprs <- dat$exprdata
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$bioconductor.probes.info$AgilentProbeID)
    probe.info$gene.name <- as.character(dat$bioconductor.probes.info$GeneName)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$clinical$Patient.ID)
    clinical$geo.sample <- as.character(dat$clinical$GSMID)
    clinical$er <- dat$clinical$ESR == "positive"
    clinical$her2 <- dat$clinical$Her2byAmplicon == "positive"
    clinical$size <- NA_real_
    clinical$stage <- NA_integer_
    clinical$grade <- as.integer(dat$clinical$GRADE)
    clinical$lymph <- dat$clinical$LN == "positive"
    clinical$age <- as.double(dat$clinical$Age)
    clinical$time <- as.double(dat$clinical$Overall.survival.months)
    clinical$event <- dat$clinical$Overall.survival.event == "DOD or DOC"
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- NA
    clinical$tamoxifen <- NA
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- NA_character_

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "19204204"
    ret$geo.series <- "GSE10886"
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- NA_character_
    ret$notes <- c("mix of 3 different agilent arrays, one is a custom array. Kept only one type of array.",
                   "time and event correspond to survival (relapse info is ambiguous w.r.t. local/distant")

    ## We have to remove the following samples:
    ## - samples from normal tissue
    ## - samples from metastases
    ## - samples from a patient with two primaries
    ## - samples hybridized on custom array
    ## - samples marked as intrinsic pairs (these are from other perou publications)
    bad.idx <- unique(c(grep("recurrence|normal|custom|primaries", dat$clinical$Tissue.Source, ignore.case=TRUE),
                        grep("replicate|intrinsic", dat$clinical$Comment, ignore.case=T),
                        grep("met", dat$clinical$Patient.ID, ignore.case=T)))
    ret$exprs <- ret$exprs[, -bad.idx]
    ret$clinical <- ret$clinical[-bad.idx, ]

    ## we keep only the profiles with the GPL1390 platform. This info was downloaded
    ## from GEO using the GEOquery package:
    ##
    ## parker.platforms <- sapply(dat$clinical$GSMID, function(d) {Meta(getGEO(d))$platform_id})
    ## names(parker.platforms) <- dat$clinical$Patient.ID
    ##
    ## the above vector was saved in an object to speed things up.
    platforms <- get(load(file.path(data.dir, "huc", "parker-platforms.rda")))
    platforms <- platforms[match(ret$clinical$id.orig, names(platforms))]
    gpl1390.idx <- which(platforms %in% "GPL1390")
    ret$clinical <- ret$clinical[gpl1390.idx, ]
    ret$exprs <- ret$exprs[, gpl1390.idx]

    return(ret)
}

huc.rearrange.servant <- function(dat, data.dir)
{
    stopifnot(!any(duplicated(dat$clinical$GEO)))

    ret <- list()

    ret$name <- "servant"

    ## exprs
    ret$exprs <- dat$exprs
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$geo.probes.info$ID)
    probe.info$gene.name <- as.character(dat$geo.probes.info$Symbol)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$clinical$GEO)
    clinical$geo.sample <- NA_character_
    clinical$er <- dat$clinical$erbyge == 1
    clinical$her2 <- dat$clinical$herbyge == 1
    clinical$size <- NA_real_
    clinical$stage <- NA_integer_
    clinical$grade <- as.integer(c(1,2,3)[match(dat$clinical$histgrade,c("EEI","EEII","EEIII"))])
    clinical$lymph <- dat$clinical$pn == "pN+"
    clinical$age <- as.double(dat$clinical$age)
    clinical$time <- as.double(dat$clinical$rec.time * 12)
    ## no non-recurrent patient has time, but according to corresponding
    ## paper, all non-recurrent patients have follow-up time at least 10 years.
    clinical$time[is.na(clinical$time)] <- 120
    clinical$event <- dat$clinical$local.relapse == 1
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- dat$clinical$ct == "ChT"
    clinical$tamoxifen <- dat$clinical$ht == "HT"
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- as.character(c("IDC", "Other", "ILC")[match(dat$clinical$histtype, c("ductal", "other type", "lobular"))])

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "20663721"
    ret$geo.series <- "GSE26304"
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- NA_character_
    ret$notes <- ""

    return(ret)
}

huc.rearrange.sotiriou <- function(dat, data.dir)
{
    ## a few consistency checks on the original data object
    stopifnot(!any(duplicated(dat$clinical$sample_name)))

    ret <- list()

    ret$name <- "sotiriou"

    ## exprs
    ret$exprs <- dat$rma.exprdata
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$geo.probes.info$ID)
    probe.info$gene.name <- as.character(dat$geo.probes.info$Gene.Symbol)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$clinical$sample_name)
    clinical$geo.sample <- as.character(dat$clinical$geo_accn)
    clinical$er <- dat$clinical$ESR == "positive"
    clinical$her2 <- dat$clinical$Her2byAmplicon == "positive"
    clinical$size <- dat$clinical$size
    clinical$stage <- NA_integer_
    clinical$grade <- as.integer(dat$clinical$GRADE)
    clinical$lymph <- dat$clinical$LN == "positive"
    clinical$age <- as.double(dat$clinical$age)
    clinical$time <- as.double(dat$clinical$time.dmfs * 12)
    clinical$event <- dat$clinical$event.dmfs == "yes"
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- NA
    clinical$tamoxifen <- dat$clinical$adjuvant.treatment == "tamoxifen"
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- "IDC"

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## we keep only the oxford lab samples due to batch effects
    oxford.ids <- grep("^OX", ret$clinical$id.orig)
    ret$exprs <- ret$exprs[, oxford.ids]
    ret$clinical <- ret$clinical[oxford.ids, ]

    ## misc
    ret$pmid <- "16478745"
    ret$geo.series <- "GSE2990"
    ret$platform.ensembl <- "affy_hg_u133a"
    ret$platform.bioconductor <- "hgu133a.db"
    ret$notes <- c("only samples from one lab (oxford) are used due to batch affects",
                   "her2 status was determined using the probes of the her2 amplicon genes",
                   "time and event correspond to distant metastasis free survival (DMFS)")

    return(ret)
}

huc.rearrange.tcga <- function(dat, data.dir, tumor.only=TRUE)
{
    ##    stopifnot(!any(duplicated(dat$clinical$Patient.ID)))
    ret <- list()

    ret$name <- "tcga"

    ## exprs
    ret$exprs <- dat$exprs
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$probe.info$ReporterID)
    probe.info$gene.name <- as.character(dat$probe.info$genes)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- colnames(dat$exprs)
    clinical$geo.sample <- NA_character_
    clinical$er <- c("Negative", "Positive")[match(as.vector(dat$clinical$breast_carcinoma_estrogen_receptor_status), c("Negative", "Positive"))] == "Positive"
    clinical$her2 <- c("Negative", "Positive")[match(as.vector(dat$clinical$lab_proc_her2_neu_immunohistochemistry_receptor_status), c("Negative", "Positive"))] == "Positive"
    clinical$size <- NA_real_
    clinical$stage <- c(1L, 1L, 1L, 2L, 3L, 4L, 4L, 4L, NA)[match(as.vector(dat$clinical$breast_tumor_pathologic_t_stage, ), c("T1", "T1b", "T1c", "T2", "T3", "T4", "T4b", "T4d"))]
    clinical$grade <- NA_integer_
    dat$clinical$breast_tumor_pathologic_n_stage[which(dat$clinical$breast_tumor_pathologic_n_stage == "null" | dat$clinical$breast_tumor_pathologic_n_stage == "pNX")] <- NA
    clinical$lymph <-  ifelse(dat$clinical$breast_tumor_pathologic_n_stage == "pN0(i-)" | dat$clinical$breast_tumor_pathologic_n_stage == "pN0(i+)", 0, 1) == 1
    clinical$age <- as.vector(dat$clinical$age_at_initial_pathologic_diagnosis)
    clinical$age[clinical$age == "null"] <- NA
    clinical$age <- as.double(clinical$age)
    ## survival
    ## clinical$time <-  as.double(as.vector(dat$clinical$days_to_last_known_alive)) / (365/12)
    ## clinical$event <- c(1, 0, NA)[match(as.vector(dat$clinical$vital_status), c("DECEASED", "LIVING", "null"))] == 1
    ## recurrence
    clinical$time <- as.vector(dat$clinical$days_to_last_followup)
    clinical$time[clinical$time == "null"] <- NA
    clinical$time <- as.double(clinical$time) / (365/12)
    clinical$event <-  c(1, 0, NA)[match(as.vector(dat$clinical$breast_tumor_clinical_m_stage), c("M1", "M0", "null"))] == 1
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & !is.na(clinical$time) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- NA
    clinical$tamoxifen <- NA
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- as.character(c("IDC", "IDC", "ILC", "IDC/ILC", NA, "Other")[match(dat$clinical$breast_cancer_optical_measurement_histologic_type, c("Infiltrating Carcinoma, NOS", "Infiltrating Ductal", "Infiltrating Lobular", "Mixed Histology", "null", "Other"))])

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- NA_character_
    ret$geo.series <- NA_character_
    ret$platform.ensembl <- NA_character_
    ret$platform.bioconductor <- NA_character_
    ret$notes <- ""

    ## remove samples with NAs in their gene expression
    rm.idx <- which(colSums(is.na(ret$exprs)) != 0)

    ## remove samples from normal tissue
    if(isTRUE(tumor.only)) {
	rm.idx <- union(rm.idx, which(dat$clinical$tumour_type != "Tumour"))
    }

    if (length(rm.idx) > 0)
    {
        ret$exprs <- ret$exprs[, -rm.idx]
        ret$clinical <- ret$clinical[-rm.idx, ]
    }

    return(ret)
}

huc.rearrange.vanvliet <- function(dat, data.dir)
{
    ## a few consistency checks on the original data object
    stopifnot(!any(duplicated(dat$clinical$Arrays_ID)))

    ret <- list()

    ret$name <- "vanvliet"

    ## exprs
    ret$exprs <- dat$exprs
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$bioconductor.probes.info$AffyProbeID)
    probe.info$gene.name <- as.character(dat$bioconductor.probes.info$GeneName)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- dat$clinical$Arrays_ID
    clinical$geo.sample <- as.character(dat$clinical$Arrays_ID)
    clinical$geo.sample[grep("^GSM", clinical$geo.sample, invert=TRUE)] <- NA
    er <- get(load(file.path(data.dir, "huc", "vanvliet-er.rda")))
    clinical$er <- er[[1]] == "green"
    clinical$her2 <- dat$Her2ByAmplicon$HER2 == "positive"
    clinical$size <- as.double(dat$clinical$Size / 10)
    clinical$stage <- NA_integer_
    clinical$grade <- as.integer(dat$clinical$Grade)
    clinical$lymph <- as.logical(dat$clinical$Node)
    clinical$age <- as.double(dat$clinical$Age)
    ## time and event information is mixed in vanvliet depending on what
    ## is available in the underlying data sets.
    miller.idx <- which(dat$clinical$Set_1 == "miller")
    pawitan.idx <- which(dat$clinical$Set_1 == "paw")
    clinical$time <- as.double(dat$clinical$t_dmfs)
    clinical$time[miller.idx] <- dat$clinical$t_sos[miller.idx]
    clinical$time[pawitan.idx] <- pmin(dat$clinical$t_rfs, dat$clinical$t_os)[pawitan.idx]
    clinical$event <- as.logical(dat$clinical$e_dmfs)
    clinical$event[miller.idx] <- as.logical(dat$clinical$e_sos[miller.idx])
    clinical$event[pawitan.idx] <- as.logical(dat$clinical$e_rfs[pawitan.idx])
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- NA
    clinical$tamoxifen <- NA
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- NA_character_

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "18684329"
    ret$geo.series <- NA_character_
    ret$platform.ensembl <- "affy_hg_u133a"
    ret$platform.bioconductor <- "hgu133a.db"
    ret$notes <- c("The bioconductor package ROCR was used to predict ER status based on the expression of the probe '205225_at'",
                   "her2 status was determined using the probes of the her2 amplicon genes",
                   "time and event correspond to different things depending on original data set",
                   "    desmedt has RFS, DMFS, and OS. We use DMFS.",
                   "    miller has SOS only, which based on the paper, might mean disease-specific OS (death from BC)",
                   "    minn has MFS, but it is unclear if local metastasis is included in this definition",
                   "    pawitan has RFS, OS, and disease-specific OS, although from the paper it actually seems RFS is really DMFS. The survival times are strange because some RFS times are greater than OS. We use RFS",
                   "    loi has RFS and DFMS. We use DFMS",
                   "    chin has all four variables but RFS and DMFS are identical and so are OS and SOS. The events for RFS and DMFS are NOT identical however. we use DMFS and the corresponding event variable.")

    return(ret)
}

huc.rearrange.wang <- function(dat, data.dir)
{
    ## a few consistency checks on the original data object
    stopifnot(!any(duplicated(dat$clinical$GSMID)))

    ret <- list()

    ret$name <- "wang"

    ## exprs
    ret$exprs <- dat$exprdata
    rownames(ret$exprs) <- NULL # Most data sets have duplicated probe
                                # names making them unsuitable as
                                # rownames of the expression data
    colnames(ret$exprs) <- paste(ret$name, 1:ncol(ret$exprs), sep=".")

    ## probe.info
    probe.info <- list()
    probe.info$probe.id <- as.character(dat$geo.probes.info$ID)
    probe.info$gene.name <- as.character(dat$geo.probes.info$Gene.Symbol)
    probe.info$gene.name[probe.info$gene.name == "NA"] <- NA
    probe.info$gene.name[grep("^\\s*$", probe.info$gene.name)] <- NA
    ret$probe.info <- as.data.frame(probe.info, stringsAsFactors=FALSE)

    ## clinical
    clinical <- list()
    clinical$id <- colnames(ret$exprs)
    clinical$id.orig <- as.character(dat$clinical$GSMID)
    clinical$geo.sample <- as.character(dat$clinical$GSMID)
    clinical$er <- dat$clinical$ESR == "positive"
    clinical$her2 <- dat$clinical$Her2byAmplicon == "positive"
    clinical$size <- NA_real_
    clinical$stage <- NA_integer_
    clinical$grade <- NA_integer_
    clinical$lymph <- dat$clinical$LN == "positive"
    clinical$age <- NA_real_
    clinical$time <- as.double(dat$clinical$time.to.relapse)
    clinical$event <- dat$clinical$Relapse == "yes"
    clinical$time.5 <- ifelse(clinical$time >= 60, 60, clinical$time)
    clinical$event.5 <- ifelse(!is.na(clinical$event) & clinical$time >= 60, FALSE, clinical$event)

    clinical$chemo <- NA
    clinical$tamoxifen <- NA
    clinical$herceptin <- NA

    clinical$intclust <- NA_character_
    clinical$type <- "IDC"

    ret$clinical <- as.data.frame(clinical, row.names=clinical$id, stringsAsFactors=FALSE)

    ## misc
    ret$pmid <- "15721472"
    ret$geo.series <- "GSE2034"
    ret$platform.ensembl <- "affy_hg_u133a"
    ret$platform.bioconductor <- "hgu133a.db"
    ret$notes <- c("her2 status was determined using the probes of the her2 amplicon genes",
                   "time and event correspond to distant metastasis",
                   "all patients are lymph node negative")

    return(ret)
}

### huc.pam50.inhouse()
###
### given a (rearranged) data set object, this function will load the
### corresponding object with the in-house pam50 assignments and returns
### a character vector with standardized names. These names are the same
### names used by the genefu package.
###
### This function relies on the pam50 objects being data frames having
### rownames that correspond to the sample IDs in the data set object.
###
### Arguments:
###
###     dat
###
###         a rearranged data set object. Specifically, this function
###         relies on the 'dat$name' variable, which is used to obtain the
###         file name of the pam50 object to load.
###
###     data.dir
###
###         the path to the project data directory
###
### Returns:
###
###     A character vector with the same number of elements as rows in
###     dat$clinical.
huc.pam50.inhouse <- function(dat, data.dir)
{
    ## the nowac.luma data set is special because they are all LumA. No
    ## need to load anything.
    if (identical(dat$name, "nowac.luma"))
        return(rep("LumA", ncol(dat$exprs)))
  
    color2pam50 <- c(blue4="LumA",
                     deepskyblue="LumB",
                     firebrick2="Basal",
                     green4="Normal",
                     hotpink2="Her2")

    ## load the inhouse pam50 table and check for strange type.
    pam50 <- get(load(file.path(data.dir, "huc", paste(dat$name, "-pam50-inhouse.rda", sep=""))))

    stopifnot(is.data.frame(pam50) && identical(names(pam50), c("Subtype", "Cor")))
    stopifnot(all(pam50$Subtype %in% c(names(color2pam50), NA)))
    stopifnot(length(intersect(rownames(pam50), dat$clinical$id.orig)) > 0)

    pam50.colors <- pam50$Subtype[match(dat$clinical$id.orig, rownames(pam50))]
    return(color2pam50[pam50.colors])
}

### huc.pam50.genefu()
###
### Computes pam50 assignments based on the pam50 centroids provided by
### the pam50 authors (using the genefu package).
###
### NOTE: the genefu functions for assigning pam50 require Entrez gene
### IDs. we map genes to Entrez gene IDs based on their symbols, which
### is not ideal, but is the best we can do right now.
###
### NOTE: this function currently uses the pam50.scale model in the
### genefu package. Two other models are also available, but this is
### probably the best generic choice.
###
### Arguments:
###
###     dat
###
###         a rearranged data set object
###
### Returns:
###
###     A character vector with the same number of elements as rows in
###     dat$clinical.
huc.pam50.genefu <- function(dat)
{
    library(genefu)
    
    blood.subsets <- regexpr("^blood", dat$name) > 0
    if (any(blood.subsets))
    return(rep(NA_character_, ncol(dat$exprs))) 

    cmap <- pam50.scale$centroids.map
    idx <- which(dat$probe.info$gene.name %in% cmap$probe)
    exprs <- t(dat$exprs[idx, ])
    colnames(exprs) <- dat$probe.info$gene.name[idx]
    annot <- data.frame(EntrezGene.ID=cmap$EntrezGene.ID[match(colnames(exprs), cmap$probe)],
                        probe=colnames(exprs))
    unname(intrinsic.cluster.predict(sbt.model=pam50.scale, data=exprs, annot=annot, do.mapping=TRUE)$subtype)
}

## huc.claudin.low()
##
## Uses centroids to call claudin-low subtypes. The centroids used in
## this code are taken from the 9 cell lines dataset and using the
## procedure describe in the paper "Phenotypic and molecular
## characterization of the claudin-low intrinsic subtype of breast
## cancer" Aleix Prat, Joel S Parker, Olga Karginova, Cheng Fan, Chad
## Livasy, Jason I Herschkowitz, Xiaping He and Charles M
## Perou, Breast Cancer Research 2010, 12:R68 [PMID = 20813035].
##
## Function written by : Eric Paquet (eric.r.paquet@gmail.com)
##
## updated July 17th, 2012 We obtained the exact procedure from Prat and this is
##                         the new version using the new procedure. Mainly
##                         the sutyping is coming from DWD.
##
## Arguments:
##
##      dat
##
##          a data set with standard layout
##
##      data.dir
##
##          path to the top-level data directory of the project
##
## Returns:
##
##      list with two elements:
##
##          claudin.low
##
##              logical vector with the same length as the number of
##              samples in 'dat' indicating which samples are "claudin-low"
##
##          dist.matrix
##
##              A 2 x n numeric matrix, where n is the number of
##              patients in 'dat'. The rows correspond to distances to
##              the "claudin-low" and "others" centroids
##              respectively. The 'claudin.low' vector is determined
##              from 'dist.matrix'.
huc.claudin.low <- function(dat, data.dir)
{
  claudin.low.subtype = NA
  claudin.low.subtypes.dist = NA
  CLAUDIN.FOLDER <- file.path(data.dir,"claudin-low")

  if (all(c("name","exprs","clinical") %in% names(dat))){
    filename <- sprintf("claudin-%s.xls",toupper(dat$name))

    if (filename %in% list.files(CLAUDIN.FOLDER)){
      # read the claudin-low subtype obtained from DWD
      claudin.low <- read.delim(file.path(CLAUDIN.FOLDER,filename),header=T,stringsAsFactors=F,sep="\t")

      # we need to convert the sample names for matching
      claudin.low[,"Sample"] <- make.names(claudin.low[,"Sample"])
      colNames <- make.names(dat$clinical$id.orig)

      # we should be gettin a subtype for all samples
      stopifnot(all(colNames %in% claudin.low[,"Sample"]))

      # Transform claudin/others to yes/no
      claudin.low[,"Predict.type"] <- ifelse(claudin.low[,"Predict.type"] == "Claudin","Yes","No")

      match.sample.idx <- match(colNames,claudin.low[,"Sample"])

      claudin.low.subtype <- claudin.low[ match.sample.idx,"Predict.type"]

      claudin.low.subtypes.dist <-  claudin.low[match.sample.idx,2:3]
      rownames(claudin.low.subtypes.dist) <- colnames(dat$exprs)
    }
    else{
      cat(sprintf("!!!!\n!!!! Sorry, but we haven't assign claudin-low subtype for %s dataset\n!!!! Setting all sample to NA\n!!!!\n",dat$name))
      claudin.low.subtype <- rep(NA_character_, ncol(dat$exprs))
    }
  }
  else{
    cat("Claudin-low subtyping needs at least $exprs and the $name of the dataset\n")
  }

  list(claudin.low=claudin.low.subtype == "Yes",
         dist.matrix=claudin.low.subtypes.dist)
}

### huc.lehmann()
###
### returns a character vector with the lehmann-defined TN subtypes
### (BL1 BL2 IM LAR M MSL UNS). Only samples that were assigned in
### the original paper will be assigned (the actual method is not
### implemented).
###
### Arguments:
###
###     dat
###
###         a rearranged data set object
###
### Returns:
###
###     A character vector with the same number of elements as rows in
###     dat$clinical.
huc.lehmann <- function(dat, data.dir)
{
    lehmann <- read.csv(file.path(data.dir, "huc", "lehmann-subtypes-auth.csv"))
    ret <- as.character(lehmann$TNBC.SUBTYPE[match(dat$clinical$geo.sample, lehmann$SAMPLE)])
    ret[ret == "UNC"] <- NA
    unname(ret)
}

### huc.cit()
###
### returns a character vector with the guedj defined CIT subtypes
### "basL","lumA","lumB","normL","lumC","mApo".
### If data.dir is defined, then the pre-computed subtypes are supplied;
### otherwise, the subtypes are computed using the citbcmst R package.
### Datasets NKI, VanVliet, Parker and TCGA are assigned mostly basL and normL
### by the algorithm likely due to the normalization of these datasets.
### Consequently, subtype assignments for these datasets are omitted.
huc.cit <- function(dat, data.dir) {
  if (!missing(data.dir)) {
    cit <- read.csv(file.path(data.dir, "huc", "cit-subtypes.csv"))
    cit$citbcmst[which(cit$dataset %in% c("nki","vanvliet","parker","tcga"))] <- NA
    ret <- unname(as.character(cit$citbcmst[match(dat$clinical$id, cit$id)]))
    return(ret)
  }

  require(citbcmst)
  data(citbcmst)

  data <-data.frame(dat$exprs, stringsAsFactors=F)
  data.annot <- dat$probe.info
  column <- "Probe.Set.ID"
  if (length(intersect(dat$probe.info$probe.id, citbcmst$data.annot$Probe.Set.ID)) < 100) {
    genes <- unique(na.omit(unlist(dat$probe.info$gene.name)))
    i <- match(genes, dat$probe.info$gene.name)
    data <- data[i,]
    data.annot <- data.annot[i,]
    data.annot$probe.id <- genes
    column <- "Gene.Symbol"
  }
  rownames(data.annot) <- rownames(data) <- data.annot$probe.id
  ret <- cit.assignBcmst(data=data,
                       data.annot=data.annot,
                       data.colId="probe.id",
                       data.colMap="probe.id",
                       citbcmst.colMap=column)
  ret$dataset <- dat$name
  ret$id <- rownames(ret)
  ret
}


### huc.cohorts()
###
### returns a list containing information about each cohort for the
### given dataset.  This function saves its result into the cache/
### directory in the data directory of the project and computes the
### coxph coefficients and p-values only if the file is missing or if
### the definition of the cohorts has changed.
###
### Arguments:
###
###     dat:
###
###         rearranged dataset object
###
###     data.dir:
###
###         path to the project data directory
###
### Returns:
###
###     a list with each element corresponding to a cohort. Each element
###     is in turn a list with 4 items:
###
###     patients:
###
###         a vector of (integer) indices pointing out the patients that belong to the cohort
###
###     coxph.coef:
###
###         a numerical vector of coxph coefficients from the survival analysis
###
###     coxph.pval:
###
###         a numerical vector of p-values from the logrank test
###
###     probes:
###
###         a vector of (integer) indices pointing out the probes that
###         are used for this cohort. If a gene has multiple probes in a
###         dataset, the probe with the best coxph p-value is used
###         (which means different probes might be used in different
###         cohorts, even within the same dataset).
huc.cohorts <- function(dat, data.dir)
{
    ret <- list()
    cl <- dat$clinical

    ## define all the patient indices for all our cohorts
    ret$all$patients <- 1:nrow(dat$clinical)
    ret$erp$patients <- which(cl$er)
    ret$ern$patients <- which(!cl$er)
    ret$her2clinicalp$patients <- which(cl$her2)
    ret$her2clinicaln$patients <- which(!cl$her2)
    ret$her2clinicalp.erp$patients <- which(cl$er & cl$her2)
    ret$her2clinicalp.ern$patients <- which(!cl$er & cl$her2)
    ret$her2clinicaln.erp$patients <- which(cl$er & !cl$her2)
    ret$her2clinicaln.ern$patients <- which(!cl$er & !cl$her2)
    ret$luma$patients <- which(cl$pam50.parker == "LumA")
    ret$luma.erp$patients <- which(cl$pam50.parker == "LumA" & cl$er)
    ret$lumb$patients <- which(cl$pam50.parker == "LumB")
    ret$lumb.erp$patients <- which(cl$pam50.parker == "LumB" & cl$er)
    ret$normal$patients <- which(cl$pam50.parker == "Normal")
    ret$normal.erp$patients <- which(cl$pam50.parker == "Normal" & cl$er)
    ret$basal$patients <- which(cl$pam50.parker == "Basal")
    ret$basal.erp$patients <- which(cl$pam50.parker == "Basal" & cl$er)
    ret$basal.ern$patients <- which(cl$pam50.parker == "Basal" & !cl$er)
    ret$her2pam50p$patients <- which(cl$pam50.parker == "Her2")
    ret$her2pam50p.erp$patients <- which(cl$pam50.parker == "Her2" & cl$er)
    ret$her2pam50p.ern$patients <- which(cl$pam50.parker == "Her2" & !cl$er)
    ret$her2pam50p.her2clinicalp$patients <- which(cl$pam50.parker == "Her2" & cl$her2)
    ret$her2pam50p.her2clinicaln$patients <- which(cl$pam50.parker == "Her2" & !cl$her2)
    ret$her2pam50p.her2clinicalp.erp$patients <- which(cl$pam50.parker == "Her2" & cl$her2 & cl$er)
    ret$her2pam50p.her2clinicalp.ern$patients <- which(cl$pam50.parker == "Her2" & cl$her2 & !cl$er)
    ret$her2pam50p.her2clinicaln.erp$patients <- which(cl$pam50.parker == "Her2" & !cl$her2 & cl$er)
    ret$her2pam50p.her2clinicaln.ern$patients <- which(cl$pam50.parker == "Her2" & !cl$her2 & !cl$er)
    ret$her2pam50n$patients <- which(cl$pam50.parker != "Her2")
    ret$her2pam50n.erp$patients <- which(cl$pam50.parker != "Her2" & cl$er)
    ret$her2pam50n.ern$patients <- which(cl$pam50.parker != "Her2" & !cl$er)
    ret$intclust1$patients <- which(cl$intclust == "intclust1")
    ret$intclust2$patients <- which(cl$intclust == "intclust2")
    ret$intclust3$patients <- which(cl$intclust == "intclust3")
    ret$intclust4$patients <- which(cl$intclust == "intclust4")
    ret$intclust5$patients <- which(cl$intclust == "intclust5")
    ret$intclust6$patients <- which(cl$intclust == "intclust6")
    ret$intclust7$patients <- which(cl$intclust == "intclust7")
    ret$intclust8$patients <- which(cl$intclust == "intclust8")
    ret$intclust9$patients <- which(cl$intclust == "intclust9")
    ret$intclust10$patients <- which(cl$intclust == "intclust10")
    if (substr(dat$name, 1, 6) == "curtis")
    {
        ret$auth.luma$patients <- which(cl$pam50.authors == "LumA")
        ret$auth.lumb$patients <- which(cl$pam50.authors == "LumB")
        ret$auth.normal$patients <- which(cl$pam50.authors == "Normal")
        ret$auth.normal.erp$patients <- which(cl$pam50.authors == "Normal" & cl$er)
        ret$auth.basal$patients <- which(cl$pam50.authors == "Basal")
        ret$auth.basal.erp$patients <- which(cl$pam50.authors == "Basal" & cl$er)
        ret$auth.basal.ern$patients <- which(cl$pam50.authors == "Basal" & !cl$er)
        ret$auth.her2pam50p$patients <- which(cl$pam50.authors == "Her2")
        ret$auth.her2pam50p.erp$patients <- which(cl$pam50.authors == "Her2" & cl$er)
        ret$auth.her2pam50p.ern$patients <- which(cl$pam50.authors == "Her2" & !cl$er)
        ret$auth.her2pam50p.her2clinicalp$patients <- which(cl$pam50.authors == "Her2" & cl$her2)
        ret$auth.her2pam50p.her2clinicaln$patients <- which(cl$pam50.authors == "Her2" & !cl$her2)
        ret$auth.her2pam50p.her2clinicalp.erp$patients <- which(cl$pam50.authors == "Her2" & cl$her2 & cl$er)
        ret$auth.her2pam50p.her2clinicalp.ern$patients <- which(cl$pam50.authors == "Her2" & cl$her2 & !cl$er)
        ret$auth.her2pam50p.her2clinicaln.erp$patients <- which(cl$pam50.authors == "Her2" & !cl$her2 & cl$er)
        ret$auth.her2pam50p.her2clinicaln.ern$patients <- which(cl$pam50.authors == "Her2" & !cl$her2 & !cl$er)
        ret$auth.her2pam50n$patients <- which(cl$pam50.authors != "Her2")
        ret$auth.her2pam50n.erp$patients <- which(cl$pam50.authors != "Her2" & cl$er)
        ret$auth.her2pam50n.ern$patients <- which(cl$pam50.authors != "Her2" & !cl$er)
    }

    ## if a cached file exists, load it and compare the patient
    ## indices. If the patient indices defined above and in the file are
    ## the same, then return the loaded object. Otherwise, we need to
    ## compute and save the object.
    fname <- file.path(data.dir, "cache", paste(dat$name, "-cohorts.rda", sep=""))
    if (file.exists(fname))
    {
        cohorts <- get(load(fname))
        if (length(cohorts) == length(ret) &&
            setequal(names(cohorts), names(ret)) &&
            all(sapply(names(cohorts), function(cname) {
                identical(ret[[cname]]$patients, cohorts[[cname]]$patients)
            })))
        {
            return(cohorts)
        }
        else
        {
            stop(paste("huc.cohorts: bad cached cohorts file detected for data set", dat$name))
        }
    }
    cat("huc.cohorts: recomputing the cohort variables for dataset", dat$name, "\n")

    ## do some consistency checks
    partitioned.by <- function(x, ..., allow.missing=FALSE)
    {
        parts <- list(...)
        u <- do.call("c", parts)
        if (isTRUE(allow.missing))
            return(!any(duplicated(u)) && all(u %in% x))
        else
            return(!any(duplicated(u)) && setequal(u, x))
    }

    stopifnot(identical(ret$all$patients, 1:ncol(dat$exprs)))
    stopifnot(partitioned.by(ret$erp$patients,
                             ret$her2clinicalp.erp$patients, ret$her2clinicaln.erp$patients,
                             allow.missing=TRUE))
    stopifnot(partitioned.by(ret$her2clinicalp$patients,
                             ret$her2clinicalp.erp$patients, ret$her2clinicalp.ern$patients,
                             allow.missing=TRUE))
    stopifnot(partitioned.by(ret$all$patients,
                             ret$luma$patients, ret$lumb$patients, ret$normal$patients,
                             ret$basal$patients, ret$her2pam50p$patients, allow.missing=TRUE))
    stopifnot(all(ret$luma.erp$patients %in% ret$luma$patients))
    stopifnot(all(ret$lumb.erp$patients %in% ret$lumb$patients))
    stopifnot(all(ret$normal.erp$patients %in% ret$normal$patients))
    stopifnot(partitioned.by(ret$basal$patients,
                             ret$basal.erp$patients, ret$basal.ern$patients, allow.missing=TRUE))
    stopifnot(partitioned.by(ret$her2pam50p$patients,
                             ret$her2pam50p.erp$patients, ret$her2pam50p.ern$patients, allow.missing=TRUE))
    stopifnot(partitioned.by(ret$her2pam50p$patients,
                             ret$her2pam50p.her2clinicalp$patients, ret$her2pam50p.her2clinicaln$patients,
                             allow.missing=TRUE))
    stopifnot(partitioned.by(ret$her2pam50p$patients,
                             ret$her2pam50p.her2clinicalp.erp$patients, ret$her2pam50p.her2clinicaln.erp$patients,
                             ret$her2pam50p.her2clinicalp.ern$patients, ret$her2pam50p.her2clinicaln.ern$patients,
                             allow.missing=TRUE))
    stopifnot(partitioned.by(ret$her2pam50p.erp$patients,
                             ret$her2pam50p.her2clinicalp.erp$patients, ret$her2pam50p.her2clinicaln.erp$patients,
                             allow.missing=TRUE))
    stopifnot(partitioned.by(ret$her2pam50p.her2clinicalp$patients,
                             ret$her2pam50p.her2clinicalp.erp$patients, ret$her2pam50p.her2clinicalp.ern$patients,
                             allow.missing=TRUE))
    stopifnot(partitioned.by(ret$all$patients,
                             ret$her2pam50p$patients, ret$her2pam50n$patients, allow.missing=TRUE))
    stopifnot(partitioned.by(ret$all$patients,
                             ret$her2pam50p.erp$patients, ret$her2pam50n.erp$patients,
                             ret$her2pam50p.ern$patients, ret$her2pam50n.ern$patients, allow.missing=TRUE))
    stopifnot(partitioned.by(ret$erp$patients,
                             ret$her2pam50p.erp$patients, ret$her2pam50n.erp$patients, allow.missing=TRUE))
    stopifnot(partitioned.by(ret$her2pam50p$patients,
                             ret$her2pam50p.erp$patients, ret$her2pam50p.ern$patients, allow.missing=TRUE))

    ## get coxph coefficient and p-value for all cohorts
    library(survival)
    library(foreach)

    ret.names <- names(ret)
    ret <- foreach(cohort=ret) %do%
    {
        pat.idx <- cohort$patients

        if (length(pat.idx) == 0)
        {
            cohort$coxph.coef <- rep(NA_real_, nrow(dat$exprs))
            cohort$coxph.pval <- rep(NA_real_, nrow(dat$exprs))
        }
        else
        {
            surv <- Surv(dat$clinical$time.5[pat.idx], dat$clinical$event.5[pat.idx])

            library(parallel)
            coxph.stuff <- mclapply(1:nrow(dat$exprs), function(ii) {
                row <- dat$exprs[ii, pat.idx]
                x <- try(coxph(surv ~ row), silent = TRUE)
                if (class(x) == "try-error")
                    return(c(coxph.coef=NA_real_, coxph.pval=NA_real_))
                x <- summary(x)
                return(c(coxph.coef=x$coefficients[ , "coef"], coxph.pval=unname(x$logtest["pvalue"])))
            })
            coxph.stuff <- simplify2array(coxph.stuff)
            cohort$coxph.coef <- unname(coxph.stuff["coxph.coef", ])
            cohort$coxph.pval <- unname(coxph.stuff["coxph.pval", ])
        }

        ## for each gene, pick the best probe to use in the current cohort
        x <- data.frame(idx=1:nrow(dat$exprs), gene.name=dat$probe.info$gene.name, pval=cohort$coxph.pval)
        x <- x[!is.na(x$gene.name), ]
        x <- x[order(x$pval),]
        cohort$probes <- x$idx[which(!duplicated(x$gene.name))]
        cohort
    }
    names(ret) <- ret.names

    cohorts <- ret
    dir.create(dirname(fname), recursive=TRUE, showWarnings=FALSE)
    save(cohorts, file=fname)

    return(cohorts)
}

### huc.clean.dataset()
###
### given a rearranged data set object, this function removes those
### patients that do not have the required clinical variables.
###
### Arguments:
###
###     dat:
###
###         A dataset object as returned by one of the
###         huc.rearrange.<>() functions. The object is assumed to
###         conform to the standard layout.
###
###     rm.na.fields:
###
###         names of clinical variables that must not be missing
###
###     rm.censored.5
###
###         logical, indicating whethor or not patients that are good
###         outcome (i.e., event == FALSE) but whose follow-up time is
###         shorter than 5 years should be removed.
###
### Returns:
###
###         a data set object with the specified patients removed.
huc.clean.dataset <- function(dat, rm.na.fields, rm.censored.5)
{
    stopifnot(all(rm.na.fields %in% names(dat$clinical)))

    rm.idx <- integer(0)
    if (length(rm.na.fields) > 0)
        rm.idx <- c(rm.idx, which(rowSums(is.na(dat$clinical[rm.na.fields])) > 0))

    if (isTRUE(rm.censored.5))
        rm.idx <- c(rm.idx, which(dat$clinical$event.5 == FALSE & dat$clinical$time < 60))

    return(huc.rm.patients(dat, rm.idx))
}

### huc.rm.patients
###
### removes patients from a data set object. Importantly, this function
### handles patient indices in the 'cohort' variables correctly.
###
### Arguments:
###
###     dat:
###
###         A data set object with standard layout
###
###     patient.idx
###
###         indices to patients that should be removed.
###
### Returns:
###
###     A data set object with specified patients removed
huc.rm.patients <- function(dat, patient.idx)
{
    stopifnot(all(patient.idx %in% 1:ncol(dat$exprs)))

    if (length(patient.idx) == 0)
        return(dat)

    ## re-index the patients in the cohorts
    if ("cohorts" %in% names(dat))
    {
        dat$cohorts <- lapply(dat$cohorts, function(cohort) {
            patients.logical <- rep(FALSE, ncol(dat$exprs))
            patients.logical[cohort$patients] <- TRUE
            patients.logical <- patients.logical[-patient.idx]
            cohort$patients <- which(patients.logical)
            return(cohort)
        })
    }

    ## remove columns from dat$exprs and rows from dat$clinical
    dat$exprs <- dat$exprs[, -patient.idx]
    dat$clinical <- dat$clinica[-patient.idx, ]

    return(dat)
}

### huc.check.dataset.layout()
###
### Checks that dataset object layout
###
### Arguments:
###
###     dat:
###
###         dataset object returned from one of the huc.rearrange.<>()
###         functions.
###
### Returns:
###
###     TRUE if 'dat' conforms to the standard layout, and FALSE
###     otherwise.
huc.check.dataset.layout <- function(dat)
{
  required.names = c("name", "exprs", "probe.info", "clinical",
                     "pmid", "geo.series", "platform.ensembl", "platform.bioconductor", "notes")
  
  identical(class(dat), "list") &&
    all(names(dat) %in% c(required.names,"extra.clinical")) && 
    identical(required.names,names(dat)[1:length(required.names)]) && 
    identical(unname(sapply(dat[1:length(required.names)], class)),
              c("character", "matrix", "data.frame", "data.frame", rep("character", 5))) &&
    is.null(rownames(dat$exprs)) &&
    identical(nrow(dat$exprs), nrow(dat$probe.info)) &&
    identical(ncol(dat$exprs), nrow(dat$clinical)) &&
    !any(is.na(dat$exprs)) &&
    identical(rownames(dat$clinical), unname(dat$clinical$id)) &&
    identical(colnames(dat$exprs), rownames(dat$clinical)) &&
    identical(unname(sapply(dat$probe.info, class)), c("character", "character")) &&
    identical(names(dat$probe.info), c("probe.id", "gene.name")) &&
    ## !any(is.na(dat$probe.info$probe.id)) && # nki apparently has no information for 15 probes
    identical(names(dat$clinical), c("id", "id.orig", "geo.sample", "er", "her2", "size", "stage", "grade", "lymph", "age",
                                     "time", "event", "time.5", "event.5", "chemo", "tamoxifen", "herceptin", "intclust",
                                     "type", "pam50.authors", "pam50.inhouse", "pam50.genefu", "claudin.low", "lehmann", "cit")) &&
    identical(unname(sapply(dat$clinical, class)), c("character", "character",
                                                     "character", "logical", "logical", "numeric", "integer", "integer",
                                                     "logical", "numeric", "numeric", "logical", "numeric",
                                                     "logical", "logical", "logical", "logical", "character",
                                                     "character", "character", "character", "character", "logical", "character", "character")) &&
    identical(unname(sapply(dat$clinical, typeof)), c("character", "character",
                                                      "character", "logical", "logical", "double", "integer", "integer",
                                                      "logical", "double", "double", "logical", "double",
                                                      "logical", "logical", "logical", "logical", "character",
                                                      "character", "character", "character", "character", "logical", "character","character")) &&
    all(dat$clinical$stage %in% c(0:4, NA)) &&
    all(dat$clinical$grade %in% c(1:4, NA)) &&
    (all(is.na(dat$clinical$age)) || min(na.omit(dat$clinical$age)) > 0) &&
    (all(is.na(dat$clinical$age)) || max(na.omit(dat$clinical$age)) > 0) &&
    (all(is.na(dat$clinical$time)) || min(dat$clinical$time, na.rm=TRUE) >= 0) &&
    all(dat$clinical$time.5[which(dat$clinical$time >= 60)] == 60) &&
    all(dat$clinical$event.5[which(dat$clinical$time >= 60)] %in% c("FALSE", NA)) &&
    length(setdiff(which(is.na(dat$clinical$event.5)), which(is.na(dat$clinical$event)))) == 0 &&
    all(dat$clinical$intclust %in% c(paste("intclust", 1:10, sep=""), NA)) &&
    all(dat$clinical$type %in% c("DCIS", "DCIS/IDC", "IDC", "ILC", "IDC/ILC", "IBC", "non-IBC", "Other", NA)) &&
    all(dat$clinical$pam50.authors %in% c("LumA", "LumB", "Her2", "Basal", "Normal", NA)) &&
    all(dat$clinical$pam50.inhouse %in% c("LumA", "LumB", "Her2", "Basal", "Normal", NA)) &&
    all(dat$clinical$pam50.genefu %in% c("LumA", "LumB", "Her2", "Basal", "Normal", NA)) &&
    all(dat$clinical$lehmann %in% c("UNS", "BL1", "BL2", "IM", "M", "MSL", "LAR", NA)) &&
    all(dat$clinical$cit %in% c("basL","lumA","lumB","normL","lumC","mApo",NA)) &&
    identical(c(length(dat$pmid), length(dat$geo), length(dat$platform.ensembl), length(dat$platform.bioconductor)),
              c(1L, 1L, 1L, 1L))
}
### huc.dataset.info
###
### returns a character vector with general information about the given
### data set. Useful for creating tables. Do
###
### data.frame(lapply(huc, get.dataset.info))
###
### to get a data frame or
###
### sapply(huc, get.dataset.info)
###
### to get a character matrix
###
### Arguments:
###
###     dat
###
###         a data set with standard layout
###
###     readable.names
###
###         logical. If TRUE, then the names of the returned character
###         vector will be human readable. Otherwise, names are as
###         indicated in the code.
###
### Returns:
###
###     a character vector with information about the given data set.
huc.dataset.info <- function(dat, readable.names=TRUE) {
  huc.clinical.info(dat$clinical, readable.names=readable.names)
}

huc.clinical.info <- function(clinical, readable.names=TRUE) {
    samples <- nrow(clinical)
    stopifnot(all(c("size","grade","er","her2","lymph","pam50.inhouse","age","event","time","event.5","time.5","tamoxifen","chemo","herceptin") %in% colnames(clinical)))

    if(all(is.na(clinical$size))==TRUE) {
        size.2 <- NA
        size.2.plus <- NA
        size.mean <- NA
        size.SD <- NA
    } else {
        size.2 <- length(which(clinical$size <= 2))
        size.2.plus <- length(which(clinical$size > 2))
        size.mean <- format(mean(clinical$size, na.rm=T), digits=4)
        size.SD <- format(sd(clinical$size, na.rm=T), digits=4)

        s <- paste("(", format(size.2/samples * 100,digits=4), "%", ")", sep="")
        s1 <- paste("(", format(size.2.plus/samples * 100, digits=4), "%", ")", sep="")

        size.2 <- paste(size.2, s, sep=" ")
        size.2.plus <- paste(size.2.plus, s1, sep=" ")
    }

    if(all(is.na(clinical$grade))==TRUE) {
        grade.1 <- NA
        grade.2 <- NA
        grade.3 <- NA
    } else {
        grade.1 <- length(which(clinical$grade == 1))
        grade.2 <- length(which(clinical$grade == 2))
        grade.3 <- length(which(clinical$grade == 3))

        g <- paste("(", format(grade.1/samples * 100, digits=4), "%", ")", sep="")
        g1 <- paste("(", format(grade.2/samples * 100, digits=4), "%", ")", sep="")
        g2 <- paste("(", format(grade.3/samples * 100, digits=4), "%", ")", sep="")

        grade.1 <- paste(grade.1, g, sep=" ")
        grade.2	<- paste(grade.2, g1, sep=" ")
        grade.3	<- paste(grade.3, g2, sep=" ")
    }

    if(all(is.na(clinical$er))==TRUE) {
        er.p <- NA
        er.n <- NA
    } else {
        er.p <- length(which(clinical$er == TRUE))
        er.n <- length(which(clinical$er == FALSE))

        e <- paste("(", format(er.p/samples * 100, digits=4), "%", ")", sep="")
        e1 <- paste("(", format(er.n/samples * 100, digits=4), "%", ")", sep="")

        er.p <- paste(er.p, e, sep=" ")
        er.n <- paste(er.n, e1, sep=" ")
    }

    if(all(is.na(clinical$her2)) == TRUE){
        her2.p <- NA
        her2.n <- NA
    } else {
        her2.p <- length(which(clinical$her2 == TRUE))
        her2.n <- length(which(clinical$her2 == FALSE))

        h <- paste("(", format(her2.p/samples * 100, digits=4), "%", ")", sep="")
        h1 <- paste("(", format(her2.n/samples * 100, digits=4), "%", ")", sep="")

        her2.p <- paste(her2.p, h, sep=" ")
        her2.n <- paste(her2.n, h1, sep=" ")
    }

    if(all(is.na(clinical$lymph)) == TRUE){
        node.p <- NA
        node.n <- NA
    } else {
        node.p <- length(which(clinical$lymph == TRUE))
        node.n <- length(which(clinical$lymph == FALSE))

        n <- paste("(", format(node.p/samples * 100, digits=4), "%", ")", sep="")
        n1 <- paste("(", format(node.n/samples * 100, digits=4), "%", ")", sep="")

        node.p <- paste(node.p, n, sep=" ")
        node.n <- paste(node.n, n1, sep=" ")
    }

    pam50 <- sapply(c("LumA","LumB","Basal","Her2","Normal"), function(type) sum(clinical$pam50.genefu == type,na.rm=T))
    pam50 <- sapply(pam50, function(count) paste0(count, " (",format(count/samples*100,digits=4),"%)"))

    if (all(is.na(clinical$cit)))
    {
        cit <- rep(NA, length(bresect.cohorts.cit))
        names(cit) <- bresect.cohorts.cit
    }
    else
    {
        cit <- sapply(bresect.cohorts.cit, function(type) sum(clinical$cit == type,na.rm=T))
        cit <- sapply(cit, function(count) paste0(count, " (",format(count/samples*100,digits=4),"%)"))
    }

#     if (all(is.na(clinical$lehmann)))
#     {
#         lehmann <- rep(NA, length(bresect.cohorts.lehmann))
#         names(lehmann) <- bresect.cohorts.lehmann
#     }
#     else
#     {
#         lehmann <- sapply(bresect.cohorts.lehmann,function(type) sum(clinical$lehmann==type,na.rm=T))
#         lehmann <- sapply(lehmann,function(count) paste0(count, " (",format(count/samples*100,digits=4),"%)"))
#     }

    if(all(is.na(clinical$age)) == TRUE){
        age.50 <- NA
        age.50.plus <- NA
        age.mean <- NA
        age.range <- c(NA,NA)
    } else {
        age <- as.numeric(clinical$age)
        age.50 <- length(which(age < 50))
        age.50.plus <- length(which(age >= 50))
        age.mean <- format(mean(age, na.rm=T), digits=4)
        age.range <- range(age, na.rm=T)

        age.mean <- paste(age.mean, " (", age.range[1], "-", age.range[2], ")", sep="")

        a <- paste("(", format(age.50/samples * 100, digits=4), "%", ")", sep="")
        a1 <- paste("(", format(age.50.plus/samples * 100, digits=4), "%", ")", sep="")

        age.50 <- paste(age.50, a, sep=" ")
        age.50.plus <- paste(age.50.plus, a1, sep=" ")
    }

    if(all(is.na(clinical$weight)) == TRUE){
      weight.70 <- NA
      weight.70.plus <- NA
      weight.mean <- NA
      weight.range <- c(NA,NA)
    } else {
      weight <- as.numeric(clinical$weight)
      weight.70 <- length(which(weight <= 70))
      weight.70.plus <- length(which(weight > 70))
      weight.mean <- format(mean(weight, na.rm=T), digits=4)
      weight.range <- range(weight, na.rm=T)
      
      weight.mean <- paste(weight.mean, " (", weight.range[1], "-", weight.range[2], ")", sep="")
      
      a <- paste("(", format(weight.70/samples * 100, digits=4), "%", ")", sep="")
      a1 <- paste("(", format(weight.70.plus/samples * 100, digits=4), "%", ")", sep="")
      
      weight.70 <- paste(weight.70, a, sep=" ")
      weight.70.plus <- paste(weight.70.plus, a1, sep=" ")
    }

    if(all(is.na(clinical$MKS)) == TRUE){
      mks.6.5 <- NA
      mks.6.5.plus <- NA
      mks.mean <- NA
      mks.range <- c(NA,NA)
    } else {
      mks <- as.numeric(clinical$MKS)
      mks.6.5 <- length(which(mks <= 6.5))
      mks.6.5.plus <- length(which(mks > 6.5))
      mks.mean <- format(mean(mks, na.rm=T), digits=4)
      mks.range <- format(range(mks, na.rm=T), digits=4)
      
      mks.mean <- paste(mks.mean, " (", mks.range[1], "-", mks.range[2], ")", sep="")
      
      a <- paste("(", format(mks.6.5/samples * 100, digits=4), "%", ")", sep="")
      a1 <- paste("(", format(mks.6.5.plus/samples * 100, digits=4), "%", ")", sep="")
      
      mks.6.5 <- paste(mks.6.5, a, sep=" ")
      mks.6.5.plus <- paste(mks.6.5.plus, a1, sep=" ")
    }

    if(all(is.na(clinical$ERS)) == TRUE){
      ers.6.5 <- NA
      ers.6.5.plus <- NA
      ers.mean <- NA
      ers.range <- c(NA,NA)
    } else {
      ers <- as.numeric(clinical$ERS)
      ers.6.5 <- length(which(ers <= 6.5))
      ers.6.5.plus <- length(which(ers > 6.5))
      ers.mean <- format(mean(ers, na.rm=T), digits=4)
      ers.range <- format(range(ers, na.rm=T), digits=4)
      
      ers.mean <- paste(ers.mean, " (", ers.range[1], "-", ers.range[2], ")", sep="")
      
      a <- paste("(", format(ers.6.5/samples * 100, digits=4), "%", ")", sep="")
      a1 <- paste("(", format(ers.6.5.plus/samples * 100, digits=4), "%", ")", sep="")
      
      ers.6.5 <- paste(ers.6.5, a, sep=" ")
      ers.6.5.plus <- paste(ers.6.5.plus, a1, sep=" ")
    }

    if(all(is.na(clinical$HER2S)) == TRUE){
      her2s.7.5 <- NA
      her2s.7.5.plus <- NA
      her2s.mean <- NA
      her2s.range <- c(NA,NA)
    } else {
      her2s <- as.numeric(clinical$HER2S)
      her2s.7.5 <- length(which(her2s <= 7.5))
      her2s.7.5.plus <- length(which(her2s > 7.5))
      her2s.mean <- format(mean(her2s, na.rm=T), digits=4)
      her2s.range <- format(range(her2s, na.rm=T), digits=4)
      
      her2s.mean <- paste(her2s.mean, " (", her2s.range[1], "-", her2s.range[2], ")", sep="")
      
      a <- paste("(", format(her2s.7.5/samples * 100, digits=4), "%", ")", sep="")
      a1 <- paste("(", format(her2s.7.5.plus/samples * 100, digits=4), "%", ")", sep="")
      
      her2s.7.5 <- paste(her2s.7.5, a, sep=" ")
      her2s.7.5.plus <- paste(her2s.7.5.plus, a1, sep=" ")
    }

    if(all(is.na(clinical$event)) == TRUE){
        num.rec <- NA
    } else {
        num.rec <- format(length(which(clinical$event == TRUE)), digits=4)
    }

    if(all(is.na(clinical$time)) == TRUE){
        within.five.rec <- NA
        mean.rec <- NA
        range.rec <- c(NA, NA)
        mean.five.rec <- NA
        follow.up.rec <- NA
    } else {
        within.five.rec <- format(length(which(clinical$event.5)), digits=4)

        mean.five.rec <- format(mean(clinical$time.5, na.rm=T), digits=4)
        range.five.rec <- range(clinical$time.5, na.rm=T)

        mean.rec <- format(mean(clinical$time, na.rm=T), digits=2)
        range.rec <- range(clinical$time, na.rm=T)

        mean.rec <- paste(mean.rec, " (", format(range.rec[1], digits=4), "-", format(range.rec[2],digits=4), ")", sep="")
        mean.five.rec <- paste(mean.five.rec, " (", format(range.five.rec[1], digits=4), "-", format(range.five.rec[2],digits=4), ")", sep="")
    }

    if(all(is.na(clinical$tamoxifen)) == TRUE){
        tamoxifen.y <- NA
        tamoxifen.n <- NA
    } else {
        tamoxifen.y  <- length(which(clinical$tamoxifen))
        tamoxifen.n  <- length(which(!clinical$tamoxifen))
    }

    if(all(is.na(clinical$chemo)) == TRUE){
        chemo.y <- NA
        chemo.n <- NA
    } else {
        chemo.y  <- length(which(clinical$chemo))
        chemo.n  <- length(which(!clinical$chemo))
    }

    if(all(is.na(clinical$herceptin)) == TRUE){
        herceptin.y <- NA
        herceptin.n <- NA
    } else {
        herceptin.y  <- length(which(clinical$herceptin))
        herceptin.n  <- length(which(!clinical$herceptin))
    }

    all <- c(samples=samples, size.2 = size.2, size.2.plus = size.2.plus, size.mean = size.mean, size.SD=size.SD, grade.1=grade.1,
             grade.2=grade.2, grade.3=grade.3, er.p=er.p, er.n=er.n, her2.p=her2.p, her2.n=her2.n,
             node.p=node.p, node.n=node.n,
             pam50.la=pam50[["LumA"]], pam50.lb=pam50[["LumB"]], pam50.b=pam50[["Basal"]],
             pam50.h=pam50[["Her2"]], pam50.n=pam50[["Normal"]],
             cit=cit, #lehmann=lehmann, ## here is an R lesson for you, how does this work?
             mks.6.5=mks.6.5, mks.6.5.plus=mks.6.5.plus, mks.mean=mks.mean,
             ers.6.5=ers.6.5, ers.6.5.plus=ers.6.5.plus, ers.mean=ers.mean,
             her2s.7.5=her2s.7.5, her2s.7.5.plus=her2s.7.5.plus, her2s.mean=her2s.mean,
             age.50=age.50, age.50.plus=age.50.plus, age.mean=age.mean,
             weight.70=weight.70, weight.70.plus=weight.70.plus, weight.mean=weight.mean,
             num.rec=num.rec, within.five.rec=within.five.rec, mean.rec=mean.rec, mean.five.rec=mean.five.rec,
             tamoxifen.y = tamoxifen.y, tamoxifen.n = tamoxifen.n, chemo.y = chemo.y, chemo.n = chemo.n, herceptin.y = herceptin.y, herceptin.n = herceptin.n)

    if (isTRUE(readable.names))
    {
        name.map <- c(samples="No. Samples",
                      size.2="Size (< 20mm)",
                      size.2.plus="Size (> 20mm)",
                      size.mean="Size: mean",
                      size.SD="Size: SD",
                      grade.1="Grade 1",
                      grade.2="Grade 2",
                      grade.3="Grade 3",
                      er.p="ER positive",
                      er.n="ER negative",
                      her2.p="HER2 positive",
                      her2.n="HER2 negative",
                      node.p="Lymph node positive",
                      node.n="Lymph node negative",
                      pam50.la="PAM50 lumA",
                      pam50.lb="PAM50 lumB",
                      pam50.b="PAM50 basal",
                      pam50.h="PAM50 her2",
                      pam50.n="PAM50 normal",
                      cit=sapply(names(cit), function(type) paste("CIT", type)), ## lesson in R pt. 2!
                      #lehmann=sapply(names(lehmann), function(type) paste("Lehmann", type)),
                      mks.6.5="Mitosis kinase score (< 6.5)",
                      mks.6.5.plus="Mitosis kinase score (> 6.5)",
                      mks.mean="Mitosis kinase score: mean (range)",
                      ers.6.5="Estrogen receptor score (< 6.5)",
                      ers.6.5.plus="Estrogen receptor score (> 6.5)",
                      ers.mean="Estrogen receptor score: mean (range)",
                      her2s.7.5="HER2 score (< 7.5)",
                      her2s.7.5.plus="HER2 score (> 7.5)",
                      her2s.mean="Estrogen receptor: mean (range)",
                      age.50="Age (< 50)",
                      age.50.plus="Age (> 50)",
                      age.mean="Age: mean (range)",
                      weight.70="Weight (< 70)",
                      weight.70.plus="weight (> 70)",
                      weight.mean="Weight: mean (range)",
                      num.rec="Total relapse",
                      within.five.rec="Total relapse (within five years)",
                      mean.rec="Total relapse: mean in months (range)",
                      mean.five.rec="Total relapse: mean in months within five years (range)",
                      tamoxifen.y="No Patients receiving Tamoxifen",
                      tamoxifen.n="No Patients NOT receiving Tamoxifen",
                      chemo.y="No Patients receiving Chemotherpy",
                      chemo.n="No Patients NOT receiving Chemotherapy",
                      herceptin.y="No Patients receiving Herception",
                      herceptin.n="No Patients NOT receiving Herceptin")

        stopifnot(setequal(names(name.map), names(all)))
        stopifnot(!any(duplicated(names(name.map))))

        names(all) <- name.map[names(all)]
    }
    return(all)
}
