### function to run limma for subtypes and create heatmap

heatmap.limma<-function(dat, n=n, subtypes=c("er", "her2","LumA","LumB", "Normal", "Her2","Basal", 
                                             "luma.erp", "lumb.erp", "basal.erp", "normal.erp", "her2pam50p.erp", "her2clinicaln.ern", "her2clinicalp.ern"),
                        scheme=c("pam50.genefu","pam50.inhouse", "cit", "hybrid.genefu"), output.name=NULL) {
  
  if (length(scheme)>1) stop("you must specify only one genomic subtyping scheme")
  if(any(table(dat$clinical[,match(scheme, colnames(dat$clinical))])<5)){
    subtype.excluded<-names(table(dat$clinical[,match(scheme, colnames(dat$clinical))]))[table(dat$clinical[,match(scheme, colnames(dat$clinical))])<5]
    print(paste(subtype.excluded, "includes less than 5 patients", sep=" "))
    subtypes<-subtypes[!table(dat$clinical[,match(scheme, colnames(dat$clinical))])<5]
  } 
  
  gene.table<-vector("list", length(subtypes)+1)
  names(gene.table)<-c(subtypes,"union")
  pop<-vector("list", length(subtypes))
  
  if (!is.null(output.name)) {
    pdf(file.path("output", paste(output.name, "pdf", sep=".")))}
  
  for (s in 1:length(subtypes)){
    if ("er" %in% subtypes[s]){
      pop[[s]]$pop1<-colnames(dat$exprs)[dat$clinical[,match("er", colnames(dat$clinical))]==T]
      pop[[s]]$pop2<-colnames(dat$exprs)[dat$clinical[,match("er", colnames(dat$clinical))]==F]}
    else if ("her2" %in% subtypes[s]){
      pop[[s]]$pop1<-colnames(dat$exprs)[dat$clinical[,match("her2", colnames(dat$clinical))]==T]
      pop[[s]]$pop2<-colnames(dat$exprs)[dat$clinical[,match("her2", colnames(dat$clinical))]==F]}
    else {
      pop[[s]]$pop1<-colnames(dat$exprs)[dat$clinical[,match(scheme, colnames(dat$clinical))]==subtypes[s]]
      pop[[s]]$pop2<-colnames(dat$exprs)[!dat$clinical[,match(scheme, colnames(dat$clinical))]==subtypes[s]]}
    
    gene.table[[s]]<-topTable(limma.simple(dat=dat, pop1=pop[[s]]$pop1, pop2=pop[[s]]$pop2), n=nrow(dat$exprs))
    gene.table[[s]]<-gene.table[[s]][,colnames(gene.table[[s]]) %in% c("row.idx","gene.name", "logFC", "P.Value", "adj.P.Val")]
    gene.table[[s]]$rank<-1:nrow(dat$exprs)
    gene.table[[(length(subtypes))+1]]<-unique(c(gene.table[[(length(subtypes))+1]],gene.table[[s]]$gene.name[1:n]))
    
    if (!is.null(output.name)){
      genes<-gene.table[[s]]$gene.name[1:n]
      exprs<-dat$exprs[gene.table[[s]]$row.idx[1:n],]
      bresat <- sig.ranksum(exprs, ns=1:nrow(exprs), full.return=TRUE)
      cc.heatmap.clinical <- huc.color.clinical(clinical=dat$clinical)
      ## order expression matrix and clinical info
      exprs <- exprs[bresat$gene.order, bresat$pat.order, drop=FALSE]
      clinical <- cc.heatmap.clinical[bresat$pat.order, ]
      plot.new()
      heatmap.simple(exprs,
                     row.labels=rownames(exprs),
                     row.clust=FALSE,
                     col.clust=FALSE,
                     title=paste(print(subtypes[s]), "markers in blood", sep=" "),
                     clinical=clinical)
    }  
  }
  if (!is.null(output.name)){
    exprs<-dat$exprs[rownames(dat$exprs) %in% gene.table$union,]
    bresat <- sig.ranksum(exprs, ns=1:nrow(exprs), full.return=TRUE)
    cc.heatmap.clinical <- huc.color.clinical(clinical=idc$clinical)
    ## order expression matrix and clinical info
    exprs<- exprs[bresat$gene.order, bresat$pat.order, drop=FALSE]
    clinical <- cc.heatmap.clinical[bresat$pat.order, ]
    grid.newpage()
    heatmap.simple(exprs,
                   row.labels=rownames(exprs),
                   row.clust=FALSE,
                   col.clust=FALSE,
                   title="union all markers in blood",
                   clinical=clinical)}
  
  return(gene.table)
}



### function to compute go enrichment analyses

go.enrichment<- function (gene.list, bckg.genes){
  library (topGO)
  #library (lumiHumanAll.db)
  
  xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
  head(xx)
  xx <- unique(unlist(xx))
  xx<-xx[xx %in% bckg.genes]
  
  gene.list<- factor(as.integer(xx %in% gene.list))
  names(gene.list) <- xx
  
  GO.genes <- new("topGOdata",
                  ontology = "BP",
                  allGenes = gene.list,
                  nodeSize = 5,
                  annot = annFUN.org, 
                  mapping = "org.Hs.eg.db",
                  ID = "symbol")
  
  resultFisher <- runTest(GO.genes, algorithm = "classic", statistic = "fisher")
  print(resultFisher)
  resultKS <- runTest(GO.genes, algorithm = "classic", statistic = "ks")
  print(resultKS)
  resultKS.elim <- runTest(GO.genes, algorithm = "elim", statistic = "ks")
  print(resultKS.elim) 
  
  allRes <- GenTable(GO.genes, classicFisher = resultFisher, elimKS = resultKS.elim,
                     orderBy = "classicFisher", ranksOf = "elimKS", topNodes = 50) 
  
  return (list(GO.data=GO.genes, GO.table=allRes, resultFisher=resultFisher))

  #printGraph(GO.genes, resultFisher, firstSigNodes=20, fn.prefix=group, useInfo="all", pdfSW=TRUE)
}



load.MSigDB.diag <- function(file=NULL) {
  if (is.null(file)) {
    file = "../../data/diag_sig.gmt"
  }
  inputFile <- file
  con  <- file(inputFile, open = "r")
  
  MSigDB.diag <- list()
  current.line = 1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line.vector <- strsplit(line, "\t")[[1]]
    name = line.vector[1]
    line <- paste(c(paste(line.vector[3:(length(line.vector)-1)],'\t',sep=""),line.vector[length(line.vector)]),collapse="")
    line.vector <- strsplit(gsub("/"," ",gsub("\t"," ", line,fixed=TRUE),fixed=TRUE)," ")[[1]]
    line.vector = line.vector[line.vector != ""]
    MSigDB.diag[[current.line]] <- line.vector
    names(MSigDB.diag)[current.line] = name
    
    current.line = current.line + 1
  }
  
  close(con)
  
  return (MSigDB.diag)
}  

load.MSigDB.I <- function(file=NULL) {
  if (is.null(file)) {
    file = "../../data/chaussabel_sets.gmt"
  }
  inputFile <- file
  con  <- file(inputFile, open = "r")
  
  MSigDB.i <- list()
  current.line = 1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line.vector <- strsplit(line, "\t")[[1]]
    name = line.vector[1]
    line <- paste(c(paste(line.vector[3:(length(line.vector)-1)],'\t',sep=""),line.vector[length(line.vector)]),collapse="")
    line.vector <- strsplit(gsub("/"," ",gsub("\t"," ", line,fixed=TRUE),fixed=TRUE)," ")[[1]]
    line.vector = line.vector[line.vector != ""]
    MSigDB.i[[current.line]] <- line.vector
    names(MSigDB.i)[current.line] = name
    
    current.line = current.line + 1
  }
  
  close(con)
  
  return (MSigDB.i)
}  

load.MSigDB.H <- function(file=NULL) {
  if (is.null(file)) {
    file = "../../data/h.all.v5.0.symbols.gmt"
  }
  
  inputFile <- file
  con  <- file(inputFile, open = "r")
  
  MSigDB.h <- list()
  current.line = 1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line.vector <- strsplit(line, "\t")[[1]]
    name = line.vector[1]
    line <- paste(c(paste(line.vector[3:(length(line.vector)-1)],'\t',sep=""),line.vector[length(line.vector)]),collapse="")
    line.vector <- strsplit(gsub("/"," ",gsub("\t"," ", line,fixed=TRUE),fixed=TRUE)," ")[[1]]
    line.vector = line.vector[line.vector != ""]
    MSigDB.h[[current.line]] <- line.vector
    names(MSigDB.h)[current.line] = name
    
    current.line = current.line + 1
  }
  
  close(con)
  
  return (MSigDB.h)
}

load.MSigDB.C2 <- function(file=NULL) {
  if (is.null(file)) {
    file = "../../data/c2.all.v5.0.symbols.gmt"
  }
  inputFile <- file
  con  <- file(inputFile, open = "r")
  
  MSigDB.c2 <- list()
  current.line = 1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line.vector <- strsplit(line, "\t")[[1]]
    name = line.vector[1]
    line <- paste(c(paste(line.vector[3:(length(line.vector)-1)],'\t',sep=""),line.vector[length(line.vector)]),collapse="")
    line.vector <- strsplit(gsub("/"," ",gsub("\t"," ", line,fixed=TRUE),fixed=TRUE)," ")[[1]]
    line.vector = line.vector[line.vector != ""]
    MSigDB.c2[[current.line]] <- line.vector
    names(MSigDB.c2)[current.line] = name
    
    current.line = current.line + 1
  }
  
  close(con)
  
  return (MSigDB.c2)
}  

load.MSigDB.C2.cp <- function(file=NULL) {
  if (is.null(file)) {
    file = "../../data/c2.cp.v5.0.symbols.gmt"
  }
  inputFile <- file
  con  <- file(inputFile, open = "r")
  
  MSigDB.c2.cp <- list()
  current.line = 1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line.vector <- strsplit(line, "\t")[[1]]
    name = line.vector[1]
    line <- paste(c(paste(line.vector[3:(length(line.vector)-1)],'\t',sep=""),line.vector[length(line.vector)]),collapse="")
    line.vector <- strsplit(gsub("/"," ",gsub("\t"," ", line,fixed=TRUE),fixed=TRUE)," ")[[1]]
    line.vector = line.vector[line.vector != ""]
    MSigDB.c2.cp[[current.line]] <- line.vector
    names(MSigDB.c2.cp)[current.line] = name
    
    current.line = current.line + 1
  }
  
  close(con)
  
  return (MSigDB.c2.cp)
}  

load.MSigDB.C2.cgp <- function(file=NULL) {
  if (is.null(file)) {
    file = "../../data/c2.cgp.v5.0.symbols.gmt"
  }
  inputFile <- file
  con  <- file(inputFile, open = "r")
  
  MSigDB.c2.cgp <- list()
  current.line = 1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line.vector <- strsplit(line, "\t")[[1]]
    name = line.vector[1]
    line <- paste(c(paste(line.vector[3:(length(line.vector)-1)],'\t',sep=""),line.vector[length(line.vector)]),collapse="")
    line.vector <- strsplit(gsub("/"," ",gsub("\t"," ", line,fixed=TRUE),fixed=TRUE)," ")[[1]]
    line.vector = line.vector[line.vector != ""]
    MSigDB.c2.cp[[current.line]] <- line.vector
    names(MSigDB.c2.cp)[current.line] = name
    
    current.line = current.line + 1
  }
  
  close(con)
  
  return (MSigDB.c2.cgp)
} 

load.MSigDB.C3 <- function(file=NULL) {
  if (is.null(file)) {
    file = "../../data/c3.all.v5.0.symbols.gmt"
  }
  inputFile <- file
  con  <- file(inputFile, open = "r")
  
  MSigDB.c3 <- list()
  current.line = 1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line.vector <- strsplit(line, "\t")[[1]]
    name = line.vector[1]
    line <- paste(c(paste(line.vector[3:(length(line.vector)-1)],'\t',sep=""),line.vector[length(line.vector)]),collapse="")
    line.vector <- strsplit(gsub("/"," ",gsub("\t"," ", line,fixed=TRUE),fixed=TRUE)," ")[[1]]
    line.vector = line.vector[line.vector != ""]
    MSigDB.c3[[current.line]] <- line.vector
    names(MSigDB.c3)[current.line] = name
    
    current.line = current.line + 1
  }
  
  close(con)
  
  return (MSigDB.c3)
}  

load.MSigDB.C5 <- function(file=NULL) {
  if (is.null(file)) {
    file = "../../data/c5.all.v5.0.symbols.gmt"
  }
  
  inputFile <- file
  con  <- file(inputFile, open = "r")
  
  MSigDB.c5 <- list()
  current.line = 1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line.vector <- strsplit(line, "\t")[[1]]
    name = line.vector[1]
    line <- paste(c(paste(line.vector[3:(length(line.vector)-1)],'\t',sep=""),line.vector[length(line.vector)]),collapse="")
    line.vector <- strsplit(gsub("/"," ",gsub("\t"," ", line,fixed=TRUE),fixed=TRUE)," ")[[1]]
    line.vector = line.vector[line.vector != ""]
    MSigDB.c5[[current.line]] <- line.vector
    names(MSigDB.c5)[current.line] = name
    
    current.line = current.line + 1
  }
  
  close(con)
  
  return (MSigDB.c5)
}

load.MSigDB.C6 <- function(file=NULL) {
  if (is.null(file)) {
    file = "../../data/c6.all.v5.0.symbols.gmt"
  }
  
  inputFile <- file
  con  <- file(inputFile, open = "r")
  
  MSigDB.c6 <- list()
  current.line = 1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line.vector <- strsplit(line, "\t")[[1]]
    name = line.vector[1]
    line <- paste(c(paste(line.vector[3:(length(line.vector)-1)],'\t',sep=""),line.vector[length(line.vector)]),collapse="")
    line.vector <- strsplit(gsub("/"," ",gsub("\t"," ", line,fixed=TRUE),fixed=TRUE)," ")[[1]]
    line.vector = line.vector[line.vector != ""]
    MSigDB.c6[[current.line]] <- line.vector
    names(MSigDB.c6)[current.line] = name
    
    current.line = current.line + 1
  }
  
  close(con)
  
  return (MSigDB.c6)
}

load.MSigDB.C7 <- function(file=NULL) {
  if (is.null(file)) {
    file = "../../data/c7.all.v5.0.symbols.gmt"
  }
  
  inputFile <- file
  con  <- file(inputFile, open = "r")
  
  MSigDB.c7 <- list()
  current.line = 1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line.vector <- strsplit(line, "\t")[[1]]
    name = line.vector[1]
    line <- paste(c(paste(line.vector[3:(length(line.vector)-1)],'\t',sep=""),line.vector[length(line.vector)]),collapse="")
    line.vector <- strsplit(gsub("/"," ",gsub("\t"," ", line,fixed=TRUE),fixed=TRUE)," ")[[1]]
    line.vector = line.vector[line.vector != ""]
    MSigDB.c7[[current.line]] <- line.vector
    names(MSigDB.c7)[current.line] = name
    
    current.line = current.line + 1
  }
  
  close(con)
  
  return (MSigDB.c7)
}



