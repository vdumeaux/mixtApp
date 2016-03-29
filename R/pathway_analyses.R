
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
  weight01Fisher <- runTest(GO.genes, algorithm = "weight01", statistic = "fisher")
  print(weight01Fisher)
#   resultKS <- runTest(GO.genes, algorithm = "classic", statistic = "ks")
#   print(resultKS)
#   resultKS.elim <- runTest(GO.genes, algorithm = "elim", statistic = "ks")
#   print(resultKS.elim) 
  
  allRes <- GenTable(GO.genes, classicFisher = resultFisher, weight01Fisher= weight01Fisher, orderBy = "weight01Fisher", topNodes=length(usedGO(GO.genes))) 
  
  return (list(GO.data=GO.genes, GO.table=allRes, resultFisher=resultFisher, weight01Fisher= weight01Fisher))

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
    file = "../../data/h.all.v5.1.symbols.gmt"
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
    file = "../../data/c2.all.v5.1.symbols.gmt"
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
    file = "../../data/c2.cp.v5.1.symbols.gmt"
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
    file = "../../data/c2.cgp.v5.1.symbols.gmt"
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
    MSigDB.c2.cgp[[current.line]] <- line.vector
    names(MSigDB.c2.cgp)[current.line] = name
    
    current.line = current.line + 1
  }
  
  close(con)
  
  return (MSigDB.c2.cgp)
} 

load.MSigDB.C3 <- function(file=NULL) {
  if (is.null(file)) {
    file = "../../data/c3.all.v5.1.symbols.gmt"
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
    file = "../../data/c5.all.v5.1.symbols.gmt"
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
    file = "../../data/c6.all.v5.1.symbols.gmt"
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
    file = "../../data/c7.all.v5.1.symbols.gmt"
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



