# Script to import datasets from the larger MIxT Repository into this R package.

# Directory where the MIxT repository is cloned. Note that you may need to run
# git raw fix within the mixt_dir/data/mixt folder to load the datasets.
mixt_dir = "~/mixt"

# Datasets we can include directly.
load(paste0(mixt_dir,"/data/mixt/bresat.RData"))
devtools::use_data(bresat,  overwrite=TRUE, compress="gzip")

load(paste0(mixt_dir,"/data/mixt/moduleColors.RData"))
devtools::use_data(moduleColors,  overwrite=TRUE, compress="gzip")

load(paste0(mixt_dir,"/data/mixt/msigdb.RData"))
devtools::use_data(msigdb.enrichment,  overwrite=TRUE, compress="gzip")

load(paste0(mixt_dir,"/data/mixt/combat_data.RData"))
devtools::use_data(dat,  overwrite=TRUE, compress="gzip")

load(paste0(mixt_dir,"/data/mixt/mod_clinical_fdr.RData"))
# the mod_clinical_fdr file contains an object named fdr not mod_clinical_fdr
# so we'll create a new object with the correct name
mod_clinical_fdr = fdr
devtools::use_data(mod_clinical_fdr,  overwrite=TRUE, compress="gzip")

load(paste0(mixt_dir,"/data/mixt/perm_cor_p.RData"))
devtools::use_data(perm.cor.p,  overwrite=TRUE, compress="gzip")

# Datasets we need to wrangle a bit to include
load(paste0(mixt_dir,"/data/mixt/TOM.RData"))

tom<-NULL
for (tissue in c("blood","biopsy", "nblood")){
  tom[[tissue]] <- TOM[[tissue]][moduleColors[[tissue]] != "grey", moduleColors[[tissue]] !="grey"]
}

tom$nblood<- TOM$nblood[moduleColors$blood != "grey", moduleColors$blood !="grey"]

net<-NULL
for (tissue in c("blood", "biopsy")){
  net[[tissue]]<-WGCNA::exportNetworkToCytoscape(
    tom[[tissue]],
    edgeFile = NULL, #paste("data/", tissue, "_smod_tom_01_edge.txt", sep=""),
    nodeFile = NULL, #paste("data/", tissue, "_smod_tom_01_node.txt", sep=""),
    weighted = TRUE,
    threshold = 0.1,
    nodeAttr = moduleColors[[tissue]][moduleColors[[tissue]] !="grey"])
}

devtools::use_data(net, overwrite=TRUE, compress="gzip")


load(paste0(mixt_dir,"/data/mixt/topGO_mod.RData"))
load(paste0(mixt_dir,"/data/mixt/go_common.RData"))
goterms <- all.single
for (tissue in names(goterms)){
  for(module in names(goterms[[tissue]])){
    goterms[[tissue]][[module]]$GO.data <- NULL
    goterms[[tissue]][[module]]$resultFisher <- NULL
    goterms[[tissue]][[module]]$common <- NULL
    goterms[[tissue]][[module]]$common <- go.common[[tissue]][[module]]
  }
}

devtools::use_data(goterms, overwrite=TRUE, compress="gzip")