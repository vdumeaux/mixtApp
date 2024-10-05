# Fetch Stroma Datasets

dir = getwd()

# Datasets we can include directly.
load(file.path(dir,"data/bresat.rda"))
# devtools::use_data(bresat,  overwrite=TRUE, compress="gzip")

load(paste0(dir,"/data/moduleColors.rda"))
# devtools::use_data(moduleColors,  overwrite=TRUE, compress="gzip")

load(paste0(dir,"/data/msigdb.enrichment.rda"))
# msigdb.enrichment = msigdb_enrichment
# devtools::use_data(msigdb.enrichment,  overwrite=TRUE, compress="gzip")

load(paste0(dir,"/data/dat.rda"))
# devtools::use_data(dat,  overwrite=TRUE, compress="gzip")

load(paste0(dir,"/data/mod_clinical_fdr.rda"))
# devtools::use_data(mod_clinical_fdr,  overwrite=TRUE, compress="gzip")

load(paste0(dir,"/data/perm_cor_p.rda"))
# devtools::use_data(perm_cor_p,  overwrite=TRUE, compress="gzip")

load(paste0(dir,"/data/net.rda"))
# net = TOM_net
# devtools::use_data(net, overwrite=TRUE, compress="gzip")

# Datasets we need to wrangle a bit to include
load(paste0(dir,"/data/goterms.rda"))
# devtools::use_data(goterms, overwrite=TRUE, compress="gzip")