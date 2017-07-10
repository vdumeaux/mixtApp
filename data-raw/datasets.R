# Fetch Stroma Datasets

dir = 'https://rawgit.com/vdumeaux/mixtData/master'

# Datasets we can include directly.
load(url(paste0(dir,"/data/bresat.RData")))
devtools::use_data(bresat,  overwrite=TRUE, compress="gzip")

load(url(paste0(dir,"/data/moduleColors.RData")))
devtools::use_data(moduleColors,  overwrite=TRUE, compress="gzip")

load(url(paste0(dir,"/data/msigdb_enrichment.RData")))
msigdb.enrichment = msigdb_enrichment
devtools::use_data(msigdb.enrichment,  overwrite=TRUE, compress="gzip")

load(url(paste0(dir,"/data/dat.RData")))
devtools::use_data(dat,  overwrite=TRUE, compress="gzip")

load(url(paste0(dir,"/data/mod_clinical_fdr.RData")))
# the mod_clinical_fdr file contains an object named fdr not mod_clinical_fdr
# so we'll create a new object with the correct name
mod_clinical_fdr = fdr
devtools::use_data(mod_clinical_fdr,  overwrite=TRUE, compress="gzip")

load(url(paste0(dir,"/data/perm_cor_p.RData")))
perm.cor.p = perm_cor_p
devtools::use_data(perm.cor.p,  overwrite=TRUE, compress="gzip")

load(url(paste0(dir,"/data/TOM_net.RData")))
net = TOM_net
devtools::use_data(net, overwrite=TRUE, compress="gzip")

# Datasets we need to wrangle a bit to include
load(url(paste0(dir,"/data/goterms.RData")))
devtools::use_data(goterms, overwrite=TRUE, compress="gzip")