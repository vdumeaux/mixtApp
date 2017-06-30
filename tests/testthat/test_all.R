library(testthat)

# Note: We do not check for correctness, just that all the functions run without
# crashing. The rational behind this is that the data in the package is
# continuously changing which would make it necessary to update the tests 
# whenever we updated the data.

test_that("heatmap works", {
  tissue = "blood"
  module = "blue"
  expect_silent(cohort_heatmap(tissue, module)) 
  expect_silent(cohort_heatmap(tissue,module,"erp"))
})


test_that("scatterplot works", {
  x.tissue = "blood"
  y.tissue = "biopsy"
  x.module = "blue"
  y.module = "blue"
  expect_silent(cohort_scatterplot(x.tissue,x.module, y.tissue, y.module))
})

test_that("boxplot works", {
  blood.module="blue"
  orderByTissue="blood"
  orderByModule="blue"
  expect_silent(cohort_boxplot(blood.module, orderByTissue, orderByModule))
  orderByTissue="biopsy"
  expect_silent(cohort_boxplot(blood.module, orderByTissue, orderByModule))
})

test_that("get modules works", {
  expect_silent(getModules("blood"))
  expect_silent(getModules("biopsy"))
  expect_silent(getModules("bnblood"))
  expect_silent(getModules("nblood"))
})

test_that("get all genes works", {
  expect_silent(getAllGenes())
})

test_that("get all modules works", {
  expect_silent(getAllModules("BRCA1"))
})

# test_that("get all genes and modules works", {
#   expect_silent(getAllGenesAndModules())
# })

test_that("get all tissues works", {
  expect_equal(getAllTissues(), c("blood","biopsy","nblood","bnblood"))
})

test_that("get gene list works", {
  expect_silent(getGeneList("blood", "blue"))
})

test_that("get en richment scores works", {
  expect_silent(getEnrichmentScores("blood","green")) 
  expect_silent(getEnrichmentScores("biopsy","green")) 
})

test_that("get gene set names works", {
  expect_silent(getGeneSetNames())
})

test_that("get nerichment for tissue works", {
  expect_silent(getEnrichmentForTissue("blood"))
  expect_silent(getEnrichmentForTissue("biopsy"))
})

test_that("get go terms works", {
  expect_silent(getGOTerms("blood","green"))
  expect_silent(getGOTerms("biopsy","green"))
})

test_that("get common genes works", {
  expect_silent(getCommonGenes("blood","green", "REACTOME_RNA_POL_I_TRANSCRIPTION"))
})

test_that("get common go terms works", {
  expect_silent(getCommonGOTermGenes("blood","green","GO:0070848"))
})

test_that("user enrichment scores works", {
  expect_silent(userEnrichmentScores("blood",c("BRCA1", "ESR1")))
})

test_that("common enrichment score genes works", {
  expect_silent(commonEnrichmentScoreGenes("blood", "green", c("BRCA1", "ESR1")))
})

test_that("get go term names works ", {
  expect_silent(getGOTermNames())
})

test_that("get go scores for tissue work", {
  expect_silent(getGOScoresForTissue("biopsy", "B cell receptor signaling pathway"))
})

test_that("gene overlap test works", {
  expect_silent(geneOverlapTest("blood", "biopsy"))
})


test_that("patient ranksum works", {
  expect_silent(patientRankSum("blood","biopsy", "all"))
  expect_silent(patientRankSum("biopsy","blood", "all"))
})

test_that("comparison analyses work", {
  expect_silent(comparisonAnalyses("blood", "biopsy", "green","blue"))
})

test_that("clinical ranksum works", {
  expect_silent(clinicalRanksum("blood"))
  expect_silent(clinicalRanksum("biopsy"))
})

# Old but compatible versions of ggplo2 will produce a warning. We well accept a
# warning message for now.
test_that("get tom nodes works", {
  expect_message(getTOMGraphNodes("blood"))
  expect_silent(getTOMGraphNodes("biopsy"))
})

test_that("get tom edges works", {
  expect_silent(getTOMGraphEdges("blood"))
  expect_silent(getTOMGraphEdges("biopsy"))
})


