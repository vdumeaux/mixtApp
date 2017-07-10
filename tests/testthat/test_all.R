library(testthat)

# First get all available tissues in the dataset
tissues = names(bresat)

# Note: We do not check for correctness, just that all the functions run without
# crashing. The rational behind this is that the data in the package is
# continuously changing which would make it necessary to update the tests
# whenever we updated the data.

test_that("heatmap works", {
  for (tissue in tissues) {
  	module = names(bresat[[tissue]])[1]
  	expect_silent(cohort_heatmap(tissue, module))
  	expect_silent(cohort_heatmap(tissue,module,"erp"))
  }
})


test_that("scatterplot works", {
  x_tissue = tissues[1]
  y_tissue = tissues[2]
  x_module = names(bresat[[x_tissue]])[1]
  y_module =  names(bresat[[y_tissue]])[1]
  expect_silent(cohort_scatterplot(x_tissue,x_module, y_tissue, y_module))
})

test_that("boxplot works", {
  tissue = tissues[1]
  module = names(bresat[[tissue]])[1]
  orderByTissue= tissue
  orderByModule=module

  expect_silent(cohort_boxplot(tissue, module, orderByTissue, orderByModule))
  orderByTissue= tissues[2]
  expect_silent(cohort_boxplot(tissue, module, orderByTissue, orderByModule))
})

test_that("get modules works", {
	for(tissue in tissues) {
		  expect_silent(getModules(tissue))
	}
})

test_that("get all genes works", {
  expect_silent(getAllGenes())
})

test_that("get all modules works", {
	genes = getAllGenes()
	gene = genes[1]
  expect_silent(getAllModules(gene))
})

# test_that("get all genes and modules works", {
#   expect_silent(getAllGenesAndModules())
# })

test_that("get all tissues works", {
  expect_equal(getAllTissues(), tissues)
})

test_that("get gene list works", {
	tissue = tissues[1]
	module = names(bresat[[tissue]])[1]
  expect_silent(getGeneList(tissue,module))
})

test_that("get en richment scores works", {
	tissue = tissues[1]
	module = names(bresat[[tissue]])[1]
  expect_silent(getEnrichmentScores(tissue, module))

	tissue = tissues[2]
	module = names(bresat[[tissue]])[1]
  expect_silent(getEnrichmentScores(tissue, module))
})

test_that("get gene set names works", {
  expect_silent(getGeneSetNames())
})

test_that("get nerichment for tissue works", {
	for(tissue in tissues){
  expect_silent(getEnrichmentForTissue(tissue))
	}
})

# only check first two tissues
test_that("get go terms works", {
	for(tissue in tissues[c(1,2)]){
		module = names(bresat[[tissue]])[1]
  	expect_silent(getGOTerms(tissue, module))
	}
})

test_that("get common genes works", {
	tissue = tissues[1]
	module = names(bresat[[tissue]])[1]
  geneset ="REACTOME_RNA_POL_I_TRANSCRIPTION"
  expect_silent(getCommonGenes(tissue, module, geneset))
})

test_that("get common go terms works", {
	tissue = tissues[1]
	module = names(bresat[[tissue]])[1]
  expect_silent(getCommonGOTermGenes(tissue, module, "GO:0070848"))
})

test_that("user enrichment scores works", {
	tissue = tissues[1]
	genes = getAllGenes()[c(1,2)]
  expect_silent(userEnrichmentScores(tissue, genes))
})

test_that("common enrichment score genes works", {
	tissue = tissues[1]
	module = names(bresat[[tissue]])[1]
	genes = getAllGenes()[c(1,2)]
  expect_silent(commonEnrichmentScoreGenes(tissue, module, genes))
})

test_that("get go term names works ", {
  expect_silent(getGOTermNames())
})

test_that("get go scores for tissue work", {
	tissue = tissues[1]
  expect_silent(getGOScoresForTissue(tissue, "B cell receptor signaling pathway"))
})

test_that("gene overlap test works", {
  x_tissue = tissues[1]
  y_tissue = tissues[2]
  expect_silent(geneOverlapTest(x_tissue, y_tissue))
})


test_that("patient ranksum works", {
  x_tissue = tissues[1]
  y_tissue = tissues[2]
  cohort = "all"
  expect_silent(patientRankSum(x_tissue, y_tissue, cohort))
  expect_silent(patientRankSum(y_tissue, x_tissue, cohort))
})

test_that("comparison analyses work", {
	x_tissue = tissues[1]
  y_tissue = tissues[2]

  x_module = names(bresat[[x_tissue]])[1]
  y_module =  names(bresat[[y_tissue]])[1]
  expect_silent(comparisonAnalyses(x_tissue, y_tissue, x_module, y_module))
})

test_that("clinical ranksum works", {
	for(tissue in tissues[c(1,2)]){
  	expect_silent(clinicalRanksum(tissue))
	}
})

# Old but compatible versions of ggplo2 will produce a warning. We well accept a
# warning message for now.
test_that("get tom nodes works", {
	tissue = tissues[1]
	expect_message(getTOMGraphNodes(tissue))
	tissue = tissues[2]
	expect_silent(getTOMGraphNodes(tissue))
})

test_that("get tom edges works", {
	for(tissue in tissues[c(1,2)]){
	  expect_silent(getTOMGraphEdges(tissue))
	}
})


