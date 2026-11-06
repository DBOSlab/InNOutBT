utils::data(vatsTree)
utils::data(vatsData)

todaydate <- format(Sys.time(), "%d%b%Y")
folder_name <- paste0("results_BayesTraits_phyreg_input", "/", todaydate)

test_that("phyreg.inputs works for all arguments", {
  phyreg.inputs(tree = vatsTree,
                data = vatsData,
                tipscol = "terminal",
                NodeCount = TRUE,
                logtransf = c("DBH", "NodeCount"),
                traitcols = c("bio12", "bio15"),
                addtraits = c("log10NodeCount", "log10DBH"),
                ordtraits = c("log10DBH", "bio12", "bio15", "log10NodeCount"),
                dir_create = "results_BayesTraits_phyreg_input",
                fileDistData = "BayesTraits_linked_data_bio12_bio15_nnodes.txt",
                fileMeanData = "BayesTraits_mean_data_bio12_bio15_nnodes.txt",
                fileOrigData = "vataireoids_1610_25May2022_BayesTraits_netnodes_logtransf.csv")

  vatsData_new <- utils::read.csv("results_BayesTraits_phyreg_input/vataireoids_1610_25May2022_BayesTraits_netnodes_logtransf.csv")

  expect_equal(class(vatsData_new), "data.frame")
  expect_equal(names(vatsData_new)[16], "log10DBH")
  expect_equal(names(vatsData_new)[23], "log10NodeCount")
  expect_equal(is.numeric(vatsData_new[,16]), TRUE)
  expect_equal(is.numeric(vatsData_new[,23]), TRUE)
  expect_equal(ncol(vatsData_new), 83)
  expect_equal(nrow(vatsData_new), 1610)
})

test_that("phyreg.inputs works without counting net nodes, NodeCount = FALSE", {
  phyreg.inputs(tree = vatsTree,
                data = vatsData,
                tipscol = "terminal",
                NodeCount = FALSE,
                logtransf = c("DBH"),
                traitcols = c("bio12", "bio15"),
                addtraits = "log10DBH",
                ordtraits = c("log10DBH", "bio12", "bio15"),
                dir_create = "results_BayesTraits_phyreg_input",
                fileDistData = "BayesTraits_linked_data_bio12_bio15.txt",
                fileMeanData = "BayesTraits_mean_data_bio12_bio15.txt",
                fileOrigData = "vataireoids_1610_25May2022_BayesTraits_logtransf.csv")

  vatsData_new <- utils::read.csv("results_BayesTraits_phyreg_input/vataireoids_1610_25May2022_BayesTraits_logtransf.csv")

  expect_equal(class(vatsData_new), "data.frame")
  expect_equal(names(vatsData_new)[16], "log10DBH")
  expect_equal(names(vatsData_new)[23], "latitude")
  expect_equal(is.numeric(vatsData_new[,16]), TRUE)
  expect_equal(is.numeric(vatsData_new[,23]), TRUE)
  expect_equal(ncol(vatsData_new), 81)
  expect_equal(nrow(vatsData_new), 1610)
})

test_that("phyreg.inputs works without counting net nodes and log-transform any column, NodeCount = FALSE and logtransf = NULL", {
  phyreg.inputs(tree = vatsTree,
                data = vatsData,
                tipscol = "terminal",
                NodeCount = FALSE,
                logtransf = NULL,
                traitcols = c("bio12", "bio15"),
                addtraits = NULL,
                ordtraits = c("bio12", "bio15"),
                dir_create = "results_BayesTraits_phyreg_input",
                fileDistData = "BayesTraits_linked_data_bio12_bio15.txt",
                fileMeanData = "BayesTraits_mean_data_bio12_bio15.txt",
                fileOrigData = "vataireoids_1610_25May2022_BayesTraits.csv")

  vatsData_new <- utils::read.csv("results_BayesTraits_phyreg_input/vataireoids_1610_25May2022_BayesTraits.csv")

  expect_equal(class(vatsData_new), "data.frame")
  expect_equal(names(vatsData_new)[16], "Sweetia")
  expect_equal(names(vatsData_new)[23], "elev")
  expect_equal(is.numeric(vatsData_new[,16]), FALSE)
  expect_equal(is.numeric(vatsData_new[,23]), TRUE)
  expect_equal(ncol(vatsData_new), 80)
  expect_equal(nrow(vatsData_new), 1610)
})

unlink("results_BayesTraits_phyreg_input", recursive = TRUE)

test_that("phyreg.inputs output errors for wrong inputs", {
  expect_error(phyreg.inputs(tree = vatsTree,
                             data = vatsData,
                             tipscol = "bio12",
                             NodeCount = FALSE,
                             logtransf = NULL,
                             traitcols = c("bio12", "bio15"),
                             addtraits = NULL,
                             ordtraits = c("bio12", "bio15"),
                             dir_create = "results_BayesTraits_phyreg_input",
                             fileDistData = "BayesTraits_linked_data_bio12_bio15.txt",
                             fileMeanData = "BayesTraits_mean_data_bio12_bio15.txt",
                             fileOrigData = "vataireoids_1610_25May2022_BayesTraits.csv"))
})

# A small example where the expected values are easy to check by hand
small_tree <- ape::read.tree(text = "((Aaa:1,Bbb:1):1,Ccc:2);")
small_data <- data.frame(terminal = c("Aaa", "Aaa", "Bbb", "Ccc"),
                         trait1 = c(10, 100, 1000, 10),
                         trait2 = c(100, 300, 500, 700),
                         longitude = 1:4)

read_outputs <- function(folder, mean_file = "BayesTraits_mean_input.txt",
                         linked_file = "BayesTraits_linked_input.txt") {
  list(mean = utils::read.delim(file.path(folder, mean_file), header = FALSE,
                                strip.white = TRUE),
       linked = readLines(file.path(folder, linked_file)))
}

test_that("phyreg.inputs calculates means and links the samples", {
  res <- in_tempdir({
    phyreg.inputs(tree = small_tree, data = small_data, tipscol = "terminal",
                  traitcols = c("trait1", "trait2"))
    read_outputs(file.path("results_BayesTraits_phyreg_input", todaydate()))
  })

  expect_equal(res$mean$V1, c("Aaa", "Bbb", "Ccc"))
  expect_equal(res$mean$V2, c(55, 1000, 10))
  expect_equal(res$mean$V3, c(200, 500, 700))
  expect_equal(res$linked, c("Aaa linked 10,100 100,300", "Bbb linked 1000 500", "Ccc linked 10 700"))
})

test_that("phyreg.inputs counts net nodes and transforms the data", {
  res <- in_tempdir({
    phyreg.inputs(tree = small_tree, data = small_data, tipscol = "terminal",
                  NodeCount = TRUE,
                  logtransf = "trait1",
                  sqrtransf = "trait1",
                  traitcols = "trait2",
                  addtraits = c("log10trait1", "sqrtrait1"),
                  ordtraits = c("log10trait1", "sqrtrait1", "trait2"),
                  fileOrigData = "orig.csv")
    list(out = read_outputs(file.path("results_BayesTraits_phyreg_input", todaydate())),
         orig = utils::read.csv(file.path("results_BayesTraits_phyreg_input", "orig.csv")))
  })

  # The order of the columns follows ordtraits
  expect_equal(res$out$mean$V2, c(1.5, 3, 1))
  expect_equal(res$out$mean$V3, c(2.5, 9, 1))
  expect_equal(res$out$mean$V4, c(200, 500, 700))

  expect_equal(res$orig$NodeCount, c(1, 1, 1, 0))
  expect_equal(res$orig$log10trait1, c(1, 2, 3, 1))
  expect_equal(res$orig$sqrtrait1, c(1, 4, 9, 1))
})

test_that("phyreg.inputs accepts both phylo and treedata trees", {
  phylo_res <- in_tempdir({
    phyreg.inputs(tree = vatsTree@phylo, data = vatsData, tipscol = "terminal",
                  traitcols = "bio12")
    read_outputs(file.path("results_BayesTraits_phyreg_input", todaydate()))
  })
  treedata_res <- in_tempdir({
    phyreg.inputs(tree = vatsTree, data = vatsData, tipscol = "terminal",
                  traitcols = "bio12")
    read_outputs(file.path("results_BayesTraits_phyreg_input", todaydate()))
  })

  expect_equal(phylo_res, treedata_res)
})

test_that("phyreg.inputs stops and reports NA values", {
  bad <- small_data
  bad$trait1[3] <- NA

  expect_error(in_tempdir(phyreg.inputs(tree = small_tree, data = bad, tipscol = "terminal",
                                        traitcols = c("trait1", "trait2"))),
               "trait trait1: rows 3")
})

test_that("phyreg.inputs stops and reports names missing from the tree", {
  bad <- small_data
  bad$terminal[4] <- "Z"

  expect_error(in_tempdir(phyreg.inputs(tree = small_tree, data = bad, tipscol = "terminal",
                                        traitcols = "trait1")),
               "does not match with the tiplabels")
})
