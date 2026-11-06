all_transf <- c("Kappa", "Delta", "Lambda", "OU", "UNI", "VR", "Fabric")

run_models <- function(...) {
  phyreg.outputs(logst_dir = fixture_dir("logs_models"),
                 responvar = "DBH",
                 explanvar = c("bio12", "bio15"),
                 dir_create = "out",
                 ...)
}

model_row <- function(res, model) res[res$Model == model, ]

test_that("phyreg.outputs builds one row per model, ordered by likelihood", {
  res <- in_tempdir(run_models(treetransf = all_transf, outformat = "CSV"))

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 10)
  expect_equal(names(res), c("Model", "bio12", "bio15", "Kappa", "Delta", "Lambda",
                             "OU", "Rsquared", "LL", "meanDBH"))
  expect_false(is.unsorted(rev(as.numeric(res$LL))))
  expect_setequal(res$Model,
                  c("DBH ~ bio12 + λ", "DBH ~ bio12", "DBH ~ bio12 + VR",
                    "DBH ~ bio12 + bio15", "DBH ~ bio12 + bio15 + VR",
                    "DBH ~ bio15 + δ", "DBH ~ bio15 + Fabric", "DBH ~ bio15 + κ",
                    "DBH ~ bio15 + λ", "DBH ~ bio15 + OU"))
})

test_that("phyreg.outputs reports the values found in the Log and Stones files", {
  res <- in_tempdir(run_models(treetransf = all_transf, outformat = "CSV"))

  # The header is skipped correctly for each type of model
  for (m in list(c("DBH ~ bio12 + λ", "bio12.txt.Log.Lambda.txt", "bio12.txt.Stones.Lambda.txt"),
                 c("DBH ~ bio12", "bio12.txt.Log.UNI.txt", "bio12.txt.Stones.UNI.txt"),
                 c("DBH ~ bio12 + VR", "bio12.txt.Log.VR.txt", "bio12.txt.Stones.VR.txt"),
                 c("DBH ~ bio15 + Fabric", "bio15.txt.Log.Fabric.txt", "bio15.txt.Stones.Fabric.txt"))) {
    s <- read_samples("logs_models", m[2])
    row <- model_row(res, m[1])
    expect_equal(as.numeric(row$Rsquared), round(mean(s$R.2), 5))
    expect_equal(as.numeric(row$meanDBH), round(10^mean(s$Alpha), 2))
    expect_equal(as.numeric(row$LL), read_ll("logs_models", m[3]))
  }
})

test_that("phyreg.outputs reports each beta with its pMCMC", {
  res <- in_tempdir(run_models(treetransf = all_transf, outformat = "CSV"))

  s <- read_samples("logs_models", "bio12.txt.Log.UNI.txt")
  beta <- round(mean(s$Beta.1), 5)
  p <- if (beta >= 0) mean(s$Beta.1 <= 0) else mean(s$Beta.1 >= 0)
  expect_equal(model_row(res, "DBH ~ bio12")$bio12, paste0(beta, " (", p, ")"))
  expect_true(is.na(model_row(res, "DBH ~ bio12")$bio15))
})

test_that("phyreg.outputs reports estimated and default tree transformations", {
  res <- in_tempdir(run_models(treetransf = all_transf, outformat = "CSV"))

  # Estimated parameter comes from the samples, the others keep their defaults
  for (m in list(c("DBH ~ bio15 + κ", "Kappa"), c("DBH ~ bio15 + δ", "Delta"),
                 c("DBH ~ bio15 + λ", "Lambda"), c("DBH ~ bio15 + OU", "OU"))) {
    s <- read_samples("logs_models", paste0("bio15.txt.Log.", m[2], ".txt"))
    row <- model_row(res, m[1])
    expect_equal(as.numeric(row[[m[2]]]), round(mean(s[[m[2]]]), 5))
    others <- setdiff(c("Kappa", "Delta", "Lambda"), m[2])
    expect_equal(unlist(row[others], use.names = FALSE), rep("1.0", length(others)))
    if (m[2] != "OU") expect_equal(row$OU, "0")
  }

  # Models without those transformations get only the default values
  for (m in c("DBH ~ bio12", "DBH ~ bio12 + VR", "DBH ~ bio15 + Fabric")) {
    row <- model_row(res, m)
    expect_equal(unlist(row[c("Kappa", "Delta", "Lambda", "OU")], use.names = FALSE),
                 c("1.0", "1.0", "1.0", "0"))
  }
})

test_that("phyreg.outputs keeps fixed values of a tree transformation", {
  res <- in_tempdir(phyreg.outputs(logst_dir = fixture_dir("logs_values"),
                                   responvar = "DBH",
                                   explanvar = c("bio12", "bio15"),
                                   treetransf = c("Lambda", "VR"),
                                   outformat = "CSV",
                                   dir_create = "out"))

  expect_equal(nrow(res), 7)
  lambda <- setNames(res$Lambda, res$Model)
  expect_equal(unname(lambda["DBH ~ bio12 + bio15 + λ0.001"]), "0.001")
  expect_equal(unname(lambda["DBH ~ bio12 + bio15 + λ1.0"]), "1.0")
  expect_equal(unname(lambda["DBH ~ bio12 + bio15 + VRλ0.001"]), "0.001")
  expect_equal(unname(lambda["DBH ~ bio12 + bio15 + VRλ1.0"]), "1.0")
  s <- read_samples("logs_values", "bio12_bio15.txt.Log.LambdaEstimate.txt")
  expect_equal(as.numeric(lambda["DBH ~ bio12 + bio15 + λEstimate"]), round(mean(s$Lambda), 5))
  s <- read_samples("logs_values", "bio12_bio15.txt.Log.VRLambdaEstimate.txt")
  expect_equal(as.numeric(lambda["DBH ~ bio12 + bio15 + VRλEstimate"]), round(mean(s$Lambda), 5))
})

test_that("phyreg.outputs drops uniform models when unirates = FALSE", {
  res <- in_tempdir(run_models(treetransf = all_transf, unirates = FALSE, outformat = "CSV"))

  expect_equal(nrow(res), 8)
  expect_false(any(res$Model %in% c("DBH ~ bio12", "DBH ~ bio12 + bio15")))
})

test_that("phyreg.outputs works with only VR, Fabric and UNI models", {
  dir <- fixture_subset("logs_models", "bio12[.]txt[.](Log|Stones)[.](UNI|VR)")
  res <- in_tempdir(phyreg.outputs(logst_dir = dir,
                                   responvar = "DBH",
                                   explanvar = "bio12",
                                   treetransf = c("UNI", "VR"),
                                   outformat = "CSV",
                                   dir_create = "out"))

  expect_equal(names(res), c("Model", "bio12", "Rsquared", "LL", "meanDBH"))
  expect_setequal(res$Model, c("DBH ~ bio12", "DBH ~ bio12 + VR"))
})

test_that("phyreg.outputs calculates Bayes Factors and writes all its files", {
  dir <- tempfile("out_")
  res <- in_tempdir({
    r <- run_models(treetransf = all_transf, bayesfactor = TRUE, outformat = "CSV",
                    outfile = "table", height = 12, width = 15)
    file.rename("out", dir)
    r
  })
  folder <- file.path(dir, todaydate())

  # Compared with the uniform model: twice the difference in log likelihood
  ll <- function(m) read_ll("logs_models", paste0("bio12.txt.Stones.", m, ".txt"))
  expect_equal(as.numeric(model_row(res, "DBH ~ bio12 + VR")$BF),
               round(2 * abs(ll("VR") - ll("UNI")), 5))
  expect_equal(as.numeric(model_row(res, "DBH ~ bio12 + λ")$BF),
               round(2 * abs(ll("Lambda") - ll("UNI")), 5))

  # Without a uniform model to compare with
  expect_match(model_row(res, "DBH ~ bio15 + κ")$BF, "were not compared!", fixed = TRUE)
  # A uniform model is not compared with itself
  expect_true(is.na(model_row(res, "DBH ~ bio12")$BF))

  expect_setequal(list.files(folder),
                  c("table.csv", "Bayes_Factor_phyreg_heatmap.pdf",
                    "Bayes_Factor_phyreg_heatmap_cluster.pdf",
                    "Bayes_Factor_phyreg_heatmap_matrix.csv"))
  expect_gt(file.size(file.path(folder, "Bayes_Factor_phyreg_heatmap.pdf")), 0)
  expect_gt(file.size(file.path(folder, "Bayes_Factor_phyreg_heatmap_cluster.pdf")), 0)

  csv <- utils::read.csv(file.path(folder, "table.csv"), check.names = FALSE)
  expect_equal(nrow(csv), nrow(res))
  heat <- utils::read.csv(file.path(folder, "Bayes_Factor_phyreg_heatmap_matrix.csv"),
                          row.names = 1, check.names = FALSE)
  expect_equal(dim(heat), c(10, 10))
  expect_true(all(c("bio12+λ", "bio15+κ", "bio15+δ") %in% rownames(heat)))
})

test_that("phyreg.outputs writes a Word table", {
  skip_if_not(rmarkdown::pandoc_available())

  dir <- tempfile("out_")
  in_tempdir({
    run_models(treetransf = all_transf, outformat = "Word", outfile = "table",
               tableleg = "My legend")
    expect_false(file.exists("temp.Rmd"))
    file.rename("out", dir)
  })
  folder <- file.path(dir, todaydate())

  expect_equal(list.files(folder), "table.docx")
  expect_gt(file.size(file.path(folder, "table.docx")), 0)
})

test_that("greekLetters replaces the names of the tree transformations", {
  expect_equal(greekLetters(c("a + Kappa", "a + Delta", "a + Lambda", "a + Lmda", "a + OU")),
               c("a + κ", "a + δ", "a + λ", "a + λ", "a + OU"))
})
