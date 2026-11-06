run_lambda <- function(...) {
  run_shell(treetransf = "Lambda", Lambda = "Estimate", ...)
}

test_that("phyreg.shell uses BayesTraitsV5 as the default executable", {
  unix <- shell_lines(run_lambda())
  expect_true(any(grepl("./BayesTraitsV5 $treefile", unix, fixed = TRUE)))
  expect_false(any(grepl("BayesTraitsV4", unix, fixed = TRUE)))
})

test_that("phyreg.shell writes the executable set in BayesTraits_exe", {
  unix <- shell_lines(run_lambda(BayesTraits_exe = "BayesTraitsV4"))
  expect_true(any(grepl("./BayesTraitsV4 $treefile", unix, fixed = TRUE)))
  expect_false(any(grepl("BayesTraitsV5", unix, fixed = TRUE)))

  windows <- shell_lines(run_lambda(syst = "windows", BayesTraits_exe = "BayesTraitsV6"), "windows")
  expect_true(any(grepl("| .\\BayesTraitsV6 $treefile", windows, fixed = TRUE)))
  expect_false(any(grepl("./BayesTraits", windows, fixed = TRUE)))
})

test_that("phyreg.shell keeps only the file name of BayesTraits_exe", {
  unix <- shell_lines(run_lambda(BayesTraits_exe = "/some/folder/BayesTraitsV5"))
  expect_true(any(grepl("^./BayesTraitsV5 [$]treefile", unix)))
})

test_that("phyreg.shell writes a complete Unix script", {
  dirs <- run_shell(vars = c("bio12", "bio12_bio15"),
                    treetransf = c("UNI", "VR"))
  unix <- trimws(shell_lines(dirs))

  expect_equal(unix[1], "#!/bin/zsh")
  expect_true("treefile=example.con.tree" %in% unix)
  expect_true("meanfile1=BayesTraits_mean_data_bio12.txt" %in% unix)
  expect_true("meanfile2=BayesTraits_mean_data_bio12_bio15.txt" %in% unix)
  expect_true("paramfile1=run_DBH_bio12_UNI.txt" %in% unix)
  expect_true("paramfile4=run_DBH_bio12_bio15_VR.txt" %in% unix)

  # Each model runs with its own data file and renames its own outputs
  expect_true("## 2. Run DBH ~ bio12 + VR" %in% unix)
  expect_true("./BayesTraitsV5 $treefile $meanfile1 < $paramfile2" %in% unix)
  expect_true("./BayesTraitsV5 $treefile $meanfile2 < $paramfile4" %in% unix)
  expect_true("mv BayesTraits_mean_data_bio12.txt.Log.txt BayesTraits_mean_data_bio12.txt.Log.UNI.txt" %in% unix)
  expect_true("mv BayesTraits_mean_data_bio12.txt.Stones.txt BayesTraits_mean_data_bio12.txt.Stones.VR.txt" %in% unix)
  expect_true("mv BayesTraits_mean_data_bio12.txt.Output.trees BayesTraits_mean_data_bio12.txt.Output.VR.trees" %in% unix)
  expect_true("mv BayesTraits_mean_data_bio12.txt.VarRates.txt BayesTraits_mean_data_bio12.txt.VarRates.VR.txt" %in% unix)
  expect_false("mv BayesTraits_mean_data_bio12.txt.Output.trees BayesTraits_mean_data_bio12.txt.Output.UNI.trees" %in% unix)

  # Moving the results and cleaning the folder
  expect_true("mv *Log* BayesTraits_phyreg_outputs_log_stone/" %in% unix)
  expect_true("mv *Schedule* BayesTraits_phyreg_outputs_sch_trees/" %in% unix)
  expect_true("rm $treefile" %in% unix)
  expect_true("rm $meanfile2" %in% unix)
  expect_true("rm $paramfile4" %in% unix)
  expect_true("rm *linked_data*" %in% unix)
  expect_true("rm run_BayesTraits_phyreg_shell_unix.sh" %in% unix)
  expect_false(any(grepl("Remove-Item|Rename-Item|Move-Item", unix)))
})

test_that("phyreg.shell writes a complete Windows script", {
  dirs <- run_shell(vars = c("bio12", "bio15"), syst = "windows",
                    treetransf = c("UNI", "VR"))
  win <- trimws(shell_lines(dirs, "windows"))

  expect_true("$treefile='example.con.tree'" %in% win)
  expect_true("$meanfile2='BayesTraits_mean_data_bio15.txt'" %in% win)
  expect_true("$paramfile1='run_DBH_bio12_UNI.txt'" %in% win)
  expect_true("Get-Content $paramfile4 | .\\BayesTraitsV5 $treefile $meanfile2" %in% win)
  expect_true("Rename-Item -Path BayesTraits_mean_data_bio12.txt.Log.txt -NewName BayesTraits_mean_data_bio12.txt.Log.UNI.txt" %in% win)
  expect_true("Rename-Item -Path BayesTraits_mean_data_bio12.txt.VarRates.txt -NewName BayesTraits_mean_data_bio12.txt.VarRates.VR.txt" %in% win)
  expect_true("Rename-Item -Path BayesTraits_mean_data_bio12.txt.Output.trees -NewName BayesTraits_mean_data_bio12.txt.Output.VR.trees" %in% win)
  expect_true("Move-Item -Path *Log* -Destination BayesTraits_phyreg_outputs_log_stone/" %in% win)
  expect_true("Remove-Item $treefile" %in% win)
  expect_true("Remove-Item $meanfile2" %in% win)
  expect_true("Remove-Item $paramfile4" %in% win)
  expect_true("Remove-Item *linked_data*" %in% win)
  expect_true("Remove-Item run_BayesTraits_phyreg_shell_windows.ps1" %in% win)
  expect_false(any(grepl("^(mv|rm) ", win)))
})

test_that("phyreg.shell saves the same files in the results folder", {
  dirs <- run_shell(treetransf = "Lambda", Lambda = c("0.001", "Estimate"))
  folder <- file.path(dirs$root, "results", todaydate())

  expected <- c("run_BayesTraits_phyreg_shell_unix.sh",
                "run_DBH_bio12_Lambda0.001.txt", "run_DBH_bio12_LambdaEstimate.txt")
  expect_setequal(list.files(folder), expected)
  expect_setequal(list.files(dirs$bt), expected)
  expect_equal(readLines(file.path(folder, expected[2]), warn = FALSE),
               param_lines(dirs, "bio12_Lambda0.001"))
})

test_that("phyreg.shell writes the MCMC settings in each run file", {
  dirs <- run_shell(treetransf = "UNI", bi = "10", it = "20", sa = "5", st = c("2", "3"))

  expect_equal(param_lines(dirs, "bio12_UNI"),
               c("9", "2", "bi 10", "it 20", "sa 5", "stones 2 3",
                 "DistData BayesTraits_linked_data_bio12.txt", "Run"))
})

test_that("phyreg.shell has default MCMC settings", {
  dirs <- run_shell(treetransf = "UNI")

  expect_equal(param_lines(dirs, "bio12_UNI")[3:6],
               c("bi 30000000", "it 130000000", "sa 100000", "stones 500 50000"))
})

test_that("phyreg.shell writes the command of each model without fixed values", {
  dirs <- run_shell(treetransf = c("UNI", "VR", "Fabric", "Kappa", "Delta", "Lambda", "OU"))

  expect_false("vr" %in% param_lines(dirs, "bio12_UNI"))
  expect_equal(param_lines(dirs, "bio12_VR")[3], "vr")
  expect_equal(param_lines(dirs, "bio12_Fabric")[3], "Fabric")
  for (m in c("Kappa", "Delta", "Lambda", "OU")) {
    expect_equal(param_lines(dirs, paste0("bio12_", m))[3], m)
  }
})

test_that("phyreg.shell writes estimated and fixed values", {
  dirs <- run_shell(treetransf = c("UNI", "VR", "Lambda", "Kappa"),
                    Lambda = c("0.001", "Estimate"),
                    Kappa = "1.0")

  expect_equal(param_lines(dirs, "bio12_Lambda0.001")[3], "Lambda 0.001")
  expect_equal(param_lines(dirs, "bio12_LambdaEstimate")[3], "Lambda")
  expect_equal(param_lines(dirs, "bio12_Kappa1.0")[3], "Kappa 1.0")
  expect_equal(param_lines(dirs, "bio12_VR")[3], "vr")
  expect_equal(param_lines(dirs, "bio12_UNI")[3], "bi 30000000")
})

test_that("phyreg.shell combines VR with Lambda when VRLambda = TRUE", {
  dirs <- run_shell(treetransf = c("VR", "Lambda", "Kappa"),
                    Lambda = c("0.001", "Estimate"),
                    Kappa = c("0.5", "Estimate"),
                    VRLambda = TRUE)

  expect_equal(param_lines(dirs, "bio12_VRLambda0.001")[3:4], c("vr", "Lambda 0.001"))
  expect_equal(param_lines(dirs, "bio12_VRLambdaEstimate")[3:4], c("vr", "Lambda"))
  expect_equal(param_lines(dirs, "bio12_Lambda0.001")[3], "Lambda 0.001")
  expect_equal(param_lines(dirs, "bio12_LambdaEstimate")[3], "Lambda")
  expect_equal(param_lines(dirs, "bio12_Kappa0.5")[3], "Kappa 0.5")
  expect_equal(param_lines(dirs, "bio12_KappaEstimate")[3], "Kappa")
  expect_false(file.exists(file.path(dirs$bt, "run_DBH_bio12_VR.txt")))

  unix <- trimws(shell_lines(dirs))
  expect_true("mv BayesTraits_mean_data_bio12.txt.VarRates.txt BayesTraits_mean_data_bio12.txt.VarRates.VRLambda0.001.txt" %in% unix)
})

test_that("phyreg.shell copies the data and tree only when cc_DataTree = TRUE", {
  dirs <- run_lambda(cc_DataTree = TRUE)
  expect_true(all(c("BayesTraits_mean_data_bio12.txt", "BayesTraits_linked_data_bio12.txt",
                    "example.con.tree") %in% list.files(dirs$bt)))

  dirs <- run_lambda(cc_DataTree = FALSE)
  expect_false(any(c("BayesTraits_mean_data_bio12.txt", "example.con.tree") %in% list.files(dirs$bt)))
})

test_that("phyreg.shell detects the operating system when syst = NULL", {
  dirs <- make_shell_inputs()
  phyreg.shell(meanfile_dir = dirs$inputs, linkfile_dir = dirs$inputs,
               treefile_dir = dirs$inputs, BayesTraits_dir = dirs$bt,
               responvar = "DBH", treetransf = "UNI",
               dir_create = file.path(dirs$root, "results"), cc_DataTree = FALSE)

  ext <- if (.Platform$OS.type == "unix") "unix.sh" else "windows.ps1"
  expect_true(file.exists(file.path(dirs$bt, paste0("run_BayesTraits_phyreg_shell_", ext))))
})
