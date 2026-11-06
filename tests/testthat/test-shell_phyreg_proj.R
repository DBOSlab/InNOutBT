project_files <- function(syst) {
  dir <- tempfile("proj_")
  dir.create(dir)
  shell.phyreg.proj(dir_create_proj = dir, syst = syst)
  file.path(dir, "BayesTraits_shell_phyreg_project")
}

test_that("shell.phyreg.proj creates a Unix project", {
  proj <- project_files("unix")

  expect_setequal(list.files(proj),
                  c("1_producing_BayesTraits_phyreg_inputs.R",
                    "2_producing_BayesTraits_phyreg_shell.R",
                    "3_processing_BayesTraits_phyreg_outputs.R",
                    "BayesTraits_full_phyreg_shell_unix.sh",
                    "Data"))
  expect_setequal(list.files(file.path(proj, "Data")),
                  c("vataireoids_1610_25May2022_BayesTraits.csv",
                    "vatcomb30_24May22_Yule.con.tree"))
  script <- readLines(file.path(proj, "2_producing_BayesTraits_phyreg_shell.R"))
  expect_true(any(grepl('BayesTraits_exe = "BayesTraitsV5"', script, fixed = TRUE)))
  # The shell script expects the run script of its own operating system
  expect_true(any(grepl('syst = "unix"', script, fixed = TRUE)))
})

test_that("shell.phyreg.proj creates a Windows project", {
  proj <- project_files("windows")

  expect_setequal(list.files(proj),
                  c("1_producing_BayesTraits_phyreg_inputs.R",
                    "2_producing_BayesTraits_phyreg_shell.R",
                    "3_processing_BayesTraits_phyreg_outputs.R",
                    "BayesTraits_full_phyreg_shell_windows.ps1",
                    "Data"))
  script <- readLines(file.path(proj, "2_producing_BayesTraits_phyreg_shell.R"))
  expect_true(any(grepl('syst = "windows"', script, fixed = TRUE)))
})

test_that("shell.phyreg.proj detects the operating system when syst = NULL", {
  proj <- project_files(NULL)
  shell <- if (.Platform$OS.type == "unix") "BayesTraits_full_phyreg_shell_unix.sh" else "BayesTraits_full_phyreg_shell_windows.ps1"

  expect_true(file.exists(file.path(proj, shell)))
})

test_that("the example R scripts of the project can be parsed", {
  for (syst in c("unix", "windows")) {
    proj <- project_files(syst)
    for (f in list.files(proj, pattern = "[.]R$", full.names = TRUE)) {
      expect_no_error(parse(f))
    }
  }
})

test_that("the project templates use the current version of BayesTraits", {
  for (syst in c("unix", "windows")) {
    proj <- project_files(syst)
    files <- list.files(proj, pattern = "[.](R|sh|ps1)$", full.names = TRUE)
    text <- unlist(lapply(files, readLines, warn = FALSE))
    expect_false(any(grepl("BayesTraitsV4", text, fixed = TRUE)))
    expect_true(any(grepl('BayesTraits_exe = "BayesTraitsV5"', text, fixed = TRUE)))
  }
})

test_that("the main script opens the folder where the results are saved", {
  for (syst in c("unix", "windows")) {
    proj <- project_files(syst)
    outputs <- readLines(file.path(proj, "3_processing_BayesTraits_phyreg_outputs.R"))
    dir_create <- sub('.*dir_create = "([^"]+)".*', "\\1", outputs[grepl("dir_create", outputs)])
    shell <- list.files(proj, pattern = "[.](sh|ps1)$", full.names = TRUE)
    main <- readLines(shell, warn = FALSE)
    opened <- main[grepl("results_BayesTraits", main, fixed = TRUE)]

    expect_gt(length(opened), 0)
    expect_true(all(grepl(paste0(dir_create, "/"), opened, fixed = TRUE)))
  }
})
