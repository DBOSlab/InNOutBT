# Run code in a new temporary working directory, without printing its messages,
# so that tests leave no files behind in the package.
in_tempdir <- function(code) {
  dir <- tempfile("InNOutBT_")
  dir.create(dir)
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  utils::capture.output(res <- force(code))
  res
}

todaydate <- function() format(Sys.time(), "%d%b%Y")

# The fixtures are real BayesTraits V4 Log and Stones files, reduced to their full
# header plus six MCMC samples (Log) or their last lines (Stones).
# The root is found here, because tests later change the working directory.
fixture_root <- normalizePath(testthat::test_path("fixtures"))
fixture_dir <- function(set) file.path(fixture_root, set)

# Copy only some models of a fixture set, e.g. pattern = "bio12[.]txt"
fixture_subset <- function(set, pattern) {
  dir <- tempfile("logs_")
  dir.create(dir)
  from <- list.files(fixture_dir(set), pattern = pattern, full.names = TRUE)
  file.copy(from, dir)
  dir
}

# Read the MCMC samples by locating their header, independently of the way
# phyreg.outputs counts the lines to skip.
read_samples <- function(set, stem) {
  path <- file.path(fixture_dir(set), paste0("BayesTraits_mean_data_", stem))
  lines <- readLines(path)
  skip <- grep("^Iteration\t", lines) - 1
  utils::read.delim(path, skip = skip, sep = "\t", header = TRUE)
}

read_ll <- function(set, stem) {
  path <- file.path(fixture_dir(set), paste0("BayesTraits_mean_data_", stem))
  lines <- readLines(path)
  as.numeric(sub(".*\t", "", lines[grepl("^Log marginal likelihood:", lines)]))
}

# A small set of empty BayesTraits input files
make_shell_inputs <- function(vars = "bio12") {
  root <- tempfile("phyreg_shell_")
  dirs <- list(root = root, bt = file.path(root, "BT"), inputs = file.path(root, "inputs"))
  dir.create(dirs$bt, recursive = TRUE)
  dir.create(dirs$inputs)
  file.create(file.path(dirs$inputs, c(paste0("BayesTraits_mean_data_", vars, ".txt"),
                                       paste0("BayesTraits_linked_data_", vars, ".txt"),
                                       "example.con.tree")))
  dirs
}

# Run phyreg.shell on those files and return its folders
run_shell <- function(vars = "bio12", syst = "unix", cc_DataTree = FALSE, ...) {
  dirs <- make_shell_inputs(vars)
  phyreg.shell(meanfile_dir = dirs$inputs,
               linkfile_dir = dirs$inputs,
               treefile_dir = dirs$inputs,
               BayesTraits_dir = dirs$bt,
               responvar = "DBH",
               syst = syst,
               dir_create = file.path(dirs$root, "results"),
               cc_DataTree = cc_DataTree,
               ...)
  dirs
}

shell_lines <- function(dirs, syst = "unix") {
  ext <- if (syst == "unix") ".sh" else ".ps1"
  readLines(file.path(dirs$bt, paste0("run_BayesTraits_phyreg_shell_", syst, ext)))
}

param_lines <- function(dirs, name) {
  readLines(file.path(dirs$bt, paste0("run_DBH_", name, ".txt")), warn = FALSE)
}
