#' Create shell scripts to run BayesTraits phylogenetic regression analyses using samples of traits data
#'
#' @author Domingos Cardoso & Matt Lavin
#'
#' @description By reporting the directory paths to the tree and data files,
#' this function produces shell scripts in both Unix and Windows PowerShell
#' formats for running one or multiple phylogenetic regression analyses with
#' Meade & Pagel's (2022) [BayesTraits](http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/BayesTraitsV4.0.0.html)
#' program. In addition to generating the shell script directly in the BayesTraits
#' directory, the function also creates BayesTraits batch files with specified
#' parameters (e.g., bi, it, sa, and st) and model of tree transformations.
#' See [BayesTraits V4.0.0 Manual](http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/Files/BayesTraitsV4.0.0-Manual.pdf)
#' for further details on the available BayesTraits parameters for phylogenetic
#' regression so as to better understand how to set the arguments of the current
#' function \code{BayesTraits.shell}. You might want to consider using
#' the function \code{\link{BayesTraits.inputs}} to generate properly formatted input files
#' of mean and linked samples of trait data prior to using \code{BayesTraits.shell}.
#'
#' @usage
#' BayesTraits.shell(meanfile_dir = NULL,
#'                   linkfile_dir = NULL,
#'                   treefile_dir = NULL,
#'                   BayesTraits_dir = NULL,
#'                   responvar = NULL,
#'                   treetransf = NULL,
#'                   bi = "30000000",
#'                   it = "130000000",
#'                   sa = "100000",
#'                   st = c("500", "50000"),
#'                   syst = NULL,
#'                   dir_create = "results_BayesTraits_shell",
#'                   cc_DataTree = TRUE)
#'
#' @param meanfile_dir Path to the folder directory where the mean data file(s)
#' of samples of trait data is(are) stored.
#'
#' @param linkfile_dir Path to the folder directory where the linked data file(s)
#' of samples of trait data is(are) stored.
#'
#' @param treefile_dir Path to the folder directory where the tree file is stored.
#'
#' @param BayesTraits_dir Path to the folder directory where the BayesTraits
#' executable file is stored.
#'
#' @param responvar Report the name of the response variable of your analysis.
#'
#' @param treetransf In vector format, report tree transformation models to be
#' estimated, such as *Kappa*, *Delta*, *Lambda*, *OU* (Ornstein–Uhlenbeck).
#' A model of variable rates *VR* or the *Fabric* model can be reported here,
#' the latter of which is a statistical model accommodating an uneven evolutionary
#' landscape, as described by [Pagel et al. (2022)](https://doi.org/10.1038/s41467-022-28595-z).
#' We use *UNI* as a designation for uniform rates of evolution (no *VR* estimated)
#' and no other tree transformation models estimated.
#'
#' @param bi Set the number of iterations for the MCMC burnin. Report this number
#' as a character (i.e., in quotes).
#'
#' @param it Set the number of total iterations for the MCMC chain. Report this number
#' as a character (i.e., in quotes).
#'
#' @param sa Set the sample frequency. Report this number as a character (i.e., in quotes).
#'
#' @param st The stepping stone sampler to estimate the marginal likelihood. Report
#' the number of stones to use and the number of iterations for each stone as
#' characters (i.e., a vector of numbers each in quotes).
#'
#' @param syst Report the operating system ("unix" or "windows").
#'
#' @param dir_create Path to the computer's working directory or where the file(s) will
#' be saved. The default setting creates a directory "results_BayesTraits_output" in
#' which the results will be saved within a subfolder named by the current date.
#'
#' @param cc_DataTree Logical, the default is \code{TRUE}. A copy of the tree file and
#' the mean and linked data files will be made in the BayesTraits folder where the
#' BayesTraits executable file is located.
#'
#' @seealso \code{\link{BayesTraits.inputs}}
#' @seealso \code{\link{BayesTraits.outputs}}
#'
#' @export
#'

BayesTraits.shell <- function (meanfile_dir = NULL,
                               linkfile_dir = NULL,
                               treefile_dir = NULL,
                               BayesTraits_dir = NULL,
                               responvar = NULL,
                               treetransf = NULL,
                               bi = "30000000",
                               it = "130000000",
                               sa = "100000",
                               st = c("500", "50000"),
                               syst = NULL,
                               dir_create = "results_BayesTraits_shell",
                               cc_DataTree = TRUE) {

  #_____________________________________________________________________________
  # Create a new directory with current date to save the shell script and run command
  # files at the R directory.
  # If there is no directory... make one
  todaydate <- gsub(" ", "", format(Sys.time(), "%d %b %Y"))
  if (!dir.exists(paste0(dir_create, "/"))) {
    dir.create(paste0(dir_create, "/"))
  }
  if (!dir.exists(paste0(dir_create, "/", todaydate))) {
    dir.create(paste0(dir_create, "/", todaydate))
  }
  # If directory was created during a previous search, get its name to save results
  folder_name <- paste0(paste0(dir_create, "/"), todaydate)


  if (is.null(syst)) {
    syst <- .Platform$OS.type
  }

  #_____________________________________________________________________________
  # Start creating the shell script
  if (syst == "unix"){
    title = "#!/bin/zsh"
    shellfile_bt = paste0(BayesTraits_dir, "/run_BayesTraits_shell_", syst, ".sh")
    if (!is.null(dir_create)) {
      shellfile_wd = paste0(folder_name, "/run_BayesTraits_shell_", syst, ".sh")
    }
  }
  if (syst == "windows"){
    title = "#!/bin/bash"
    shellfile_bt = paste0(BayesTraits_dir, "/run_BayesTraits_shell_", syst, ".ps1")
    if (!is.null(dir_create)) {
      shellfile_wd = paste0(folder_name, "/run_BayesTraits_shell_", syst, ".ps1")
    }
  }

  # Get the file names of input data
  tf <- grepl("mean_data", list.files(meanfile_dir))
  meanfile_name <- list.files(meanfile_dir)[tf]
  tf <- grepl("linked_data", list.files(linkfile_dir))
  linkfile_name <- list.files(linkfile_dir)[tf]
  tf <- grepl("con[.]tree$|[.]tree$|[.]tre$", list.files(treefile_dir))
  treefile_name <- list.files(treefile_dir)[tf]

  cat(title, "\n", "\n", file = shellfile_bt, sep="")
  cat("## input files\n", file = shellfile_bt, append=T)
  if (syst == "unix") {
    cat(paste0("treefile=", treefile_name), "\n", "\n", file = shellfile_bt, append=T)
  }
  if (syst == "windows") {
    cat(paste0("$treefile=\'", treefile_name, "\'"), "\n", "\n", file = shellfile_bt, append=T)
  }
  if (!is.null(dir_create)) {
    cat(title, "\n", "\n", file = shellfile_wd, sep="")
    cat("## input files\n", file = shellfile_wd, append=T)
    if (syst == "unix") {
      cat(paste0("treefile=", treefile_name), "\n", "\n", file = shellfile_wd, append=T)
    }
    if (syst == "windows") {
      cat(paste0("$treefile=\'", treefile_name, "\'"), "\n", "\n", file = shellfile_wd, append=T)
    }
  }


  meanfiles_nbrs <- vector()
  for(i in seq_along(meanfile_name)) {
    meanfiles_nbrs[i] <- paste0("meanfile", i)
    if (syst == "unix") {
      cat(paste0("meanfile", i, "=", meanfile_name[i]), "\n", file = shellfile_bt, append=T)
    }
    if (syst == "windows") {
      cat(paste0("$meanfile", i, "=\'", meanfile_name[i], "\'"), "\n", file = shellfile_bt, append=T)
    }
    if (!is.null(dir_create)) {
      if (syst == "unix") {
        cat(paste0("meanfile", i, "=", meanfile_name[i]), "\n", file = shellfile_wd, append=T)
      }
      if (syst == "windows") {
        cat(paste0("$meanfile", i, "=\'", meanfile_name[i], "\'"), "\n", file = shellfile_wd, append=T)
      }
    }
  }

  cat("\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("\n", file = shellfile_wd, append=T)
  }

  cat("mkdir", "BayesTraits_outputs_log_stone", "\n", file = shellfile_bt, append=T)
  cat("mkdir", "BayesTraits_outputs_sch_trees", "\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("mkdir", "BayesTraits_outputs_log_stone", "\n", file = shellfile_wd, append=T)
    cat("mkdir", "BayesTraits_outputs_sch_trees", "\n", file = shellfile_wd, append=T)
  }

  cat("\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("\n", file = shellfile_wd, append=T)
  }

  param_temp <- list()
  temp <- vector()
  run_temp <- list()
  temp_run <- vector()
  for (l in seq_along(meanfile_name)) {
    explanvar <- gsub("BayesTraits_mean_data_", "", meanfile_name[l])
    explanvar <- gsub("[.]txt", "", explanvar)
    for(i in seq_along(treetransf)) {
      temp[i] <- paste0("paramfile", "=", "run_", responvar, "_", explanvar, "_", treetransf[i], ".txt")

      if (treetransf[i] == "UNI") {
        temp_run[i] <- paste0("9", "\n",
                              "2", "\n",
                              "bi ", bi, "\n",
                              "it ", it, "\n",
                              "sa ", sa, "\n",
                              "stones ", st[1], " ", st[2], "\n",
                              "DistData ",
                              linkfile_name[l], "\n",
                              "Run")
      }
      if (treetransf[i] == "Kappa") {
        temp_run[i] <- paste0("9", "\n",
                              "2", "\n",
                              "Kappa", "\n",
                              "bi ", bi, "\n",
                              "it ", it, "\n",
                              "sa ", sa, "\n",
                              "stones ", st[1], " ", st[2], "\n",
                              "DistData ",
                              linkfile_name[l], "\n",
                              "Run")
      }
      if (treetransf[i] == "Delta") {
        temp_run[i] <- paste0("9", "\n",
                              "2", "\n",
                              "Delta", "\n",
                              "bi ", bi, "\n",
                              "it ", it, "\n",
                              "sa ", sa, "\n",
                              "stones ", st[1], " ", st[2], "\n",
                              "DistData ",
                              linkfile_name[l], "\n",
                              "Run")
      }
      if (treetransf[i] == "Lambda") {
        temp_run[i] <- paste0("9", "\n",
                              "2", "\n",
                              "Lambda", "\n",
                              "bi ", bi, "\n",
                              "it ", it, "\n",
                              "sa ", sa, "\n",
                              "stones ", st[1], " ", st[2], "\n",
                              "DistData ",
                              linkfile_name[l], "\n",
                              "Run")
      }
      if (treetransf[i] == "OU") {
        temp_run[i] <- paste0("9", "\n",
                              "2", "\n",
                              "OU", "\n",
                              "bi ", bi, "\n",
                              "it ", it, "\n",
                              "sa ", sa, "\n",
                              "stones ", st[1], " ", st[2], "\n",
                              "DistData ",
                              linkfile_name[l], "\n",
                              "Run")
      }
      if (treetransf[i] == "VR") {
        temp_run[i] <- paste0("9", "\n",
                              "2", "\n",
                              "vr", "\n",
                              "bi ", bi, "\n",
                              "it ", it, "\n",
                              "sa ", sa, "\n",
                              "stones ", st[1], " ", st[2], "\n",
                              "DistData ",
                              linkfile_name[l], "\n",
                              "Run")
      }
      if (treetransf[i] == "Fabric") {
        temp_run[i] <- paste0("9", "\n",
                              "2", "\n",
                              "Fabric", "\n",
                              "bi ", bi, "\n",
                              "it ", it, "\n",
                              "sa ", sa, "\n",
                              "stones ", st[1], " ", st[2], "\n",
                              "DistData ",
                              linkfile_name[l], "\n",
                              "Run")
      }
    }
    param_temp[[l]] <- temp
    run_temp[[l]] <- temp_run
  }

  param_temp <- unlist(param_temp)
  run_temp <- unlist(run_temp)
  for(i in seq_along(param_temp)) {
    if (syst == "unix") {
      temp <- gsub("paramfile", paste0("paramfile", i), param_temp[i])
    }
    if (syst == "windows") {
      temp <- paste0(paste0("$paramfile", i, "="), "\'", gsub("paramfile=", "", param_temp[i]), "\'")
    }
    cat(temp, "\n", file = shellfile_bt, append=T)
    if (!is.null(dir_create)) {
      cat(temp, "\n", file = shellfile_wd, append=T)
    }
    # Writing the individual files with blocks of BayesTraits run commands
    if (!is.null(dir_create)) {
      runfiles = paste0(folder_name, "/", gsub("paramfile[=]", "", param_temp[i]))
      cat(file = runfiles, sep="")
      cat(run_temp[i], file = runfiles, append=T)
    }
    runfiles = paste0(BayesTraits_dir, "/", gsub("paramfile[=]", "", param_temp[i]))
    cat(file = runfiles, sep="")
    cat(run_temp[i], file = runfiles, append=T)
  }

  # Back to writing the script into the shell file
  cat("\n", "\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("\n", "\n", file = shellfile_wd, append=T)
  }

  commands <- vector()
  command_heads <- vector()
  for(i in seq_along(param_temp)) {
    temp <- gsub(paste0("paramfile=run_", responvar, "_|[.]txt"), "", param_temp[i])
    temp <- gsub("_", " + ", temp)
    command_heads[i] <- paste0("## ", i, ".", " Run ", responvar, " ~ ", temp)

    commands[i] <- paste0("./", gsub(".+[/]", "" , "BayesTraitsV4"), " $treefile", " $meanfile", " <", paste0(" $paramfile", i))

    names(commands)[i] <- temp
  }

  temp <- gsub(" [^ ]+$", "", names(commands))
  temp <- gsub("\\s[+]$", "", temp)

  explanvar <- gsub("BayesTraits_mean_data_", "", meanfile_name)
  explanvar <- gsub("[.]txt", "", explanvar)
  explanvar <- gsub("_", " + ", explanvar)

  for (i in seq_along(explanvar)) {
    commands[temp %in% explanvar[i]] <- gsub("meanfile", paste0("meanfile", i), commands[temp %in% explanvar[i]])
  }

  if (syst == "windows") {
    commands <- gsub("[.]/BayesTraitsV4", "", commands)
    temp <- gsub(" .+\\s<", "", commands)
    commands <- paste0("Get-Content", temp, " | ", ".\\BayesTraitsV4", gsub(" <.+", "", commands))
  }

  #_____________________________________________________________________________
  # Getting each block of files names to be changed after each BayesTraits run
  logs <- list()
  stones <- list()
  scheds <- list()
  outrees <- list()
  varates <- list()
  temp_logs <- vector()
  temp_stones <- vector()
  temp_scheds <- vector()

  if (any(treetransf == "VR")) {
    temp_outrees <- vector()
    temp_varates <- vector()
  }

  for(i in seq_along(meanfile_name)) {
    for (l in seq_along(treetransf)) {
      if (syst == "unix") {
        temp_logs[l] <- paste0("mv ", meanfile_name[i], ".Log.txt ", meanfile_name[i], ".Log.", treetransf[l], ".txt")
        temp_stones[l] <- paste0("mv ", meanfile_name[i], ".Stones.txt ", meanfile_name[i], ".Stones.", treetransf[l], ".txt")
        temp_scheds[l] <- paste0("mv ", meanfile_name[i], ".Schedule.txt ", meanfile_name[i], ".Schedule.", treetransf[l], ".txt")
      }
      if (syst == "windows") {
        temp_logs[l] <- paste0("Rename-Item -Path ", meanfile_name[i], ".Log.txt ", "-NewName ", meanfile_name[i], ".Log.", treetransf[l], ".txt")
        temp_stones[l] <- paste0("Rename-Item -Path ", meanfile_name[i], ".Stones.txt ", "-NewName ", meanfile_name[i], ".Stones.", treetransf[l], ".txt")
        temp_scheds[l] <- paste0("Rename-Item -Path ", meanfile_name[i], ".Schedule.txt ", "-NewName ", meanfile_name[i], ".Schedule.", treetransf[l], ".txt")
      }
      if (treetransf[l] == "VR") {
        if (syst == "unix") {
          temp_outrees[l] <- paste0("mv ", meanfile_name[i], ".Output.trees ", meanfile_name[i], ".Output.", treetransf[l], ".trees")
          temp_varates[l] <- paste0("mv ", meanfile_name[i], ".VarRates.txt ", meanfile_name[i], ".VarRates.", treetransf[l], ".txt")
        }
        if (syst == "windows") {
          temp_outrees[l] <- paste0("Rename-Item -Path ", meanfile_name[i], ".Output.trees ", "-NewName ", meanfile_name[i], ".Output.", treetransf[l], ".trees")
          temp_varates[l] <- paste0("Rename-Item -Path ", meanfile_name[i], ".VarRates.txt ", "-NewName ", meanfile_name[i], ".VarRates.", treetransf[l], ".txt")
        }
      }
    }
    logs[[i]] <- temp_logs
    stones[[i]] <- temp_stones
    scheds[[i]] <- temp_scheds
    if (any(treetransf == "VR")) {
      outrees[[i]] <- temp_outrees
      varates[[i]] <- temp_varates
    }
  }

  temp <- append(unlist(logs), unlist(stones))
  temp <- append(temp, unlist(scheds))
  if (any(treetransf == "VR")) {
    temp <- append(temp, unlist(outrees))
    temp <- append(temp, unlist(varates))
  }
  allnames_tochange <- temp[!is.na(temp)]

  if (syst == "unix") {
    vars <- gsub("mv BayesTraits_mean_data_", "", allnames_tochange)
  }
  if (syst == "windows") {
    vars <- gsub("Rename-Item -Path BayesTraits_mean_data_", "", allnames_tochange)
  }
  vars <- gsub("[.]txt.+", "", vars)
  mods <- gsub("[.]txt$|[.]trees", "", allnames_tochange)
  mods <- gsub(".+[.]", "", mods)

  change_names <- list()
  temp <- list()
  for (i in seq_along(unique(vars))) {
    tf <- vars %in% unique(vars)[i]
    for (l in seq_along(treetransf)) {
      temp[[l]] <- allnames_tochange[tf][mods[tf] %in% treetransf[l]]
    }
    change_names[[i]] <- temp
  }

  paste_names <- list()
  for (n in seq(change_names)) {
    for (l in seq(change_names[[n]])) {
      change_names[[n]][[l]] <- paste0(change_names[[n]][[l]], collapse = "\n")
    }
  }

  chunck_change_names <- unlist(change_names[[1]])
  for (i in 2:length(change_names)) {
    chunck_change_names <- append(chunck_change_names, unlist(change_names[[i]]))
  }

  for(i in seq_along(param_temp)) {
    if (!is.null(dir_create)) {
      cat(command_heads[i], "\n", file = shellfile_wd, append=T)
      cat(commands[i], "\n", file = shellfile_wd, append=T)
      cat(chunck_change_names[i], "\n", "\n", "\n", file = shellfile_wd, append=T)
    }
    cat(command_heads[i], "\n", file = shellfile_bt, append=T)
    cat(commands[i], "\n", file = shellfile_bt, append=T)
    cat(chunck_change_names[i], "\n", "\n", "\n", file = shellfile_bt, append=T)
  }

  #_____________________________________________________________________________
  # Moving specific files to specific folders after analysis is complete
  if (syst == "unix") {
    cat("mv", "*Log*", "BayesTraits_outputs_log_stone/", "\n", file = shellfile_bt, append=T)
    cat("mv", "*Stones*", "BayesTraits_outputs_log_stone/", "\n", "\n", file = shellfile_bt, append=T)
    cat("mv", "*Schedule*", "BayesTraits_outputs_sch_trees/", "\n", file = shellfile_bt, append=T)
    cat("mv", "*Output.VR.trees*", "BayesTraits_outputs_sch_trees/", "\n", file = shellfile_bt, append=T)
    cat("mv", "*VarRates*", "BayesTraits_outputs_sch_trees/", file = shellfile_bt, append=T)
  }
  if (syst == "windows") {
    cat("Move-Item -Path", "*Log*", "-Destination", "BayesTraits_outputs_log_stone/", "\n", file = shellfile_bt, append=T)
    cat("Move-Item -Path", "*Stones*", "-Destination", "BayesTraits_outputs_log_stone/", "\n", "\n", file = shellfile_bt, append=T)
    cat("Move-Item -Path", "*Schedule*", "-Destination", "BayesTraits_outputs_sch_trees/", "\n", file = shellfile_bt, append=T)
    cat("Move-Item -Path", "*Output.VR.trees*", "-Destination", "BayesTraits_outputs_sch_trees/", "\n", file = shellfile_bt, append=T)
    cat("Move-Item -Path", "*VarRates*", "-Destination", "BayesTraits_outputs_sch_trees/", file = shellfile_bt, append=T)
  }
  if (!is.null(dir_create)) {
    if (syst == "unix") {
      cat("mv", "*Log*", "BayesTraits_outputs_log_stone/", "\n", file = shellfile_wd, append=T)
      cat("mv", "*Stones*", "BayesTraits_outputs_log_stone/", "\n", "\n", file = shellfile_wd, append=T)
      cat("mv", "*Schedule*", "BayesTraits_outputs_sch_trees/", "\n", file = shellfile_wd, append=T)
      cat("mv", "*Output.VR.trees*", "BayesTraits_outputs_sch_trees/", "\n", file = shellfile_wd, append=T)
      cat("mv", "*VarRates*", "BayesTraits_outputs_sch_trees/", file = shellfile_wd, append=T)
    }
    if (syst == "windows") {
      cat("Move-Item -Path", "*Log*", "-Destination", "BayesTraits_outputs_log_stone/", "\n", file = shellfile_wd, append=T)
      cat("Move-Item -Path", "*Stones*", "-Destination", "BayesTraits_outputs_log_stone/", "\n", "\n", file = shellfile_wd, append=T)
      cat("Move-Item -Path", "*Schedule*", "-Destination", "BayesTraits_outputs_sch_trees/", "\n", file = shellfile_wd, append=T)
      cat("Move-Item -Path", "*Output.VR.trees*", "-Destination", "BayesTraits_outputs_sch_trees/", "\n", file = shellfile_wd, append=T)
      cat("Move-Item -Path", "*VarRates*", "-Destination", "BayesTraits_outputs_sch_trees/", file = shellfile_wd, append=T)
    }
  }

  #_____________________________________________________________________________
  # Making a copy of the tree, mean, and linked data files within the BayesTraits folder
  if (cc_DataTree) {
    for (i in seq_along(meanfile_name)) {
      invisible(file.copy(from = paste0(meanfile_dir, "/", meanfile_name[i]),
                          to = paste0(BayesTraits_dir, "/", meanfile_name[i])))
      invisible(file.copy(from = paste0(linkfile_dir, "/", linkfile_name[i]),
                          to = paste0(BayesTraits_dir, "/", linkfile_name[i])))
    }
    invisible(file.copy(from = paste0(treefile_dir, "/", treefile_name),
                        to = paste0(BayesTraits_dir, "/", treefile_name)))
  }

}
