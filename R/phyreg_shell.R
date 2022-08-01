#' Create shell scripts to run BayesTraits phylogenetic regression analyses using
#' samples of trait data
#'
#' @author Domingos Cardoso & Matt Lavin
#'
#' @description By reporting the directory paths to the tree and data files, this
#' function produces shell scripts in both Unix and Windows PowerShell formats
#' for running one or multiple phylogenetic regression analyses with Meade & Pagel's
#' (2022) [BayesTraits](http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/BayesTraitsV4.0.0.html)
#' program. In addition to generating the shell script directly in the BayesTraits
#' directory, the function also creates BayesTraits batch files with specified
#' parameters, including generations set for burnin, iterations, and the stepping
#' stones sampler (e.g., bi, it, sa, and st). In addition, models involving alternative
#' tree transformations can be specified. See [BayesTraits V4.0.0 Manual](http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/Files/BayesTraitsV4.0.0-Manual.pdf)
#' for further details involving BayesTraits parameters and default settings for
#' independent contrast regression analysis so as to better understand how to set
#' arguments of this function \code{phyreg.shell}. You might consider using
#' the function \code{\link{phyreg.inputs}} to generate properly formatted
#' input files of mean and linked samples of trait data prior to using \code{phyreg.shell}.
#' See also a more complete article on how to create a [shell script](https://dboslab.github.io/InNOutBT/articles/phyreg_shell.html)
#' with \code{phyreg.shell} for multiple phylogenetic regression analyses.
#'
#' @usage
#' phyreg.shell(meanfile_dir = NULL,
#'              linkfile_dir = NULL,
#'              treefile_dir = NULL,
#'              BayesTraits_dir = NULL,
#'              responvar = NULL,
#'              treetransf = NULL,
#'              Kappa = NULL,
#'              Delta = NULL,
#'              Lambda = NULL,
#'              OU = NULL,
#'              VRLambda = FALSE,
#'              bi = "30000000",
#'              it = "130000000",
#'              sa = "100000",
#'              st = c("500", "50000"),
#'              syst = NULL,
#'              dir_create = "results_BayesTraits_phyreg_shell",
#'              cc_DataTree = TRUE)
#'
#' @param meanfile_dir Path to the folder directory where the data file(s)
#' containing trait means for samples of trait data is(are) stored.
#'
#' @param linkfile_dir Path to the folder directory where the data file(s)
#' containing samples of linked trait data is(are) stored.
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
#' We use *UNI* as a designation for uniform rates of evolution where *VR* is not
#' estimated and no other tree transformation models are estimated (e.g., Lambda
#' is assumed to be 1.0). In this case, "UNI" here set a complete default BayesTraits
#'  regression analysis.
#'
#' @param Kappa If set as default \code{NULL}, this tree transformation model will
#' be estimated, provided that it was also defined in the previous parameter \code{treetransf}.
#' But here you can also set a vector with any specific value between 0 and 1. If
#' you want to run analyses with both specific and estimated values for for *Kappa*,
#' then you have to set here a vector such as c("0.001", "1.0", "Estimate"). Report
#' any number as a character (i.e., in quotes).
#'
#' @param Delta If set as default \code{NULL}, this tree transformation model will
#' be estimated, provided that it was also defined in the previous parameter \code{treetransf}.
#' But here you can also set a vector with any specific value between 0 and 1. If
#' you want to run analyses with both specific and estimated values for for *Delta*,
#' then you have to set here a vector such as c("0.001", "1.0", "Estimate"). Report
#' any number as a character (i.e., in quotes).
#'
#' @param Lambda If set as default \code{NULL}, this tree transformation model will
#' be estimated, provided that it was also defined in the previous parameter \code{treetransf}.
#' But here you can also set a vector with any specific value between 0 and 1. If
#' you want to run analyses with both specific and estimated values for for *Lambda*,
#' then you have to set here a vector such as c("0.001", "1.0", "Estimate"). Report
#' any number as a character (i.e., in quotes).
#'
#' @param OU If set as default \code{NULL}, this tree transformation model will
#' be estimated, provided that it was also defined in the previous parameter \code{treetransf}.
#' But here you can also set a vector with any specific value between 0 and 1. If
#' you want to run analyses with both specific and estimated values for for *OU*,
#' then you have to set here a vector such as c("0.001", "1.0", "Estimate"). Report
#' any number as a character (i.e., in quotes).
#'
#' @param VRLambda Logical, the default is \code{FALSE}. If set as \code{TRUE},
#' then the specific values set for the tree transformation models *Kappa*, *Delta*,
#' *Lambda* or *OU* will be considered as separate runs when using the model of
#' variable rates *VR*, provided that you have also chosen this model in the parameter
#' \code{treetransf}.
#'
#' @param bi Set the number of iterations for the MCMC burnin. Report this number
#' as a character (i.e., in quotes).
#'
#' @param it Set the number of total iterations for the MCMC chain. Report this
#' number as a character (i.e., in quotes).
#'
#' @param sa Set the sample frequency. Report this number as a character (i.e., in quotes).
#'
#' @param st The stepping stone sampler to estimate the marginal likelihood. Report
#' the number of stones to use and the number of iterations for each stone as
#' characters (i.e., a vector of numbers each in quotes).
#'
#' @param syst Report the operating system ("unix" or "windows").
#'
#' @param dir_create Path to the computer's working directory or where the file(s)
#' will be saved. The default setting creates a directory "results_BayesTraits_phyreg_output"
#' in which the results will be saved within a subfolder named by the current date.
#'
#' @param cc_DataTree Logical, the default is \code{TRUE}. A copy of the tree file
#' and the mean and linked data files will be made in the BayesTraits folder where
#' the BayesTraits executable file is located.
#'
#' @seealso \code{\link{phyreg.inputs}}
#' @seealso \code{\link{phyreg.outputs}}
#'
#' @export
#'

phyreg.shell <- function (meanfile_dir = NULL,
                          linkfile_dir = NULL,
                          treefile_dir = NULL,
                          BayesTraits_dir = NULL,
                          responvar = NULL,
                          treetransf = NULL,
                          Kappa = NULL,
                          Delta = NULL,
                          Lambda = NULL,
                          OU = NULL,
                          VRLambda = FALSE,
                          bi = "30000000",
                          it = "130000000",
                          sa = "100000",
                          st = c("500", "50000"),
                          syst = NULL,
                          dir_create = "results_BayesTraits_phyreg_shell",
                          cc_DataTree = TRUE) {

  #_____________________________________________________________________________
  # Create a new directory with current date to save the shell script and run command
  # files at the R directory.
  # If there is no directory... make one
  todaydate <- format(Sys.time(), "%d%b%Y")
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
    shellfile_bt = paste0(BayesTraits_dir, "/run_BayesTraits_phyreg_shell_", syst, ".sh")
    if (!is.null(dir_create)) {
      shellfile_wd = paste0(folder_name, "/run_BayesTraits_phyreg_shell_", syst, ".sh")
    }
  }
  if (syst == "windows"){
    title = "#!/bin/bash"
    shellfile_bt = paste0(BayesTraits_dir, "/run_BayesTraits_phyreg_shell_", syst, ".ps1")
    if (!is.null(dir_create)) {
      shellfile_wd = paste0(folder_name, "/run_BayesTraits_phyreg_shell_", syst, ".ps1")
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

  cat("mkdir", "BayesTraits_phyreg_outputs_log_stone", "\n", file = shellfile_bt, append=T)
  cat("mkdir", "BayesTraits_phyreg_outputs_sch_trees", "\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("mkdir", "BayesTraits_phyreg_outputs_log_stone", "\n", file = shellfile_wd, append=T)
    cat("mkdir", "BayesTraits_phyreg_outputs_sch_trees", "\n", file = shellfile_wd, append=T)
  }

  cat("\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("\n", file = shellfile_wd, append=T)
  }

  treetransf_value <- list(Kappa, Delta, Lambda, OU)
  names(treetransf_value) <- c("Kappa", "Delta", "Lambda", "OU")

  treetransf_complete <- list()
  for (t in seq_along(treetransf_value)) {
    if (!is.null(treetransf_value[[t]])) {
      treetransf_complete[[t]] <- paste0(names(treetransf_value)[t], treetransf_value[[t]])
    }
  }

  tf <- treetransf %in% unique(gsub("Estimate|[[:digit:]]|[.]", "", unlist(treetransf_complete)))
  treetransf_complete <- append(unlist(treetransf_complete), treetransf[!tf])

  if(VRLambda) {
    treetransf_complete <- append(treetransf_complete[(!treetransf_complete %in% "VR")],
                                  paste0("VR", treetransf_complete[grepl("Lambda", treetransf_complete)]))
  }

  param_temp <- list()
  temp <- vector()
  run_temp <- list()
  temp_run <- vector()

  for (l in seq_along(meanfile_name)) {
    explanvar <- gsub("BayesTraits_mean_data_", "", meanfile_name[l])
    explanvar <- gsub("[.]txt", "", explanvar)
    for(i in seq_along(treetransf_complete)) {
      temp[i] <- paste0("paramfile", "=", "run_", responvar, "_", explanvar, "_", treetransf_complete[i], ".txt")

      if (is.null(unlist(treetransf_value))) {
        temp_run[[i]] <- run_commands(treetransf_complete[i], bi, it, sa, st,
                                      linkfile_name[l],
                                      ttransf_value = NULL,
                                      VRLambda=VRLambda)
      }

      if (!is.null(unlist(treetransf_value))) {

        if (treetransf_complete[i] %in% c("UNI", "VR", "Fabric")) {
          temp_run[i] <- run_commands(treetransf_complete[i], bi, it, sa, st,
                                      linkfile_name[l],
                                      ttransf_value = NULL,
                                      VRLambda=VRLambda)
        } else {
          temp_run[i] <- run_commands(treetransf_complete[i], bi, it, sa, st,
                                      linkfile_name[l],
                                      gsub("VRLambda|Kappa|Delta|Lambda|OU", "",
                                           treetransf_complete[i]),
                                      VRLambda=VRLambda)
        }

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

  if (any(grepl("VR", treetransf_complete))) {
    temp_outrees <- vector()
    temp_varates <- vector()
  }

  for(i in seq_along(meanfile_name)) {
    for (l in seq_along(treetransf_complete)) {
      if (syst == "unix") {
        temp_logs[l] <- paste0("mv ", meanfile_name[i], ".Log.txt ", meanfile_name[i], ".Log.", treetransf_complete[l], ".txt")
        temp_stones[l] <- paste0("mv ", meanfile_name[i], ".Stones.txt ", meanfile_name[i], ".Stones.", treetransf_complete[l], ".txt")
        temp_scheds[l] <- paste0("mv ", meanfile_name[i], ".Schedule.txt ", meanfile_name[i], ".Schedule.", treetransf_complete[l], ".txt")
      }
      if (syst == "windows") {
        temp_logs[l] <- paste0("Rename-Item -Path ", meanfile_name[i], ".Log.txt ", "-NewName ", meanfile_name[i], ".Log.", treetransf_complete[l], ".txt")
        temp_stones[l] <- paste0("Rename-Item -Path ", meanfile_name[i], ".Stones.txt ", "-NewName ", meanfile_name[i], ".Stones.", treetransf_complete[l], ".txt")
        temp_scheds[l] <- paste0("Rename-Item -Path ", meanfile_name[i], ".Schedule.txt ", "-NewName ", meanfile_name[i], ".Schedule.", treetransf_complete[l], ".txt")
      }
      if (treetransf_complete[l] == "VR" |
          grepl("VRLambda", treetransf_complete[l])) {
        if (syst == "unix") {
          temp_outrees[l] <- paste0("mv ", meanfile_name[i], ".Output.trees ", meanfile_name[i], ".Output.", treetransf_complete[l], ".trees")
          temp_varates[l] <- paste0("mv ", meanfile_name[i], ".VarRates.txt ", meanfile_name[i], ".VarRates.", treetransf_complete[l], ".txt")
        }
        if (syst == "windows") {
          temp_outrees[l] <- paste0("Rename-Item -Path ", meanfile_name[i], ".Output.trees ", "-NewName ", meanfile_name[i], ".Output.", treetransf_complete[l], ".trees")
          temp_varates[l] <- paste0("Rename-Item -Path ", meanfile_name[i], ".VarRates.txt ", "-NewName ", meanfile_name[i], ".VarRates.", treetransf_complete[l], ".txt")
        }
      }
    }
    logs[[i]] <- temp_logs
    stones[[i]] <- temp_stones
    scheds[[i]] <- temp_scheds
    if (any(treetransf_complete == "VR") |
        any(grepl("VRLambda", treetransf_complete))) {
      outrees[[i]] <- temp_outrees
      varates[[i]] <- temp_varates
    }
  }

  temp <- append(unlist(logs), unlist(stones))
  temp <- append(temp, unlist(scheds))
  if (any(treetransf_complete == "VR") |
      any(grepl("VRLambda", treetransf_complete))) {
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
  mods <- gsub(".+Stones[.]|.+Log[.]|.+Schedule[.]|.+Output[.]|.+VarRates[.]", "", mods)

  change_names <- list()
  temp <- list()
  for (i in seq_along(unique(vars))) {
    tf <- vars %in% unique(vars)[i]
    for (l in seq_along(treetransf_complete)) {
      temp[[l]] <- allnames_tochange[tf][mods[tf] %in% treetransf_complete[l]]
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
    cat("mv", "*Log*", "BayesTraits_phyreg_outputs_log_stone/", "\n", file = shellfile_bt, append=T)
    cat("mv", "*Stones*", "BayesTraits_phyreg_outputs_log_stone/", "\n", "\n", file = shellfile_bt, append=T)
    cat("mv", "*Schedule*", "BayesTraits_phyreg_outputs_sch_trees/", "\n", file = shellfile_bt, append=T)
    cat("mv", "*Output.VR*", "BayesTraits_phyreg_outputs_sch_trees/", "\n", file = shellfile_bt, append=T)
    cat("mv", "*VarRates*", "BayesTraits_phyreg_outputs_sch_trees/", file = shellfile_bt, append=T)
  }
  if (syst == "windows") {
    cat("Move-Item -Path", "*Log*", "-Destination", "BayesTraits_phyreg_outputs_log_stone/", "\n", file = shellfile_bt, append=T)
    cat("Move-Item -Path", "*Stones*", "-Destination", "BayesTraits_phyreg_outputs_log_stone/", "\n", "\n", file = shellfile_bt, append=T)
    cat("Move-Item -Path", "*Schedule*", "-Destination", "BayesTraits_phyreg_outputs_sch_trees/", "\n", file = shellfile_bt, append=T)
    cat("Move-Item -Path", "*Output.VR*", "-Destination", "BayesTraits_phyreg_outputs_sch_trees/", "\n", file = shellfile_bt, append=T)
    cat("Move-Item -Path", "*VarRates*", "-Destination", "BayesTraits_phyreg_outputs_sch_trees/", file = shellfile_bt, append=T)
  }
  if (!is.null(dir_create)) {
    if (syst == "unix") {
      cat("mv", "*Log*", "BayesTraits_phyreg_outputs_log_stone/", "\n", file = shellfile_wd, append=T)
      cat("mv", "*Stones*", "BayesTraits_phyreg_outputs_log_stone/", "\n", "\n", file = shellfile_wd, append=T)
      cat("mv", "*Schedule*", "BayesTraits_phyreg_outputs_sch_trees/", "\n", file = shellfile_wd, append=T)
      cat("mv", "*Output.VR*", "BayesTraits_phyreg_outputs_sch_trees/", "\n", file = shellfile_wd, append=T)
      cat("mv", "*VarRates*", "BayesTraits_phyreg_outputs_sch_trees/", file = shellfile_wd, append=T)
    }
    if (syst == "windows") {
      cat("Move-Item -Path", "*Log*", "-Destination", "BayesTraits_phyreg_outputs_log_stone/", "\n", file = shellfile_wd, append=T)
      cat("Move-Item -Path", "*Stones*", "-Destination", "BayesTraits_phyreg_outputs_log_stone/", "\n", "\n", file = shellfile_wd, append=T)
      cat("Move-Item -Path", "*Schedule*", "-Destination", "BayesTraits_phyreg_outputs_sch_trees/", "\n", file = shellfile_wd, append=T)
      cat("Move-Item -Path", "*Output.VR*", "-Destination", "BayesTraits_phyreg_outputs_sch_trees/", "\n", file = shellfile_wd, append=T)
      cat("Move-Item -Path", "*VarRates*", "-Destination", "BayesTraits_phyreg_outputs_sch_trees/", file = shellfile_wd, append=T)
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

  #_____________________________________________________________________________
  # Remove running files from the folder where the BayesTraits executable file is stored
  cat("\n", "\n", "\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("\n", "\n", "\n", file = shellfile_wd, append=T)
  }

  cat("## remove running files from the BayesTraits folder\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("## remove running files from the BayesTraits folder\n", file = shellfile_wd, append=T)
  }

  if (syst == "unix") {
    cat(paste0("rm $treefile"), "\n", file = shellfile_bt, append=T)
  }
  if (syst == "windows") {
    cat(paste0("Remove-Item $treefile"), "\n", file = shellfile_bt, append=T)
  }
  if (!is.null(dir_create)) {
    if (syst == "unix") {
      cat(paste0("rm $treefile"), "\n", file = shellfile_wd, append=T)
    }
    if (syst == "windows") {
      cat(paste0("Remove-Item $treefile"), "\n", file = shellfile_wd, append=T)
    }
  }

  cat("\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("\n", file = shellfile_wd, append=T)
  }

  for(i in seq_along(meanfiles_nbrs)) {
    if (syst == "unix") {
      cat(paste0("rm $", meanfiles_nbrs[i]), "\n", file = shellfile_bt, append=T)
    }
    if (syst == "windows") {
      cat(paste0("Remove-Item $", meanfiles_nbrs[i]), "\n", file = shellfile_bt, append=T)
    }
    if (!is.null(dir_create)) {
      if (syst == "unix") {
        cat(paste0("rm $", meanfiles_nbrs[i]), "\n", file = shellfile_wd, append=T)
      }
      if (syst == "windows") {
        cat(paste0("Remove-Item $", meanfiles_nbrs[i]), "\n", file = shellfile_wd, append=T)
      }
    }
  }

  cat("\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("\n", file = shellfile_wd, append=T)
  }

  for(i in seq_along(param_temp)) {
    if (syst == "unix") {
      cat(paste0("rm $paramfile", i), "\n", file = shellfile_bt, append=T)
    }
    if (syst == "windows") {
      cat(paste0("Remove-Item $paramfile", i), "\n", file = shellfile_bt, append=T)
    }
    if (!is.null(dir_create)) {
      if (syst == "unix") {
        cat(paste0("rm $paramfile", i), "\n", file = shellfile_wd, append=T)
      }
      if (syst == "windows") {
        cat(paste0("Remove-Item $paramfile", i), "\n", file = shellfile_wd, append=T)
      }
    }
  }

  cat("\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("\n", file = shellfile_wd, append=T)
  }

  if (syst == "unix") {
    cat(paste0(" rm *linked_data*"), "\n", file = shellfile_bt, append=T)
  }
  if (syst == "windows") {
    cat(paste0("Remove-Item *linked_data*"), "\n", file = shellfile_bt, append=T)
  }
  if (!is.null(dir_create)) {
    if (syst == "unix") {
      cat(paste0("rm  *linked_data*"), "\n", file = shellfile_wd, append=T)
    }
    if (syst == "windows") {
      cat(paste0("Remove-Item *linked_data*"), "\n", file = shellfile_wd, append=T)
    }
  }

  cat("\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("\n", file = shellfile_wd, append=T)
  }

  if (syst == "unix") {
    cat(paste0("rm run_BayesTraits_phyreg_shell_unix.sh"), "\n", file = shellfile_bt, append=T)
  }
  if (syst == "windows") {
    cat(paste0("Remove-Item run_BayesTraits_phyreg_shell_windows.ps1"), "\n", file = shellfile_bt, append=T)
  }
  if (!is.null(dir_create)) {
    if (syst == "unix") {
      cat(paste0("rm run_BayesTraits_phyreg_shell_unix.sh"), "\n", file = shellfile_wd, append=T)
    }
    if (syst == "windows") {
      cat(paste0("Remove-Item run_BayesTraits_phyreg_shell_windows.ps1"), "\n", file = shellfile_wd, append=T)
    }
  }

  cat("\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat("\n", file = shellfile_wd, append=T)
  }

  cat(paste0("echo ", "\"", "The input tree file $treefile and all associated running linked and mean trait data were removed from BayesTraits folder.", "\""), "\n", file = shellfile_bt, append=T)
  if (!is.null(dir_create)) {
    cat(paste0("echo ", "\"", "The input tree file $treefile and all associated running linked and mean trait data were removed from BayesTraits folder.", "\""), "\n", file = shellfile_wd, append=T)
  }

}


# Auxiliary function to generate the run commands
run_commands <- function (ttransf, bi, it, sa, st, linkfile, ttransf_value, VRLambda=VRLambda) {

  ttransf <- gsub("Estimate|[[:digit:]]|[.]", "", ttransf)

  if (ttransf == "VR") {
    ttransf <- tolower(ttransf)
  }

  if (ttransf == "UNI") {
    run_res <- paste0("9", "\n",
                      "2", "\n",
                      "bi ", bi, "\n",
                      "it ", it, "\n",
                      "sa ", sa, "\n",
                      "stones ", st[1], " ", st[2], "\n",
                      "DistData ",
                      linkfile, "\n",
                      "Run")
  }

  if (ttransf != "UNI" & is.null(ttransf_value)) {
    run_res <- paste0("9", "\n",
                      "2", "\n",
                      ttransf, "\n",
                      "bi ", bi, "\n",
                      "it ", it, "\n",
                      "sa ", sa, "\n",
                      "stones ", st[1], " ", st[2], "\n",
                      "DistData ",
                      linkfile, "\n",
                      "Run")
  }

  if (ttransf != "UNI" & !is.null(ttransf_value)) {

    if (VRLambda) {

      if (ttransf == "VRLambda") {

        if (ttransf_value == "Estimate") {
          run_res <- paste0("9", "\n",
                            "2", "\n",
                            "vr", "\n",
                            "Lambda", "\n",
                            "bi ", bi, "\n",
                            "it ", it, "\n",
                            "sa ", sa, "\n",
                            "stones ", st[1], " ", st[2], "\n",
                            "DistData ",
                            linkfile, "\n",
                            "Run")
        } else {

          run_res <- paste0("9", "\n",
                            "2", "\n",
                            "vr", "\n",
                            "Lambda", " ", ttransf_value, "\n",
                            "bi ", bi, "\n",
                            "it ", it, "\n",
                            "sa ", sa, "\n",
                            "stones ", st[1], " ", st[2], "\n",
                            "DistData ",
                            linkfile, "\n",
                            "Run")
        }

      } else {

        if (ttransf_value == "Estimate") {
          run_res <- paste0("9", "\n",
                            "2", "\n",
                            ttransf, "\n",
                            "bi ", bi, "\n",
                            "it ", it, "\n",
                            "sa ", sa, "\n",
                            "stones ", st[1], " ", st[2], "\n",
                            "DistData ",
                            linkfile, "\n",
                            "Run")
        } else {
          run_res <- paste0("9", "\n",
                            "2", "\n",
                            ttransf, " ", ttransf_value, "\n",
                            "bi ", bi, "\n",
                            "it ", it, "\n",
                            "sa ", sa, "\n",
                            "stones ", st[1], " ", st[2], "\n",
                            "DistData ",
                            linkfile, "\n",
                            "Run")
        }

      }
    }

    if (!VRLambda) {

      if (ttransf_value == "Estimate") {
        run_res <- paste0("9", "\n",
                          "2", "\n",
                          ttransf, "\n",
                          "bi ", bi, "\n",
                          "it ", it, "\n",
                          "sa ", sa, "\n",
                          "stones ", st[1], " ", st[2], "\n",
                          "DistData ",
                          linkfile, "\n",
                          "Run")
      } else {
        run_res <- paste0("9", "\n",
                          "2", "\n",
                          ttransf, " ", ttransf_value, "\n",
                          "bi ", bi, "\n",
                          "it ", it, "\n",
                          "sa ", sa, "\n",
                          "stones ", st[1], " ", st[2], "\n",
                          "DistData ",
                          linkfile, "\n",
                          "Run")
      }

    }
  }

  return(run_res)
}

