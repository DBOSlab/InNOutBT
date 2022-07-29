#' Create a BayesTraits project folder with associated example scripts and input
#' data for a full shell script run and log processing of BayesTraits phylogenetic
#' regression analyses
#'
#' @author Domingos Cardoso & Matt Lavin
#'
#' @description By reporting a directory path and the computer's operating system
#' ("unix" or "windows"), this function creates a BayesTraits project folder with
#' associated example scripts and input data for a full shell script run and log
#' processing multiple phylogenetic regression analyses with Meade & Pagel's
#' (2022) [BayesTraits](http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/BayesTraitsV4.0.0.html)
#' program. The function produces an example of a shell script in Unix or Windows PowerShell
#' that helps to fully automatize both the complete BayesTraits regression analysis
#' and processing of the resulting log and stones files, by automatically running
#' the three main \code{InNOutBT} functions \code{\link{BayesTraits.inputs}},
#' \code{\link{BayesTraits.shell}}, and \code{\link{BayesTraits.outputs}}, without
#' directly opening the R or R Studio programs. Note that these files are just
#' examples from an specific study of the Vataireoid legumes, and so you will have
#' to edit the folder paths and analysis parameters in such example scripts according
#' to your desired regression analyses. See also a more complete article
#' [describing a BayesTraits shell run project](https://dboslab.github.io/InNOutBT/articles/independent_contrast_regression_shell_project.html)
#' that can be created with \code{BayesTraits.shell.reg.proj} for fully automated
#' multiple phylogenetic regression analyses.
#'
#' @usage
#' BayesTraits.shell.reg.proj(dir_create_proj = NULL,
#'                            syst = NULL)
#'
#' @param dir_create_proj Path to any computer's directory where an example
#' of a BayesTraits shell project folder and associated files for a full run with
#' shell script should be saved.
#'
#' @param syst Report the operating system ("unix" or "windows").
#'
#' @seealso \code{\link{BayesTraits.inputs}}
#' @seealso \code{\link{BayesTraits.shell}}
#' @seealso \code{\link{BayesTraits.outputs}}
#'
#' @examples
#' \dontrun{
#' library(InNOutBT)
#'
#' BayesTraits.shell.reg.proj(dir_create_proj = "/Users/domingoscardoso/Downloads",
#'                            syst = "unix")
#'}
#'
#' @export
#'

BayesTraits.shell.reg.proj <- function (dir_create_proj = NULL,
                                        syst = NULL) {

  if (is.null(syst)) {
    syst <- .Platform$OS.type
  }


  #_____________________________________________________________________________
  # Create new folders
  dir.create(paste0(dir_create_proj, "/", "BayesTraits_shell_run_project"))
  dir.create(paste0(dir_create_proj, "/", "BayesTraits_shell_run_project/Data"))


  #_____________________________________________________________________________
  # Get name of the files
  systdir <- system.file("extdata", package = "InNOutBT")
  allfiles <- list.files(system.file("extdata", package = "InNOutBT"))

  tree <- allfiles[grepl("[.]tree", allfiles)]
  dataset <- allfiles[grepl("[.]csv", allfiles)]
  uniscripts <- allfiles[grepl("unix[.]R", allfiles)]
  unishell <- allfiles[grepl("[.]sh$", allfiles)]
  winscripts <- allfiles[grepl("windows[.]R", allfiles)]
  winshell <- allfiles[grepl("[.]ps1$", allfiles)]


  #_____________________________________________________________________________
  # Making a copy of the files into the newly created BayesTraits_shell_run_project

  invisible(file.copy(from = paste0(systdir, "/", tree),
                      to = paste0(dir_create_proj, "/",
                                  "BayesTraits_shell_run_project/Data/", tree)))
  invisible(file.copy(from = paste0(systdir, "/", dataset),
                      to = paste0(dir_create_proj, "/",
                                  "BayesTraits_shell_run_project/Data/", dataset)))

  if (syst == "unix") {
    for (i in seq_along(uniscripts)) {
      invisible(file.copy(from = paste0(systdir, "/", uniscripts[i]),
                          to = paste0(dir_create_proj, "/", "BayesTraits_shell_run_project/",
                                      gsub("_unix", "", uniscripts[i]))))
    }
    invisible(file.copy(from = paste0(systdir, "/", unishell),
                        to = paste0(dir_create_proj, "/", "BayesTraits_shell_run_project/",
                                    unishell)))
  }

  if (syst == "windows") {
    for (i in seq_along(winscripts)) {
      invisible(file.copy(from = paste0(systdir, "/", winscripts[i]),
                          to = paste0(dir_create_proj, "/", "BayesTraits_shell_run_project/",
                                      gsub("_windows", "", winscripts[i]))))
    }
    invisible(file.copy(from = paste0(systdir, "/", unishell),
                        to = paste0(dir_create_proj, "/", "BayesTraits_shell_run_project/",
                                    winshell)))
  }

}

