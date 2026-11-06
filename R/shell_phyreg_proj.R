#' Create a BayesTraits project folder with associated example scripts and input
#' data for a full shell script run and log processing of BayesTraits phylogenetic
#' regression analyses
#'
#' @author Domingos Cardoso & Matt Lavin
#'
#' @description By reporting a directory path and the computer's operating system
#' ("unix" or "windows"), this function creates a BayesTraits project folder with
#' associated example scripts and input data for a full shell script run and log
#' processing of multiple phylogenetic regression analyses with Meade & Pagel's
#' (2022) [BayesTraits](https://www.evolution.reading.ac.uk/BayesTraitsV5.0.3/BayesTraitsV5.0.3.html)
#' program. The function produces an example of a shell script in Unix or Windows PowerShell
#' that helps to fully automate both the complete BayesTraits regression analysis
#' and the processing of the resulting log and stones files, by automatically running
#' the three main \code{InNOutBT} functions \code{\link{phyreg.inputs}},
#' \code{\link{phyreg.shell}}, and \code{\link{phyreg.outputs}}, without
#' directly opening the R or RStudio programs. Note that these files are just
#' examples from a specific study of the Vataireoid legumes, and so you will have
#' to edit the folder paths and analysis parameters in these example scripts according
#' to your desired regression analyses. See also a more complete article
#' [describing a BayesTraits shell run project](https://dboslab.github.io/InNOutBT/articles/phyreg_shell_project.html)
#' that can be created with \code{shell.phyreg.proj} for fully automated
#' multiple phylogenetic regression analyses.
#'
#' @details
#' After creating the project folder with this function, follow these steps to
#' run the project. These steps are for macOS (\code{syst = "unix"}); see the
#' article linked above for the steps in Windows PowerShell.
#'
#' **Before the first run**
#'
#' 1. Install the BayesTraits program in any folder and note the path to that
#'    folder and the name of the executable file (e.g., \code{BayesTraitsV5}).
#' 2. Check that the executable file starts by running its full path in the Terminal.
#'    If macOS stops it (e.g., the message \code{killed}), remove the block that macOS
#'    adds to downloaded files, provided that you trust the source of the download:
#'    \code{xattr -d com.apple.quarantine /path/to/BayesTraits/BayesTraitsV5}.
#' 3. Check that the R packages \code{InNOutBT} and \code{treeio} are installed, and
#'    that \code{pandoc} is available in the Terminal, as the Word table needs it:
#'    \code{Rscript -e 'rmarkdown::pandoc_available()'}.
#'
#' **Edit the project files in any text editor**
#'
#' 4. In \code{BayesTraits_full_phyreg_shell_unix.sh}, set \code{PROJECTFOLDER} to the
#'    path of the project folder, \code{BAYESTRFOLDER} to the path of the folder that
#'    contains the BayesTraits executable file, and \code{EMAIL} to your e-mail address.
#' 5. In \code{2_producing_BayesTraits_phyreg_shell.R}, set \code{BayesTraits_dir} to the
#'    same path as \code{BAYESTRFOLDER}, and \code{BayesTraits_exe} to the name of the
#'    executable file.
#' 6. To run your own study, also edit the three numbered R scripts and replace the
#'    files in the \code{Data} folder.
#'
#' **Run the project in the Terminal**
#'
#' 7. Move to the project folder, not to the BayesTraits folder, because the main shell
#'    script moves to the BayesTraits folder by itself:
#'    \code{cd /path/to/BayesTraits_shell_phyreg_project}.
#' 8. Allow the main shell script to run, which is needed only once:
#'    \code{chmod 777 BayesTraits_full_phyreg_shell_unix.sh}.
#' 9. Start the full analysis: \code{./BayesTraits_full_phyreg_shell_unix.sh}. Keep the
#'    Terminal open until it finishes, which can take a long time with the default
#'    number of iterations.
#'
#' The script then produces the input files, creates and runs the BayesTraits
#' analyses, moves the resulting log and stones files to the project folder, and
#' processes them. The final table and heatmaps are saved in the project folder, in
#' \code{results_BayesTraits_phyreg_output}, within a subfolder named with the current date.
#'
#' @usage
#' shell.phyreg.proj(dir_create_proj = NULL,
#'                   syst = NULL)
#'
#' @param dir_create_proj Path to any directory on your computer where an example
#' of a BayesTraits shell project folder and associated files for a full run with
#' shell script should be saved.
#'
#' @param syst Report the operating system ("unix" or "windows").
#'
#' @seealso \code{\link{phyreg.inputs}}
#' @seealso \code{\link{phyreg.shell}}
#' @seealso \code{\link{phyreg.outputs}}
#'
#' @examples
#' \dontrun{
#' library(InNOutBT)
#'
#' shell.phyreg.proj(dir_create_proj = "/Users/domingoscardoso",
#'                   syst = "unix")
#'}
#'
#' @export
#'

shell.phyreg.proj <- function (dir_create_proj = NULL,
                               syst = NULL) {

  if (is.null(syst)) {
    syst <- .Platform$OS.type
  }


  #_____________________________________________________________________________
  # Create new folders
  dir.create(paste0(dir_create_proj, "/", "BayesTraits_shell_phyreg_project"))
  dir.create(paste0(dir_create_proj, "/", "BayesTraits_shell_phyreg_project/Data"))


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
  # Making a copy of the files into the newly created BayesTraits_shell_phyreg_project

  invisible(file.copy(from = paste0(systdir, "/", tree),
                      to = paste0(dir_create_proj, "/",
                                  "BayesTraits_shell_phyreg_project/Data/", tree)))
  invisible(file.copy(from = paste0(systdir, "/", dataset),
                      to = paste0(dir_create_proj, "/",
                                  "BayesTraits_shell_phyreg_project/Data/", dataset)))

  if (syst == "unix") {
    for (i in seq_along(uniscripts)) {
      invisible(file.copy(from = paste0(systdir, "/", uniscripts[i]),
                          to = paste0(dir_create_proj, "/", "BayesTraits_shell_phyreg_project/",
                                      gsub("_unix", "", uniscripts[i]))))
    }
    invisible(file.copy(from = paste0(systdir, "/", unishell),
                        to = paste0(dir_create_proj, "/", "BayesTraits_shell_phyreg_project/",
                                    unishell)))
  }

  if (syst == "windows") {
    for (i in seq_along(winscripts)) {
      invisible(file.copy(from = paste0(systdir, "/", winscripts[i]),
                          to = paste0(dir_create_proj, "/", "BayesTraits_shell_phyreg_project/",
                                      gsub("_windows", "", winscripts[i]))))
    }
    invisible(file.copy(from = paste0(systdir, "/", winshell),
                        to = paste0(dir_create_proj, "/", "BayesTraits_shell_phyreg_project/",
                                    winshell)))
  }

}

