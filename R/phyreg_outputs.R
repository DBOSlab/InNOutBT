#' Process BayesTraits phylogenetic regression output data
#'
#' @author Domingos Cardoso & Matt Lavin
#'
#' @description This function processes the output reported in the log and stones .txt files,
#' which are generated from the phylogenetic regression analyses in Meade & Pagel's (2022)
#' [BayesTraits](http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/BayesTraitsV4.0.0.html)
#' program. The resulting files are tables either Word or CSV formats and these rank
#' competing models by their likelihood scores. Bayes Factors can also be reported in the
#' case where models with tree transformations (e.g., either Lambda or variable rates are
#' estimated) always have high likelihood scores than simpler models (e.g., those with neither
#' Lambda nor variable rates estimated). In this case, the reported Bayes Factors pertain
#' to the more complex model listed in the table and the simpler nested model that didn't
#' involve a tree transformation. Heatmaps based on pairwise Bayes Factor comparison
#' are also automatically produced in PDF format. These are intended to identify
#' models that are most similar to the model with the highest likelihood score.
#' See [BayesTraits V4.0.0 Manual](http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/Files/BayesTraitsV4.0.0-Manual.pdf)
#' for details of the output files produced by BayesTraits, which may facilitate an
#' understanding of how output files are processed with the function \code{phyreg.outputs}.
#' See also a more complete article on how to use \code{phyreg.outputs} to
#' [process the resulting output files](https://dboslab.github.io/InNOutBT/articles/phyreg_outputs.html)
#' from the phylogenetic regression analyses.
#'
#' @usage
#' phyreg.outputs(logst_dir = NULL,
#'                responvar = NULL,
#'                explanvar = NULL,
#'                treetransf = c("Kappa", "Delta", "Lambda", "OU", "UNI", "VR", "Fabric"),
#'                bayesfactor = FALSE,
#'                unirates = TRUE,
#'                outformat = c("Word", "CSV"),
#'                tableleg = "Phylogenetic regression in BayesTraits: models and coefficients.",
#'                dir_create = "results_BayesTraits_phyreg_output",
#'                outfile = "BayesTraits_phyreg_output_table",
#'                ...)
#'
#' @param logst_dir Path to the directory where are stored the log and stepping stones files
#' generated during the independent contrast regression analysis. The log file(s) contain
#' the model options and output. The stepping stones file(s) contain the marginal log
#' likelihood for a given model. These files are labeled with the explanatory variable(s)
#' of the relevant model, where each variable is separated by an underscore. Explanatory
#' variables include any tree transformation included in a model, and these can include
#' *Kappa*, *Delta*, *Lambda*, *OU* (Ornstein–Uhlenbeck), *UNI* (for no tree transformation),
#' *VR* (variable rates) or *Fabric* (for a model that accommodates uneven evolutionary
#' landscape, as described by [Pagel et al. (2022)](https://doi.org/10.1038/s41467-022-28595-z).
#' Thus, Log file(s) have a format like "BayesTraits_mean_data_bio12_bio15_VR.log.txt" and
#' "BayesTraits_mean_data_bio12_bio15_nnodes_Lambda.log.txt" (where nnodes signifies net
#' nodes). This indicates that the first log file contains three explanatory variables
#' ("bio12", "bio15", and variable rates), and the second log file contains four explanatory
#' variables ("bio12", "bio15", "nnodes", and Lambda). Stepping stones file(s) should be
#' named correspondingly (e.g., "BayesTraits_mean_data_bio12_bio15_VR.stones.txt" and
#' "BayesTraits_mean_data_bio12_bio15_nnodes_Lambda.stones.txt".)
#'
#' @param responvar Report the name of the response variable.
#'
#' @param explanvar Report the name(s) of the explanatory variable(s) exactly as
#' they are also written in the names of the log and stones files being processed.
#'
#' @param treetransf Select any of the available BayesTraits tree transformations
#' (*Kappa*, *Delta*, *Lambda*, *OU*, *UNI*, *VR* or *Fabric*) that are included
#' among the competing models.
#'
#' @param bayesfactor Logical, if \code{FALSE}, then no Bayes Factor pairwise comparisons
#' will performed among the models.
#'
#' @param unirates Logical, if \code{FALSE}, then models lacking tree transformations
#' (*Kappa*, *Delta*, *Lambda* or *OU*) and variable rates and Fabric (i.e., default values
#' and uniform rates of evolution are assumed) will not be displayed in the main table.
#'
#' @param outformat Define either "Word" or "CSV", or both, a vector c("Word", "CSV"),
#' for writing the results in such formats.
#'
#' @param tableleg A legend for the main table, provided that you have also chosen
#' \code{outformat = "Word"}. If you do not report a legend, then the Word-formatted
#' table will have the default legend **Phylogenetic regression in BayesTraits: models and coefficients**.
#'
#' @param dir_create Path to the directory where the file(s) will be saved. the default
#' setting creates a directory named **results_BayesTraits_phyreg_output** where the results
#' will be saved in a subfolder named by the current date.
#'
#' @param outfile Name of the resulting table-formatted files in either "Word", "CSV",
#' or both (depending on the specified argument \code{outformat}). If no name is reported,
#' the default setting creates a file named
#' *BayesTraits_phyreg_output_table.docx* and/or *BayesTraits_phyreg_output_table.csv*.
#'
#' @param ... Additional parameters passed to pdf.
#'
#' @return Table in dataframe format, which is also saved in .csv format.
#'
#' @seealso \code{\link{phyreg.inputs}}
#'
#' @importFrom dplyr arrange
#' @importFrom magrittr "%>%"
#' @importFrom tibble tibble add_column
#' @importFrom combinat combn
#' @importFrom reshape2 melt
#' @importFrom wesanderson wes_palette
#' @importFrom pheatmap pheatmap
#' @importFrom grid gpar unit
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_x_discrete scale_y_discrete scale_fill_gradientn theme_minimal theme element_text xlab ylab
#' @importFrom rmarkdown word_document render
#' @importFrom knitr kable
#'
#' @export
#'

phyreg.outputs <- function(logst_dir = NULL,
                           responvar = NULL,
                           explanvar = NULL,
                           treetransf = c("Kappa", "Delta", "Lambda", "OU", "UNI", "VR", "Fabric"),
                           bayesfactor = FALSE,
                           unirates = TRUE,
                           outformat = c("Word", "CSV"),
                           tableleg = "Phylogenetic regression in BayesTraits: models and coefficients.",
                           dir_create = "results_BayesTraits_phyreg_output",
                           outfile = "BayesTraits_phyreg_output_table",
                           ...) {

  # Input data
  tf <- grepl("Log", list.files(logst_dir))
  log_infiles <- list.files(logst_dir)[tf]
  tf <- grepl("Stones", list.files(logst_dir))
  stn_infiles <- list.files(logst_dir)[tf]

  fnames <- gsub("BayesTraits_mean_data_|[.]txt.Log[.]|[.]txt$", "_", log_infiles)
  fnames <- gsub("^_|_$", "", fnames)
  fnames <- gsub("_Lmda", "_Lambda", fnames)

  # Differentiating the models of tree transformation "Kappa", "Delta", "Lambda",
  # and "OU" from the remaining, i.e. "Fabric", "VR", and "UNI".
  treetransf <- treetransf[!treetransf %in% c("Fabric", "VR", "UNI")]
  if(length(treetransf) == 0) treetransf <- NULL


  #_____________________________________________________________________________
  # Getting the Log marginal likelihood
  print("Analyzing Log marginal likelihood within the Stones files... ")
  stones_files <- list()
  logL <- list()
  for (i in seq_along(stn_infiles)) {
    stones_files[[i]] <- readLines(paste0(logst_dir, "/", stn_infiles[i]))
    logL[[i]] <- stones_files[[i]][grepl("Log marginal likelihood[:]\t", stones_files[[i]])]
    logL[[i]] <- as.numeric(gsub(".+\t", "", logL[[i]]))
  }
  names(logL) <- fnames


  #_____________________________________________________________________________
  # Processing the Log files
  BayesTraits_outs <- list()
  for (i in seq_along(log_infiles)) {

    print(paste0("Loading... ", log_infiles[i]))

    BayesTraits_outs[[i]] <- readLines(paste0(logst_dir, "/", log_infiles[i]))

    ntaxa <- BayesTraits_outs[[i]][grepl("No Taxa", BayesTraits_outs[[i]])]
    ntaxa <- as.numeric(gsub(".+\\s", "", ntaxa))
    nrates <- BayesTraits_outs[[i]][grepl("No of Rates", BayesTraits_outs[[i]])]
    nrates <- as.numeric(gsub(".+\\s", "", nrates))

    trtransf <- BayesTraits_outs[[i]][grepl("Estimate", BayesTraits_outs[[i]])]
    trtransf <- gsub("[:].+", "", trtransf)
    trtransf_value <- ifelse(length(trtransf)>0, 0, -1)
    VR <- BayesTraits_outs[[i]][grepl("VRBL", BayesTraits_outs[[i]])]
    VR <- ifelse(length(VR)>0, 10, 0)
    Fabric <- BayesTraits_outs[[i]][grepl("VR_LS_BL", BayesTraits_outs[[i]])]
    Fabric <- ifelse(length(Fabric)>0, 10, 0)

    skiplines <- (34+ntaxa+(2*nrates))+trtransf_value+VR+Fabric

    BayesTraits_outs[[i]] <- read.delim(paste0(logst_dir, "/", log_infiles[i]),
                                        skip = skiplines, sep="\t", header=TRUE)
  }
  names(BayesTraits_outs) <- fnames


  mean_pMCMC <- list()
  for (i in seq_along(log_infiles)) {

    betas <- names(BayesTraits_outs[[i]])[grepl("^Beta[.][[:digit:]]", names(BayesTraits_outs[[i]]))]
    MCMCsamples <- length(BayesTraits_outs[[i]][[1]])

    if (!is.null(treetransf)) {
      df <- data.frame(param = c(betas, treetransf, "Rsquared", paste0("mean", responvar)),
                       mean = rep(NA, length(betas)+(length(treetransf)+2)),
                       pMCMC = NA)
    } else {
      df <- data.frame(param = c(betas, treetransf, "Rsquared", paste0("mean", responvar)),
                       mean = rep(NA, length(betas)+2),
                       pMCMC = NA)
    }


    for (l in seq_along(betas)) {
      df$mean[l] <- round(mean(BayesTraits_outs[[i]][[betas[l]]]), digits=5)
      if (df$mean[l] >= 0) {
        df$pMCMC[l] <- length(which(BayesTraits_outs[[i]][[betas[l]]] <= 0))/MCMCsamples
      }
      if (df$mean[l] <= 0) {
        df$pMCMC[l] <- length(which(BayesTraits_outs[[i]][[betas[l]]] >= 0))/MCMCsamples
      }
    }

    if (any(treetransf %in% names(BayesTraits_outs[[i]]))) {
      treetransf_temp <- gsub(".+_|[[:digit:]]|[.]|Estimate|VR", "", fnames[i])
      df$mean[df$param %in% treetransf_temp] <- round(mean(BayesTraits_outs[[i]][[treetransf_temp]]), digits=5)
      df$mean[df$param %in% "Rsquared"] <- round(mean(BayesTraits_outs[[i]][["R.2"]]), digits=5)
      df$mean[df$param %in% paste0("mean", responvar)] <- round(10^(mean(BayesTraits_outs[[i]][["Alpha"]])), digits=2)

      # Adding "0" for OU tree transformation and "1.0" for all other tree transformation if they were not estimated
      df$mean[is.na(df$mean)] <- ifelse(df$param[is.na(df$mean)][df$param[is.na(df$mean)] %in% treetransf] %in% "OU", "0", "1.0")

      # Adding back the value of tree transformation when not estimated
      if (gsub(".+_|[[:alpha:]]", "", fnames[i]) != "") {
        df$mean[df$param %in% gsub(".+_|[[:digit:]]|[.]|VR", "", fnames[i])] <- gsub(".+_|[[:alpha:]]", "", fnames[i])
      }

    } else {
      df$mean[df$param %in% "Rsquared"] <- round(mean(BayesTraits_outs[[i]][["R.2"]]), digits=5)
      df$mean[df$param %in% paste0("mean", responvar)] <- round(10^(mean(BayesTraits_outs[[i]][["Alpha"]])), digits=2)
      if (gsub(".+_|[[:digit:]]|[.]|Estimate", "", fnames)[i] == "VR") {
        if (!is.null(treetransf)) {
          df$mean[df$param %in% treetransf] <- ifelse(df$param[df$param %in% treetransf] %in% "OU", "0", "1.0")
        }
      }
      if (gsub(".+_", "", fnames)[i] == "Fabric") {
        if (!is.null(treetransf)) {
          df$mean[df$param %in% treetransf] <- ifelse(df$param[df$param %in% treetransf] %in% "OU", "0", "1.0")
        }
      }
      if (gsub(".+_", "", fnames)[i] == "UNI") {
        if (!is.null(treetransf)) {
          df$mean[df$param %in% treetransf] <- ifelse(df$param[df$param %in% treetransf] %in% "OU", "0", "1.0")
        }
      }
    }
    mean_pMCMC[[i]] <- df
  }
  names(mean_pMCMC) <- fnames

  # Writing the values in the final dataframe to be output
  output_log_df <- tibble::tibble(Model = paste(responvar, "~", gsub("_", " + ", fnames)),
                                  Rsquared = NA,
                                  LL = as.numeric(paste(unlist(logL))),
                                  !!paste0("mean", responvar) := NA)

  if (!is.null(treetransf)) {
    for (i in seq_along(treetransf)) {
      output_log_df <- tibble::add_column(output_log_df, !!treetransf[i] := NA, .before = "Rsquared")
    }
  }

  output_log_df$Model <- greekLetters(output_log_df$Model)
  output_log_df$Model <- gsub("\\s[+]\\sUNI", "", output_log_df$Model)

  if (!is.null(treetransf)) {
    for (i in seq_along(explanvar)) {
      output_log_df <- tibble::add_column(output_log_df, !!explanvar[i] := NA, .before = treetransf[1])
    }
  } else {
    for (i in seq_along(explanvar)) {
      output_log_df <- tibble::add_column(output_log_df, !!explanvar[i] := NA, .before = "Rsquared")
    }
  }
  output_log_df <- as.data.frame(output_log_df)

  for (i in seq_along(mean_pMCMC)) {
    temp <- strsplit(gsub("_+[^_]+$", "", names(mean_pMCMC)[i]), "_")[[1]]
    if (!is.null(treetransf)) {
      output_log_df[i, c(temp, treetransf, "Rsquared", paste0("mean", responvar))] <- mean_pMCMC[[i]][["mean"]]
    } else {
      output_log_df[i, c(temp, "Rsquared", paste0("mean", responvar))] <- mean_pMCMC[[i]][["mean"]]
    }
    output_log_df[i, temp] <- paste(output_log_df[i, temp], paste0("(",
                                                                   na.omit(mean_pMCMC[[i]][["pMCMC"]]),
                                                                   ")"))
  }
  output_log_df <- output_log_df %>% arrange(desc(LL))
  output_log_df$LL <- as.character(output_log_df$LL)

  # Exclude models with uniform rates
  if (!unirates) {
    output_log_df <- output_log_df[grepl("\\S[+]VR[+]|VR$|Fabric$|κ|δ|λ|OU$|Estimate$", output_log_df$Model), ]
  }


  #_____________________________________________________________________________
  # Comparing models with Bayes Factor analyses
  if (bayesfactor) {
    bfactor_pairs <- data.frame(combinat::combn(fnames, 2))

    print("Calculating Bayes Factor between all models... ")
    bfactor_res <- data.frame(model_comparison = rep(NA, length(bfactor_pairs)),
                              BF = NA)
    for (i in seq_along(bfactor_pairs)) {
      logpairs <- unlist(logL[names(logL) %in% bfactor_pairs[[i]]])
      maxlog <- logpairs[logpairs %in% max(logpairs)]
      minlog <- logpairs[logpairs %in% min(logpairs)]
      if (length(maxlog) == 2) {
        maxlog <- maxlog[1]
        minlog <- minlog[2]
      }
      temp <- 2*((as.numeric(maxlog)) - (as.numeric(minlog)))
      bfactor_res$model_comparison[i] <- paste(names(maxlog), "vs.", names(minlog))
      bfactor_res$BF[i] <- round(temp, 5)
    }

    output_bfactor <- bfactor_res
    output_bfactor$model_comparison <- gsub("_UNI", "", output_bfactor$model_comparison)
    output_bfactor$model_comparison <- greekLetters(output_bfactor$model_comparison)
    output_bfactor$model_comparison <- gsub("_", "+", output_bfactor$model_comparison)

    n <- sort(unique(unlist(strsplit(bfactor_res$model_comparison, " vs. "))))
    output_bfactor_heat <- matrix(NA,
                                  nrow = length(n),
                                  ncol = length(n))
    row.names(output_bfactor_heat) <- n
    colnames(output_bfactor_heat) <- n

    for (i in seq_along(bfactor_res$model_comparison)) {
      temp <- sort(strsplit(bfactor_res$model_comparison[i], " vs. ")[[1]])
      r <- which(row.names(output_bfactor_heat) %in% temp[2] == TRUE)
      c <- which(colnames(output_bfactor_heat) %in% temp[1] == TRUE)
      output_bfactor_heat[r,c] <- bfactor_res$BF[i]
    }

    row.names(output_bfactor_heat) <- gsub("_UNI", "", row.names(output_bfactor_heat))
    row.names(output_bfactor_heat) <- gsub("_", "+", row.names(output_bfactor_heat))
    colnames(output_bfactor_heat) <- gsub("_UNI", "", colnames(output_bfactor_heat))
    colnames(output_bfactor_heat) <- gsub("_", "+", colnames(output_bfactor_heat))

    output_log_df$BF <- NA
    tf <- grepl("\\S[+]VR[+]|VR$|Fabric$|κ|δ|λ|OU$|Estimate$", output_log_df$Model)
    n_var <- gsub("DBH\\s~\\s|\\s", "", output_log_df$Model[tf])
    n_uni <- gsub("DBH\\s~\\s|\\s|[+]\\sVR|[+]\\sFabric|[+]\\sκ|[+]\\sδ|[+]\\sλ|[+]\\sOU", "", output_log_df$Model[tf])
    n_var <- paste(n_var, "vs.", n_uni)

    for (i in seq_along(n_var)) {
      f <- output_bfactor$model_comparison %in% n_var[i]
      if (length(which(f == FALSE)) == length(f)) {
        # For loaded Log files of the same model that do not included Uniform tree transformation
        output_log_df$BF[tf][i] <- paste(n_var[i], "were not compared!")
      } else {
        output_log_df$BF[tf][i] <- output_bfactor$BF[output_bfactor$model_comparison %in% n_var[i]]
      }
    }

  }


  #_____________________________________________________________________________
  # Create a new directory to save the results with current date
  # If there is no directory... make one!
  todaydate <- format(Sys.time(), "%d%b%Y")
  if (!dir.exists(paste0(dir_create, "/"))) {
    dir.create(paste0(dir_create, "/"))
  }
  if (!dir.exists(paste0(dir_create, "/", todaydate))) {
    dir.create(paste0(dir_create, "/", todaydate))
  }
  # If directory was created during a previous search, get its name to save results
  folder_name <- paste0(paste0(dir_create, "/"), todaydate)
  print(paste0("Writing '", folder_name, "' on disk."))

  #_____________________________________________________________________________
  # Save the Bayes Factor outputs
  if (bayesfactor) {
    print(paste0("Writing the matrix and figure of Bayes Factor pairwise comparisons on the disk folder ", folder_name))

    output_bfactor_heat_csv <- output_bfactor_heat

    colnames(output_bfactor_heat_csv) <- greekLetters(colnames(output_bfactor_heat_csv))
    row.names(output_bfactor_heat_csv) <- greekLetters(row.names(output_bfactor_heat_csv))

    write.csv(output_bfactor_heat_csv,
              file = paste0(folder_name, "/Bayes_Factor_phyreg_heatmap_matrix.csv"))

    cormat <- reshape2::melt(output_bfactor_heat)
    cormat$value <- round(cormat$value, 1)
    xlabs <- unique(gsub(paste0("[+]", "Kappa"), paste0("\n+", "Kappa"), cormat$Var1))
    xlabs <- unique(gsub(paste0("[+]", "Delta"), paste0("\n+", "Delta"), xlabs))
    xlabs <- unique(gsub(paste0("[+]", "Lambda"), paste0("\n+", "Lambda"), xlabs))
    xlabs <- unique(gsub(paste0("[+]", "OU"), paste0("\n+", "OU"), xlabs))
    xlabs <- gsub("[+]VR", "\n+VR", xlabs)
    xlabs <- gsub("[+]Fabric", "\n+Fabric", xlabs)
    ylabs <- unique(gsub(paste0("[+]", "Kappa"), paste0("\n+", "Kappa"), cormat$Var2))
    ylabs <- unique(gsub(paste0("[+]", "Delta"), paste0("\n+", "Delta"), ylabs))
    ylabs <- unique(gsub(paste0("[+]", "Lambda"), paste0("\n+", "Lambda"), ylabs))
    ylabs <- unique(gsub(paste0("[+]", "OU"), paste0("\n+", "OU"), ylabs))
    ylabs <- gsub("[+]VR", "\n+VR", ylabs)
    ylabs <- gsub("[+]Fabric", "\n+Fabric", ylabs)

    pal <- rev(wesanderson::wes_palette("Zissou1", 100, type = "continuous"))
    #scales::show_col(pal)


    #___________________________________________________________________________
    # Plotting a general heatmap with values inside
    pdf(paste0(folder_name, "/Bayes_Factor_phyreg_heatmap.pdf"), ...)
    p <- ggplot2::ggplot(data = cormat, aes(Var1, Var2, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = value), color = "white",
                size = ifelse(length(log_infiles)>10, 3, 6),
                na.rm = TRUE) +
      scale_x_discrete(labels = xlabs) +
      scale_y_discrete(labels = ylabs) +
      scale_fill_gradientn(colours = pal,
                           name="BF", na.value="gray80") +
      theme_minimal() +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
      theme(plot.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                       colour = "black", face = "bold"),
            text = element_text(size = ifelse(length(log_infiles)>10, 12, 16))) +
      theme(axis.text.y = element_text(colour = "black", face = "bold")) +
      xlab("") +
      ylab("")
    print(p)
    dev.off()


    #___________________________________________________________________________
    # Plotting the heatmap with dendrograms

    # Make the asymmetric matrix as symmetric so as to plot with the pheatmap package
    lower <- output_bfactor_heat
    upper <- output_bfactor_heat
    upper <- lower  #return two from one
    upper[lower.tri(upper)] <- t(upper)[lower.tri(t(upper))]
    lower <- t(lower)
    lower[lower.tri(lower)] <- t(lower)[lower.tri(lower)]
    cormatpheat <- lower

    pcor <- pheatmap::pheatmap(cormatpheat,
                               color = pal,
                               treeheight_col = ifelse(length(log_infiles)>10, 40, 80),
                               treeheight_row = ifelse(length(log_infiles)>10, 40, 80),
                               cutree_cols = 4,
                               cutree_rows = 4,
                               silent = TRUE)

    pcor$gtable$grobs[[1]]$gp=grid::gpar(lwd=1.0) # control for dendrogram
    pcor$gtable$grobs[[2]]$gp=grid::gpar(lwd=1.0) # control for dendrogram
    pcor$gtable$grobs[[3]]$gp=grid::gpar(lwd=0.5) # control for rects width
    pcor$gtable$grobs[[4]]$gp=grid::gpar(font="1") # control for site labels in columns
    pcor$gtable$grobs[[5]]$gp=grid::gpar(font="1") # control for site labels in columns

    pdf(paste0(folder_name, "/Bayes_Factor_phyreg_heatmap_cluster.pdf"), ...)
    print(pcor)
    dev.off()
  }


  #_____________________________________________________________________________
  # Save the output table in CSV format for the processed log files
  if (any(outformat == "CSV")) {
    print(paste0("Writing the table-formatted output ", paste0(outfile, ".csv"),
                 " file on the disk folder ", folder_name))
    write.csv(output_log_df,
              file = paste0(folder_name, "/", outfile, ".csv"))
  }

  #_____________________________________________________________________________
  # Save the output table in Word .docx format using rmarkdown
  if (any(outformat == "Word")) {
    title = "Processed BayesTraits phylogenetic regression outputs"
    cat(  "---\n", "title: ", "\"", title, "\"\n", file = "temp.Rmd", sep="")

    fpath <- system.file("extdata", "word-styles-reference-01.docx", package="InNOutBT")

    cat("output:\n", file = "temp.Rmd", append=T)
    cat("  rmarkdown::word_document:\n", file = "temp.Rmd", append=T)
    cat(paste0("    reference_docx: ", fpath, "\n"), file = "temp.Rmd", append=T)

    cat(
      "vignette: >\n",
      " %\\VignetteIndexEntry{Table skeleton}\n",
      " %\\VignetteEngine{knitr::rmarkdown}\n",
      " %\\VignetteEncoding{UTF-8}\n",
      "---\n",
      file = "temp.Rmd", sep="", append = T)
    cat("&nbsp;\n", fill=T, file = "temp.Rmd", sep="", append=T)

    # Adding the table legend
    cat(paste("__TABLE 1.__ ", tableleg), "\n", fill=T, file = "temp.Rmd", sep="", append=T)
    cat("&nbsp;\n", fill=T, file = "temp.Rmd", sep="", append=T)

    # Adding the resulting table
    cat(knitr::kable(output_log_df), fill=T, file = "temp.Rmd", sep="", append=T)

    suppressWarnings(rmarkdown::render("temp.Rmd",
                                       output_file = paste0(folder_name, "/", outfile, ".docx"), quiet=T))

    # Deleting the temp R markdown file
    unlink("temp.Rmd")
    cat("The Table skeleton of the processed BayesTraits phylogenetic regression outputs was saved in the disk folder:")
    cat("\n", folder_name)
  }

  # Return the resulting table as an object to the R environment
  return(output_log_df)
}

# Auxiliary function for writing Greek letters
greekLetters <- function(x) {
  x <- gsub("Kappa", "κ", x)
  x <- gsub("Delta", "δ", x)
  x <- gsub("Lmda|Lambda", "λ", x)
  return(x)
}

