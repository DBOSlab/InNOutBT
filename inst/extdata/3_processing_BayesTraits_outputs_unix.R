#!/usr/bin/env Rscript

# Script for processing BayesTraits outputs


#_______________________________________________________________________________
# Load packages

library(InNOutBT)


#_______________________________________________________________________________
# Processing all BayesTraits output data from phylogenetic regression analyses

BayesTraits.outputs(logst_dir = paste0(getwd(),"/BayesTraits_outputs_log_stone"),
                    responvar = "DBH",
                    explanvar = c("bio12", "bio15", "nnodes"),
                    treetransf = c("UNI", "Lambda", "VR"),
                    bayesfactor = TRUE,
                    unirates = FALSE,
                    outformat = c("Word", "CSV"),
                    tableleg = "Phylogenetic regression in BayesTraits: models and coefficients.",
                    dir_create = "results_BayesTraits_output",
                    outfile = "BayesTraits_output_table",
                    height = 12,
                    width = 15)