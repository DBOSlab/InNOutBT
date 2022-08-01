# Script for processing BayesTraits outputs
# Version for a full run using Windows PowerShell script

#_______________________________________________________________________________
# Load packages

library(InNOutBT)


#_______________________________________________________________________________
# Processing all BayesTraits output data from phylogenetic regression analyses

phyreg.outputs(logst_dir = paste0(getwd(),"/BayesTraits_phyreg_outputs_log_stone"),
               responvar = "DBH",
               explanvar = c("bio12", "bio15", "nnodes"),
               treetransf = c("Lambda", "VR"),
               bayesfactor = TRUE,
               unirates = FALSE,
               outformat = c("Word", "CSV"),
               tableleg = "Phylogenetic regression in BayesTraits: models and coefficients.",
               dir_create = "results_BayesTraits_phyreg_output",
               outfile = "BayesTraits_phyreg_output_table",
               height = 12,
               width = 15)

