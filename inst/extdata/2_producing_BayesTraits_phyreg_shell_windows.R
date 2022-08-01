# Script for producing BayesTraits shell script for phylogenetic regression analyses
# Version for a full run using Windows PowerShell script

#_______________________________________________________________________________
# Load packages

library(InNOutBT)


#_______________________________________________________________________________
# Producing BayesTraits shell script and associated run commands

todaydate <- format(Sys.time(), "%d%b%Y")

phyreg.shell(meanfile_dir = paste0("results_BayesTraits_phyreg_input/", todaydate),
             linkfile_dir = paste0("results_BayesTraits_phyreg_input/", todaydate),
             treefile_dir = "Data",
             BayesTraits_dir = "C:/Users/Matt_/Documents/BayesTraitsV4.0.0",
             responvar = "DBH",
             treetransf = c("Lambda", "VR"),
             Kappa = NULL,
             Delta = NULL,
             Lambda = c("0.001", "1.0", "Estimate"),
             OU = NULL,
             VRLambda = TRUE,
             bi = "30000000",
             it = "130000000",
             sa = "100000",
             st = c("500", "5000"),
             syst = "windows",
             dir_create = "results_BayesTraits_phyreg_shell",
             cc_DataTree = TRUE)
