#!/usr/bin/env Rscript

# Script for producing BayesTraits shell script for phylogenetic regression analyses


#_______________________________________________________________________________
# Load packages

library(InNOutBT)


#_______________________________________________________________________________
# Producing BayesTraits shell script and associated run commands

todaydate <- format(Sys.time(), "%d%b%Y")

BayesTraits.shell(meanfile_dir = paste0("results_BayesTraits_input/", todaydate),
                  linkfile_dir = paste0("results_BayesTraits_input/", todaydate),
                  treefile_dir = "Data",
                  BayesTraits_dir = "/Users/domingoscardoso/BayesTraitsV4",
                  responvar = "DBH",
                  treetransf = c("UNI", "Lambda", "VR"),
                  Kappa = NULL,
                  Delta = NULL,
                  Lambda = c("0.001", "1.0", "Estimate"),
                  OU = NULL,
                  VRLambda = TRUE,
                  bi = "300000", #"30000000"
                  it = "1300000", #"130000000"
                  sa = "1000", #"100000"
                  st = c("50", "500"), #c("500", "50000")
                  syst = "unix",
                  dir_create = "results_BayesTraits_shell",
                  cc_DataTree = TRUE)
