#!/usr/bin/env Rscript

# Script for producing BayesTraits inputs from csv file of collection/occurrence data


#_______________________________________________________________________________
# Load packages

library(treeio)
#install.packages("devtools")
#devtools::install_github("DBOSlab/InNOutBT")
library(InNOutBT)


#_______________________________________________________________________________
# Producing BayesTraits input data for phylogenetic regression analyses

# See original function (file BayesTraits_inputs.R) for a full description of 
# all its available arguments


# Load data
traitdata <- read.csv("Data/vataireoids_1610_25May2022_BayesTraits.csv")
phylo <- read.beast("Data/vatcomb30_24May22_Yule.con.tree")


#_______________________________________________________________________________
# Create sample of trait (distribution) data using all explanatory variables
# The first trait listed serves as the response variable
BayesTraits.inputs(tree = phylo,
                   data = traitdata,
                   tipscol = "terminal",
                   NodeCount = TRUE,
                   logtransf = c("DBH", "NodeCount"),
                   traitcols = c("bio12", "bio15"),
                   addtraits = c("log10NodeCount", "log10DBH"),
                   ordtraits = c("log10DBH", "bio12", "bio15", "log10NodeCount"),
                   dir_create = "results_BayesTraits_input",
                   fileDistData = "BayesTraits_linked_data_bio12_bio15_nnodes.txt",
                   fileMeanData = "BayesTraits_mean_data_bio12_bio15_nnodes.txt",
                   fileOrigData = "vataireoids_1610_25May2022_BayesTraits_netnodes_logtransf.csv")
                   
                   
#_______________________________________________________________________________
# Create sample of trait data using a subset of explanatory variables
BayesTraits.inputs(tree = phylo,
                   data = traitdata,
                   tipscol = "terminal",
                   NodeCount = FALSE,
                   logtransf = c("DBH"),
                   traitcols = c("bio12", "bio15"),
                   addtraits = c("log10DBH"),
                   ordtraits = c("log10DBH", "bio12", "bio15"),
                   dir_create = "results_BayesTraits_input",
                   fileDistData = "BayesTraits_linked_data_bio12_bio15.txt",
                   fileMeanData = "BayesTraits_mean_data_bio12_bio15.txt")
