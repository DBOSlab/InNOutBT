
<!-- README.md is generated from README.Rmd. Please edit that file -->

# InNOutBT

<!-- badges: start -->
<!-- badges: end -->

An R package for creating input data to and processing output data from
Meade & Pagel’s (2022)
[BayesTraits](http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/BayesTraitsV4.0.0.html)
computer package.

## Overview

The package InNOutBT is still under development. Currently it has only
three major functions that help to automatize BayesTraits phylogenetic
independent contrast regression analyses. To better understand how these
functions work, see [BayesTraits V4.0.0
Manual](http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/Files/BayesTraitsV4.0.0-Manual.pdf)
for further details on the specific format of the input files that are
required by BayesTraits, as well as the resulting output files after
analyses.

## Installation

You can install the development version of InNOutBT from
[GitHub](https://github.com/) with:

``` r
#install.packages("devtools")
devtools::install_github("DBOSlab/InNOutBT")
```

``` r
library(InNOutBT)
```

## Usage

A general description of each of the three available functions
(`BayesTraits.inputs`, `BayesTraits.shell`, and `BayesTraits.outputs`)
are provided below.

#### *1. `BayesTraits.inputs`: Creating input linked and mean data for phylogenetic regression analyses*

This function produces all BayesTraits input data of linked samples of
trait/climate data (distribution data), by automatically generating two
.txt files: one including the linked traits and the other with mean
values for each included trait. The function also includes an argument
to calculate the net node count for each species/accession given a
phylogenetic tree of reference, which has been successfully used to
explore the effect of speciation rates (e.g. [O’Donovan et
al. 2018](https://doi.org/10.1038/s41559-017-0454-6)) on the response
variable.

##### Example with `BayesTraits.inputs`:

``` r
library(InNOutBT)
library(treeio)

traitdata <- read.csv("Data/vataireoids_1610_25May2022_BayesTraits.csv")
phylo <- read.beast("Data/vatcomb30_24May22_Yule.con.tree")

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
```

#### *2. `BayesTraits.shell`: Creating shell script for running phylogenetic regression analyses*

Based on the directory paths to the tree, linked and mean data input
files, this function automatically produces shell scripts in both Unix
and Windows PowerShell formats for running one or multiple phylogenetic
regression analyses. In addition to generate the shell script directly
in BayesTraits directory, the function also creates all running files
with BayesTraits specified parameters (bi, it, sa, and st) for any run,
depending on the chosen model of tree transformation.

##### Example with `BayesTraits.shell`:

``` r
library(InNOutBT)

BayesTraits.shell(meanfile_dir = "results_BayesTraits_input/03Jun2022",
                  linkfile_dir = "results_BayesTraits_input/03Jun2022",
                  treefile_dir = "Data",
                  BayesTraits_dir = "/Users/domingoscardoso/BayesTraitsV4",
                  responvar = "DBH",
                  treetransf = c("UNI", "Lambda", "VR"),
                  bi = "30000000",
                  it = "130000000",
                  sa = "100000",
                  st = c("500", "50000"),
                  syst = "unix",
                  dir_create = "results_BayesTraits_shell",
                  cc_DataTree = TRUE)
```

#### *3. `BayesTraits.outputs`: Processing output data from phylogenetic regression analyses*

This function processes the output Log and Stones .txt files from
phylogenetic regression analyses. The resulting processed files are
tables either in Word or CSV formats summarizing each model, including
comparison with Bayes Factor. Heatmaps based on pairwise Bayes Factor
comparison are also automatically produced in PDF format.

##### Example with `BayesTraits.outputs`:

``` r
library(InNOutBT)

BayesTraits.outputs(logst_dir = paste0(getwd(),"/BayesTraits_outputs_log_stone"),
                    responvar = "DBH",
                    explanvar = c("bio12", "bio15", "nnodes", "nnodes2"),
                    treetransf = c("Kappa", "Delta", "Lambda", "OU"),
                    bayesfactor = TRUE,
                    unirates = FALSE,
                    outformat = c("Word", "CSV"),
                    tableleg = "Phylogenetic regression in BayesTraits: models and coefficients.",
                    dir_create = "results_BayesTraits_output",
                    outfile = "BayesTraits_output_table",
                    height = 17,
                    width = 20)
```

## Documentation

A detailed description of the **InNOutBT**’s full functionality is
available in [here](https://dboslab.github.io/InNOutBT/).

## Citation

Cardoso, D. & Lavin, M. (2022). InNOutBT: An R package for creating
input data to and processing output data from BayesTraits analyses.
<https://github.com/dboslab/InNOutBT>

<img src="man/figures/DBOSlab_logo.png" style="width:30.0%" />
