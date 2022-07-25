
<!-- README.md is generated from README.Rmd. Please edit that file -->

# InNOutBT <img src="man/figures/innoutbt.png" align="right" alt="" width="120" />

<!-- badges: start -->
<!-- badges: end -->

An R package for creating input data and processing output data in
formats required of Meade & Pagel’s (2022)
[BayesTraits](http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/BayesTraitsV4.0.0.html)
computer package.  
  

## Overview

The package **InNOutBT** remains under development. Currently, it has
three major functions that help to automate BayesTraits phylogenetic
independent contrast regression analyses. To better understand how these
functions work, see [BayesTraits V4.0.0
Manual](http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/Files/BayesTraitsV4.0.0-Manual.pdf)
for further details on the specific format of the input files that are
required by BayesTraits, as well as the BayesTraits strucuture of
resulting output files.  
  

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
  

#### *1. `BayesTraits.inputs`: Creating input mean and linked data for phylogenetic regression analyses*

This function produces the required BayesTraits input data of linked
samples of trait/climate data (distribution data), by generating two
text files, one of which includes the mean values associated with each
trait associated with a terminal taxon and a second of which reports
samples of trait values that are linked to each other because they were
sampled from one accession of the associated terminal taxon. The
function also includes an argument to calculate the net node count for
each species/accession given a phylogenetic tree. Net node counts have
been used to explore the effects of speciation rates on the response
variable (e.g. [O’Donovan et
al. 2018](https://doi.org/10.1038/s41559-017-0454-6)).  
  

##### Example of using `BayesTraits.inputs`:

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

Using directory paths to the tree and the input files of mean and linked
trait data, this function produces shell scripts in both Unix and
Windows PowerShell formats for running one or multiple phylogenetic
regression analyses. This function creates the necessary batch file
commands, including the desired settings for the burnin, iterations, and
stepping stones parameters (e.g., abbreviated as bi, it, sa, and st) for
a given run, in addition to incorporating any model of tree
transformation. In the example below, we use “UNI”, “Lambda”, and “VR”
in the treefransf argument to indicate we want three separate sets of
batch file commands all identical except one, designated as “UNI” will
not estimate “Lambda” or “VR” (variable rates). We use “UNI” to indicate
the assumption of uniform rates of evolution and no other tree
transformation (e.g., such as estimating “Lambda”). In this case,
“Lambda” assumes the default setting of 1.0 (see the BayesTraits
manual).  
  

##### Example of using `BayesTraits.shell`:

``` r
library(InNOutBT)

BayesTraits.shell(meanfile_dir = "results_BayesTraits_input/03Jun2022",
                  linkfile_dir = "results_BayesTraits_input/03Jun2022",
                  treefile_dir = "Data",
                  BayesTraits_dir = "/Users/domingoscardoso/BayesTraitsV4",
                  responvar = "DBH",
                  treetransf = c("UNI", "Kappa", "Lambda", "VR"),
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
comparisons of Bayes Factors. Heatmaps using pairwise Bayes Factor
comparisons are also produced in PDF format (we use heatmaps only to
identify models similar to the model with the highest likelihood
score).  
  

##### Example of using `BayesTraits.outputs`:

``` r
library(InNOutBT)

BayesTraits.outputs(logst_dir = paste0(getwd(),"/BayesTraits_outputs_log_stone"),
                    responvar = "DBH",
                    explanvar = c("bio12", "bio15", "nnodes"),
                    treetransf = c("UNI", "Kappa", "Lambda", "VR"),
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
available in different
[articles](https://dboslab.github.io/InNOutBT/).  
  

## Citation

Cardoso, D. & Lavin, M. (2022). InNOutBT: An R package for creating
input data and processing output data required for BayesTraits analyses.
<https://github.com/dboslab/InNOutBT>

<img src="man/figures/DBOSlab_logo.png" align="left" alt="" width="120" />
