# Process BayesTraits phylogenetic regression output data

This function processes the output reported in the log and stones .txt
files, which are generated from the phylogenetic regression analyses in
Meade & Pagel's (2022)
[BayesTraits](https://www.evolution.reading.ac.uk/BayesTraitsV5.0.3/BayesTraitsV5.0.3.html)
program. The resulting files are tables in either Word or CSV format,
and these rank competing models by their likelihood scores. Bayes
Factors can also be reported in the case where models with tree
transformations (e.g., either Lambda or variable rates are estimated)
always have higher likelihood scores than simpler models (e.g., those
with neither Lambda nor variable rates estimated). In this case, the
reported Bayes Factors pertain to the more complex model listed in the
table and the simpler nested model that did not involve a tree
transformation. Heatmaps based on pairwise Bayes Factor comparisons are
also automatically produced in PDF format. These are intended to
identify models that are most similar to the model with the highest
likelihood score. See the [BayesTraits 5.0.3
Manual](https://www.evolution.reading.ac.uk/BayesTraitsV5.0.3/Files/BayesTraitsV5.0.2-Manual.pdf)
for details of the output files produced by BayesTraits, which may
facilitate an understanding of how output files are processed with the
function `phyreg.outputs`. See also a more complete article on how to
use `phyreg.outputs` to [process the resulting output
files](https://dboslab.github.io/InNOutBT/articles/phyreg_outputs.html)
from the phylogenetic regression analyses.

## Usage

``` r
phyreg.outputs(logst_dir = NULL,
               responvar = NULL,
               explanvar = NULL,
               treetransf = c("Kappa", "Delta", "Lambda", "OU", "UNI", "VR", "Fabric"),
               bayesfactor = FALSE,
               unirates = TRUE,
               outformat = c("Word", "CSV"),
               tableleg = "Phylogenetic regression in BayesTraits",
               dir_create = "results_BayesTraits_phyreg_output",
               outfile = "BayesTraits_phyreg_output_table",
               ...)
```

## Arguments

- logst_dir:

  Path to the directory where the log and stepping stones files
  generated during the independent contrast regression analysis are
  stored. The log file(s) contain the model options and output. The
  stepping stones file(s) contain the marginal log likelihood for a
  given model. These files are labeled with the explanatory variable(s)
  of the relevant model, where each variable is separated by an
  underscore. Explanatory variables include any tree transformation
  included in a model, and these can include *Kappa*, *Delta*, *Lambda*,
  *OU* (Ornstein–Uhlenbeck), *UNI* (for no tree transformation), *VR*
  (variable rates) or *Fabric* (for a model that accommodates an uneven
  evolutionary landscape, as described by [Pagel et al.
  (2022)](https://doi.org/10.1038/s41467-022-28595-z)). Thus, log
  file(s) have a format like
  "BayesTraits_mean_data_bio12_bio15_VR.log.txt" and
  "BayesTraits_mean_data_bio12_bio15_nnodes_Lambda.log.txt" (where
  nnodes signifies net nodes). This indicates that the first log file
  contains three explanatory variables ("bio12", "bio15", and variable
  rates), and the second log file contains four explanatory variables
  ("bio12", "bio15", "nnodes", and Lambda). Stepping stones file(s)
  should be named correspondingly (e.g.,
  "BayesTraits_mean_data_bio12_bio15_VR.stones.txt" and
  "BayesTraits_mean_data_bio12_bio15_nnodes_Lambda.stones.txt").

- responvar:

  Report the name of the response variable.

- explanvar:

  Report the name(s) of the explanatory variable(s) exactly as they are
  written in the names of the log and stones files being processed.

- treetransf:

  Select any of the available BayesTraits tree transformations (*Kappa*,
  *Delta*, *Lambda*, *OU*, *UNI*, *VR* or *Fabric*) that are included
  among the competing models.

- bayesfactor:

  Logical, if `FALSE`, then no Bayes Factor pairwise comparisons will be
  performed among the models.

- unirates:

  Logical, if `FALSE`, then models lacking tree transformations
  (*Kappa*, *Delta*, *Lambda* or *OU*), variable rates, and Fabric
  (i.e., default values and uniform rates of evolution are assumed) will
  not be displayed in the main table.

- outformat:

  Define either "Word" or "CSV", or both as a vector c("Word", "CSV"),
  for writing the results in those formats.

- tableleg:

  A legend for the main table, provided that you have also chosen
  `outformat = "Word"`. If you do not report a legend, then the
  Word-formatted table will have the default legend **Phylogenetic
  regression in BayesTraits**.

- dir_create:

  Path to the directory where the file(s) will be saved. The default
  setting creates a directory named
  **results_BayesTraits_phyreg_output** where the results will be saved
  in a subfolder named with the current date.

- outfile:

  Name of the resulting table-formatted files in either "Word", "CSV",
  or both (depending on the specified argument `outformat`). If no name
  is reported, the default setting creates a file named
  *BayesTraits_phyreg_output_table.docx* and/or
  *BayesTraits_phyreg_output_table.csv*.

- ...:

  Additional parameters passed to pdf.

## Value

Table in data frame format, which is also saved in .csv format.

## See also

[`phyreg.inputs`](https://DBOSlab.github.io/InNOutBT/reference/phyreg.inputs.md)

## Author

Domingos Cardoso & Matt Lavin
