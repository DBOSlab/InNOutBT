# Create shell scripts to run BayesTraits phylogenetic regression analyses using samples of trait data

By reporting the directory paths to the tree and data files, this
function produces shell scripts in both Unix and Windows PowerShell
formats for running one or multiple phylogenetic regression analyses
with Meade & Pagel's (2022)
[BayesTraits](https://www.evolution.reading.ac.uk/BayesTraitsV5.0.3/BayesTraitsV5.0.3.html)
program. In addition to generating the shell script directly in the
BayesTraits directory, the function also creates BayesTraits batch files
with specified parameters, including generations set for burnin,
iterations, and the stepping stones sampler (e.g., bi, it, sa, and st).
In addition, models involving alternative tree transformations can be
specified. See the [BayesTraits 5.0.3
Manual](https://www.evolution.reading.ac.uk/BayesTraitsV5.0.3/Files/BayesTraitsV5.0.2-Manual.pdf)
for further details involving BayesTraits parameters and default
settings for independent contrast regression analysis so as to better
understand how to set the arguments of this function `phyreg.shell`. You
might consider using the function
[`phyreg.inputs`](https://DBOSlab.github.io/InNOutBT/reference/phyreg.inputs.md)
to generate properly formatted input files of mean and linked samples of
trait data prior to using `phyreg.shell`. See also a more complete
article on how to create a [shell
script](https://dboslab.github.io/InNOutBT/articles/phyreg_shell.html)
with `phyreg.shell` for multiple phylogenetic regression analyses.

## Usage

``` r
phyreg.shell(meanfile_dir = NULL,
             linkfile_dir = NULL,
             treefile_dir = NULL,
             BayesTraits_dir = NULL,
             BayesTraits_exe = "BayesTraitsV5",
             responvar = NULL,
             treetransf = NULL,
             Kappa = NULL,
             Delta = NULL,
             Lambda = NULL,
             OU = NULL,
             VRLambda = FALSE,
             bi = "30000000",
             it = "130000000",
             sa = "100000",
             st = c("500", "50000"),
             syst = NULL,
             dir_create = "results_BayesTraits_phyreg_shell",
             cc_DataTree = TRUE)
```

## Arguments

- meanfile_dir:

  Path to the folder directory where the data file(s) containing trait
  means for samples of trait data are stored.

- linkfile_dir:

  Path to the folder directory where the data file(s) containing samples
  of linked trait data are stored.

- treefile_dir:

  Path to the folder directory where the tree file is stored.

- BayesTraits_dir:

  Path to the folder directory where the BayesTraits executable file is
  stored.

- BayesTraits_exe:

  Name of the BayesTraits executable file stored in `BayesTraits_dir`,
  which depends on the installed version of the program. The default is
  "BayesTraitsV5". For example, set "BayesTraitsV4" if you are running
  BayesTraits V4.

- responvar:

  Report the name of the response variable of your analysis.

- treetransf:

  In vector format, report tree transformation models to be estimated,
  such as *Kappa*, *Delta*, *Lambda*, or *OU* (Ornstein–Uhlenbeck). A
  model of variable rates *VR* or the *Fabric* model can be reported
  here, the latter of which is a statistical model accommodating an
  uneven evolutionary landscape, as described by [Pagel et al.
  (2022)](https://doi.org/10.1038/s41467-022-28595-z). We use *UNI* as a
  designation for uniform rates of evolution where *VR* is not estimated
  and no other tree transformation models are estimated (e.g., Lambda is
  assumed to be 1.0). In this case, "UNI" sets a complete default
  BayesTraits regression analysis.

- Kappa:

  If left as the default `NULL`, this tree transformation model will be
  estimated, provided that it was also defined in the previous parameter
  `treetransf`. But here you can also set a vector with any specific
  value between 0 and 1. If you want to run analyses with both specific
  and estimated values for *Kappa*, then you have to set a vector here,
  such as c("0.001", "1.0", "Estimate"). Report any number as a
  character (i.e., in quotes).

- Delta:

  If left as the default `NULL`, this tree transformation model will be
  estimated, provided that it was also defined in the previous parameter
  `treetransf`. But here you can also set a vector with any specific
  value between 0 and 1. If you want to run analyses with both specific
  and estimated values for *Delta*, then you have to set a vector here,
  such as c("0.001", "1.0", "Estimate"). Report any number as a
  character (i.e., in quotes).

- Lambda:

  If left as the default `NULL`, this tree transformation model will be
  estimated, provided that it was also defined in the previous parameter
  `treetransf`. But here you can also set a vector with any specific
  value between 0 and 1. If you want to run analyses with both specific
  and estimated values for *Lambda*, then you have to set a vector here,
  such as c("0.001", "1.0", "Estimate"). Report any number as a
  character (i.e., in quotes).

- OU:

  If left as the default `NULL`, this tree transformation model will be
  estimated, provided that it was also defined in the previous parameter
  `treetransf`. But here you can also set a vector with any specific
  value between 0 and 1. If you want to run analyses with both specific
  and estimated values for *OU*, then you have to set a vector here,
  such as c("0.001", "1.0", "Estimate"). Report any number as a
  character (i.e., in quotes).

- VRLambda:

  Logical, the default is `FALSE`. If set as `TRUE`, then the specific
  values set for the tree transformation models *Kappa*, *Delta*,
  *Lambda* or *OU* will be considered as separate runs when using the
  model of variable rates *VR*, provided that you have also chosen this
  model in the parameter `treetransf`.

- bi:

  Set the number of iterations for the MCMC burnin. Report this number
  as a character (i.e., in quotes).

- it:

  Set the total number of iterations for the MCMC chain. Report this
  number as a character (i.e., in quotes).

- sa:

  Set the sample frequency. Report this number as a character (i.e., in
  quotes).

- st:

  The stepping stone sampler to estimate the marginal likelihood. Report
  the number of stones to use and the number of iterations for each
  stone as characters (i.e., a vector of numbers each in quotes).

- syst:

  Report the operating system ("unix" or "windows").

- dir_create:

  Path to the computer's working directory or to the directory where the
  file(s) will be saved. The default setting creates a directory named
  "results_BayesTraits_phyreg_shell", in which the results will be saved
  within a subfolder named with the current date.

- cc_DataTree:

  Logical, the default is `TRUE`. A copy of the tree file and the mean
  and linked data files will be made in the BayesTraits folder where the
  BayesTraits executable file is located.

## See also

[`phyreg.inputs`](https://DBOSlab.github.io/InNOutBT/reference/phyreg.inputs.md)

[`phyreg.outputs`](https://DBOSlab.github.io/InNOutBT/reference/phyreg.outputs.md)

## Author

Domingos Cardoso & Matt Lavin
