# Create a BayesTraits project folder with associated example scripts and input data for a full shell script run and log processing of BayesTraits phylogenetic regression analyses

By reporting a directory path and the computer's operating system
("unix" or "windows"), this function creates a BayesTraits project
folder with associated example scripts and input data for a full shell
script run and log processing of multiple phylogenetic regression
analyses with Meade & Pagel's (2022)
[BayesTraits](https://www.evolution.reading.ac.uk/BayesTraitsV5.0.3/BayesTraitsV5.0.3.html)
program. The function produces an example of a shell script in Unix or
Windows PowerShell that helps to fully automate both the complete
BayesTraits regression analysis and the processing of the resulting log
and stones files, by automatically running the three main `InNOutBT`
functions
[`phyreg.inputs`](https://DBOSlab.github.io/InNOutBT/reference/phyreg.inputs.md),
[`phyreg.shell`](https://DBOSlab.github.io/InNOutBT/reference/phyreg.shell.md),
and
[`phyreg.outputs`](https://DBOSlab.github.io/InNOutBT/reference/phyreg.outputs.md),
without directly opening the R or RStudio programs. Note that these
files are just examples from a specific study of the Vataireoid legumes,
and so you will have to edit the folder paths and analysis parameters in
these example scripts according to your desired regression analyses. See
also a more complete article [describing a BayesTraits shell run
project](https://dboslab.github.io/InNOutBT/articles/phyreg_shell_project.html)
that can be created with `shell.phyreg.proj` for fully automated
multiple phylogenetic regression analyses.

## Usage

``` r
shell.phyreg.proj(dir_create_proj = NULL,
                  syst = NULL)
```

## Arguments

- dir_create_proj:

  Path to any directory on your computer where an example of a
  BayesTraits shell project folder and associated files for a full run
  with shell script should be saved.

- syst:

  Report the operating system ("unix" or "windows").

## Details

After creating the project folder with this function, follow these steps
to run the project. These steps are for macOS (`syst = "unix"`); see the
article linked above for the steps in Windows PowerShell.

**Before the first run**

1.  Install the BayesTraits program in any folder and note the path to
    that folder and the name of the executable file (e.g.,
    `BayesTraitsV5`).

2.  Check that the executable file starts by running its full path in
    the Terminal. If macOS stops it (e.g., the message `killed`), remove
    the block that macOS adds to downloaded files, provided that you
    trust the source of the download:
    `xattr -d com.apple.quarantine /path/to/BayesTraits/BayesTraitsV5`.

3.  Check that the R packages `InNOutBT` and `treeio` are installed, and
    that `pandoc` is available in the Terminal, as the Word table needs
    it: `Rscript -e 'rmarkdown::pandoc_available()'`.

**Edit the project files in any text editor**

1.  In `BayesTraits_full_phyreg_shell_unix.sh`, set `PROJECTFOLDER` to
    the path of the project folder, `BAYESTRFOLDER` to the path of the
    folder that contains the BayesTraits executable file, and `EMAIL` to
    your e-mail address.

2.  In `2_producing_BayesTraits_phyreg_shell.R`, set `BayesTraits_dir`
    to the same path as `BAYESTRFOLDER`, and `BayesTraits_exe` to the
    name of the executable file.

3.  To run your own study, also edit the three numbered R scripts and
    replace the files in the `Data` folder.

**Run the project in the Terminal**

1.  Move to the project folder, not to the BayesTraits folder, because
    the main shell script moves to the BayesTraits folder by itself:
    `cd /path/to/BayesTraits_shell_phyreg_project`.

2.  Allow the main shell script to run, which is needed only once:
    `chmod 777 BayesTraits_full_phyreg_shell_unix.sh`.

3.  Start the full analysis: `./BayesTraits_full_phyreg_shell_unix.sh`.
    Keep the Terminal open until it finishes, which can take a long time
    with the default number of iterations.

The script then produces the input files, creates and runs the
BayesTraits analyses, moves the resulting log and stones files to the
project folder, and processes them. The final table and heatmaps are
saved in the project folder, in `results_BayesTraits_phyreg_output`,
within a subfolder named with the current date.

## See also

[`phyreg.inputs`](https://DBOSlab.github.io/InNOutBT/reference/phyreg.inputs.md)

[`phyreg.shell`](https://DBOSlab.github.io/InNOutBT/reference/phyreg.shell.md)

[`phyreg.outputs`](https://DBOSlab.github.io/InNOutBT/reference/phyreg.outputs.md)

## Author

Domingos Cardoso & Matt Lavin

## Examples

``` r
if (FALSE) { # \dontrun{
library(InNOutBT)

shell.phyreg.proj(dir_create_proj = "/Users/domingoscardoso",
                  syst = "unix")
} # }
```
