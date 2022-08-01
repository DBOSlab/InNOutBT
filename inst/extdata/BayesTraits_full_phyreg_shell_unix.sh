#!/bin/zsh


## Give the path to both  the project run folder where the R scripts are, and the BayesTraits
## folder where the executable program is installed.
PROJECTFOLDER="/Users/domingoscardoso/BayesTraits_shell_phyreg_project"
BAYESTRFOLDER="/Users/domingoscardoso/BayesTraitsV4"

EMAIL="cardosobot@gmail.com"


## If not previously done, run the following command to set permission to run this shell script
## chmod 777 BayesTraits_full_phyreg_shell_unix.sh

## To execute this shell script, run the following command
## ./BayesTraits_full_phyreg_shell_unix.sh


## ---------------------------------------------------------------------------------------
## DO NOT CHANGE ANYTHING BELOW


## Producing BayesTraits inputs
./1_producing_BayesTraits_phyreg_inputs.R


## Producing BayesTraits shell script
./2_producing_BayesTraits_phyreg_shell.R


## Sending an email message
echo "Your BayesTraits input and shell script files are done. Now running BayesTraits program..." | mail -s "BayesTraits shell run" $EMAIL


## Running BayesTraits program with the automatically produced shell script
cd $BAYESTRFOLDER
chmod 777 run_BayesTraits_phyreg_shell_unix.sh
./run_BayesTraits_phyreg_shell_unix.sh


## Moving the resulting outputs into the main project folder
mv BayesTraits_phyreg_outputs_log_stone $PROJECTFOLDER
mv BayesTraits_phyreg_outputs_sch_trees $PROJECTFOLDER


## Processing BayesTraits output
cd $PROJECTFOLDER
./3_processing_BayesTraits_phyreg_outputs.R


## Sending an email message when the analysis is fully finished
echo "Your BayesTraits run and the InNOutBT processing of associated BayesTraits output files are fully done." | mail -s "BayesTraits shell run" $EMAIL


## Opening automatically the resulting .pdf and .docx files
CURRENTDATE=$(date +%d%b%Y)
open results_BayesTraits_phyreg_output/"$CURRENTDATE"/Bayes_Factor_phyreg_heatmap_cluster.pdf results_BayesTraits_phyreg_output/"$CURRENTDATE"/Bayes_Factor_phyreg_heatmap.pdf results_BayesTraits_phyreg_output/"$CURRENTDATE"/BayesTraits_phyreg_output_table.docx
