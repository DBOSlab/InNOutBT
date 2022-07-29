#!/bin/zsh


## Give the path to both  the project run folder where the R scripts are, and the BayesTraits
## folder where the executable program is installed.
PROJECTFOLDER="/Users/domingoscardoso/BayesTraits_shell_run_project"
BAYESTRFOLDER="/Users/domingoscardoso/BayesTraitsV4"


## ---------------------------------------------------------------------------------------
## DO NOT CHANGE ANYTHING BELOW


## Producing BayesTraits inputs
./1_producing_BayesTraits_inputs.R


## Producing BayesTraits shell script
./2_producing_BayesTraits_shell.R


## Sending an email message
echo "Your BayesTraits input and shell script files are done. Now running BayesTraits program..." | mail -s "BayesTraits shell run" "cardosobot@gmail.com"


## Running BayesTraits program with the automatically produced shell script
cd $BAYESTRFOLDER
chmod 777 run_BayesTraits_shell_unix.sh
./run_BayesTraits_shell_unix.sh


## Moving the resulting outputs into the main project folder
mv BayesTraits_outputs_log_stone $PROJECTFOLDER
mv BayesTraits_outputs_sch_trees $PROJECTFOLDER


## Processing BayesTraits output
cd $PROJECTFOLDER
./3_processing_BayesTraits_outputs.R


## Sending an email message when the analysis is fully finished
echo "Your BayesTraits run and the InNOutBT processing of associated BayesTraits output files are fully done." | mail -s "BayesTraits shell run" "cardosobot@gmail.com"


## Opening automatically the resulting .pdf and .docx files
CURRENTDATE=$(date +%d%b%Y)
open results_BayesTraits_output/"$CURRENTDATE"/Bayes_Factor_heatmap_cluster.pdf results_BayesTraits_output/"$CURRENTDATE"/Bayes_Factor_heatmap.pdf results_BayesTraits_output/"$CURRENTDATE"/BayesTraits_output_table.docx
