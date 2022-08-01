#!/bin/bash


## Give the paths to both the project run folder where the R scripts are, 
## and the BayesTraits folder where the executable program is installed.
$PROJECTFOLDER="C:\Users\Matt_\Documents\BayesTraits_shell_phyreg_project"
$BAYESTRFOLDER="C:\Users\Matt_\Documents\BayesTraitsV4.0.0"

## Give the paths to the Rscript, Adobe and Microsoft Word executable files 
$RFOLDER="C:\'Program Files'\R\R-4.2.0\bin\Rscript"
$ACROBATFOLDER="C:\'Program Files (x86)'\Adobe\'Acrobat DC'\Acrobat\Acrobat.exe"
$WORDFOLDER="C:\'Program Files (x86)'\'Microsoft Office'\root\Office16\WINWORD.EXE"


## If not previously done, run the following commands to set permission to run this PowerShell script
## Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope LocalMachine
## Unblock-File -Path .\BayesTraits_full_phyreg_shell_windows.ps1

## To execute this shell script, run the following command
## .\BayesTraits_full_phyreg_shell_windows.ps1


## ---------------------------------------------------------------------------------------
## DO NOT CHANGE ANYTHING BELOW


## Producing BayesTraits inputs
Invoke-Expression "$RFOLDER .\1_producing_BayesTraits_phyreg_inputs.R"


## Producing BayesTraits shell script
Invoke-Expression "$RFOLDER .\2_producing_BayesTraits_phyreg_shell.R"


## Running BayesTraits program with the automatically produced shell script
cd $BAYESTRFOLDER

Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope LocalMachine
Unblock-File -Path .\run_BayesTraits_phyreg_shell_windows.ps1
.\run_BayesTraits_phyreg_shell_windows.ps1


## Moving the resulting outputs into the main project folder
Move-Item -Path BayesTraits_phyreg_outputs_log_stone -Destination $PROJECTFOLDER
Move-Item -Path BayesTraits_phyreg_outputs_sch_trees -Destination $PROJECTFOLDER


## Processing BayesTraits output
cd $PROJECTFOLDER
Invoke-Expression "$RFOLDER .\3_processing_BayesTraits_phyreg_outputs.R"


## Opening automatically the resulting .pdf and .docx files
$CURRENTDATE=Get-Date -UFormat "%d%b%Y"
cd results_BayesTraits_output/"$CURRENTDATE"

Invoke-Expression "$ACROBATFOLDER Bayes_Factor_phyreg_heatmap_cluster.pdf"
Invoke-Expression "$ACROBATFOLDER Bayes_Factor_phyreg_heatmap.pdf"
Invoke-Expression "$WORDFOLDER BayesTraits_phyreg_output_table.docx"

