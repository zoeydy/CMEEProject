#!/bin/bash

# Author:  Sam Turner sat19@ic.ac.uk
# Script:  run_MiniProject.sh
# Desc:    ties together other scripts into MiniProject workflow, resulting in compiled write up being written to ../report/
# Arguments: none
# Date: March 2020


printf "Install R packages if required.\n"
Rscript "install_packages.R"

printf "Preparing data.\n"
Rscript "data_preparation.R"

printf "Calculating initial values.\n"
python3 "initial_vals.py"

printf "Fitting models: may take up to 10 minites \n"
Rscript "model_fitting.R"

printf "Making demonstration plots.\n"
Rscript "demo_plots.R"


printf "Performing analysis.\n"
Rscript "analysis.R"

printf "Compiling report.\n"
bash "CompileLaTeX.sh" > junk.txt

mv write_up.pdf ../report/write_up.pdf

rm -f Rplots.pdf
rm -f *.sum 
rm -f .write_up.pdf 
rm -f .Rhistory
rm -f junk.txt
