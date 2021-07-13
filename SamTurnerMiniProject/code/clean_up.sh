#!/bin/bash

# Author:  Sam Turner sat19@ic.ac.uk
# Script:  clean_up.sh
# Desc:    removes files created by running run_MiniProject, including write up - restoring MiniProject directory to original state.
# Arguments: none
# Date: March 2020


rm -f ../results/*.pdf
rm -f ../report/write_up.pdf
rm -f ../data/LogisticGrowthDataLogClean.csv
rm -f ../data/inits.csv
rm -f ../data/*.rds