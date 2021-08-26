#!/bin/bash


pdflatex main.tex
pdflatex main.tex
name=$(echo "ZongyiHu_Disertation.tex" | cut -f 1 -d '.')
bibtex $name
bibtex $name
echo "$name"
pdflatex main.tex
pdflatex main.tex
evince $name.pdf 

##Cleanup
rm *~
rm *.aux
rm *.dvi
rm *.log
rm *.nav
rm *.out
rm *.snm
rm *.toc
rm *.bbl
rm *.blg
rm *.fdb_latexmk
rm *.fls