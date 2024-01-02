#!/bin/bash

cd "../artifacts/pwy_pdfs/$1"
pwd
input_pdfs=$(ls)
pdfunite $input_pdfs $1.pdf