#!/bin/bash

mkdir "../artifacts/pwy_pdfs/$1"
for file in $(ls "../artifacts/pwy_svgs/$1")
do
    pref=${file%.*}
    inkscape --export-pdf="../artifacts/pwy_pdfs/$1/$pref.pdf" "../artifacts/pwy_svgs/$1/$file"
done