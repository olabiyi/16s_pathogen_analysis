#!/usr/bin/env bash

#Author: Olabiyi Aderemi Obayomi
#email: obadbotanist@yahoo.com
#created: February 2018
set -e
# This script will get the sample names from a QIIME mapping file
#and store them in a varriable called samples

echo "Provide the full path to your QIIME mapping file"
read map

echo "which column number contains your samples"
read col

#get the list of samples in a QIIME mapping file
#by getting the required field, removing the first line then
#sorting it, getting the unique sample names and saving them in a file 
#called sample ids
cut -f$col $map |sed '1d'| sort -u  > sample_ids.txt 

#assign the samples to a variable
samples=$(cat sample_ids.txt)

#convert the samples variable to an array
samples=($samples)
