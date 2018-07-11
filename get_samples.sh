#!/bin/bash

#Author: Olabiyi Aderemi Obayomi

# This script will get the sample names from a QIIME mapping file
#and store them a varriable called samples

echo "Provide the full path to your QIIME mapping file"
read map

echo "which column number contains your samples"
read col

#get the list of samples in a QIIME mapping file
#by getting the required field, removing the first line then
#sorting it, getting the unique sample sampl names and saving them in a file called sample ids
cut -f$col $map |sed '1d'| sort -u  > sample_ids.txt 

#put the samples in a variable
samples=$(cat sample_ids.txt)

#convert to an array
samples=($samples)
