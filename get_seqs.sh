#!/bin/bash

#Author: Olabiyi Aderemi Obayomi
#year:2018

#This scripts generates a fasta file of whatever you are looking for from a QIIME formatted OTU table

#I suggest that you open a directory for the analysis
#if you will be analysing for pathogens for exemple make sure you have the following files

#1 human_pathogens.txt - a text file containing the list of human pathogen genera to be found
#2 Or  plant_pathogens.txt - a text file containing the list of plant pathogen genera to be found
#3 otu_table.biom - a biom table containg the list of all the OTUs and there  abundance values
#4 seqs_otus.txt - a tab delimited text file containing the map of all the OTUs to their related sequences
#5  seqs.fna - a fasta file generated after quality filtering and demultiplexing which contains the list of all sequences

#Get the user's inputs and save them in their corresponding variables 
echo "What are you looking for? :"
read search 
echo "Enter the full path to the file containing the comma separated list of what you looking for: "
read list
echo "Enter the path to the OTU.biom file: "
read  otu_biom
echo "Enter the path to the seqs_otus.txt file: "
read  seqs_otus
echo "Enter the path to the seqs.fna file: "
read seqs_fna
#pathogen_seqs=$(cat $treatment)
list=$(cat $list)
#otu_biom=$(<otu_biom)
#seqs_otus=$(<seqs_otus)
#seqs_fna=$(<seqs_fna)
#filtering out the biom tables for what you are looking for e.g. human pathogens
filter_taxa_from_otu_table.py -i $otu_biom -o $search.biom -p $list

#convert otu biom table to tsv
biom convert     -i  $search.biom   -o  $search.biom.tsv     --to-tsv     --header-key taxonomy     --output-metadata-id "Consensus Lineage"

#get the  otu names
#cut  the first field, remove the first 2 lines, substitute newline/end of line with space and then substitute space with |
#then remove the last character '|' from the last line,replace newline with space such that all the  otu names are line-up on oneline
#then replace all spaces with an empty string so the search will be compact
#otus=$(cut -f1 $search.biom.tsv|sed -e 1,2d|sed 's/$/ /'|sed 's/ /|/' | sed '$ s/.$//'|tr '\n' ' '|sed 's/ //g')

#get the  otu names
#cut  the first field, remove the first 2 lines
otus=$(cut -f1 $search.biom.tsv|sed -e 1,2d)
#convert the otus to an array
otus=($otus)

#loop over every element in the otus array
for otu in ${otus[*]}; do 
    #creating the map by greping the exact word match
    grep -wE "$otu" $seqs_otus  >> $search_map.txt 
done

#get the fasta sequences for what you are looking for e.g. human pathogens
filter_fasta.py -f $seqs_fna -o $search.fasta -m  $search_map.txt

