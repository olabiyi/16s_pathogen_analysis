#!/usr/bin/env bash

#Author: Olabiyi Aderemi Obayomi
#created: February 2018
#email: obadbotanist@yahoo.com

#This scripts generates a fasta file of whatever you are looking for, from a QIIME formatted OTU table

#I suggest that you open a directory for the analysis
#if you will be analysing for pathogens for example make sure you have the following files

#1a human_pathogens.txt - a text file containing the list of human pathogen genera to be found
# Or  1b. plant_pathogens.txt - a text file containing the list of plant pathogen genera to be found
#2 otu_table.biom - a biom table containg the list of all the OTUs and there  abundance values
#3 seqs_otus.txt - a tab delimited text file containing the map of all the OTUs to their related sequences
#4  seqs.fna - a fasta file generated after quality filtering and demultiplexing which contains the list of all sequences

#Get the user's inputs and save them in their respective variables 
echo "What are you looking for? :"
read search 
echo "Enter the full path to the file containing the comma separated list of what you looking for: "
read list
echo "Enter the full path to the OTU.biom file: "
read  otu_biom
echo "Enter the full path to the seqs_otus.txt file: "
read  seqs_otus
echo "Enter the full path to the seqs.fna file: "
read seqs_fna

list=$(cat $list)

#filter out the OTUs of what you are looking for from the given OTU table e.g. human pathogens
filter_taxa_from_otu_table.py -i $otu_biom -o $search.biom -p $list

#convert biom formattted OTU table to tsv
biom convert     -i  $search.biom   -o  $search.biom.tsv     --to-tsv     --header-key taxonomy     --output-metadata-id "Consensus Lineage"



#get the OTU names
#cut the first field, remove the first 2 lines
otus=$(cut -f1 $search.biom.tsv|sed -e 1,2d)
#convert the otus to an array
otus=($otus)

#find each OTU in the OTU map
for otu in ${otus[*]}; do 
    #create the new OTU map by searching for the exact word match
    grep -wE "$otu" $seqs_otus  >> $search_map.txt 
done

#get the fasta sequences for what you are looking for e.g. human pathogens fasta sequences
filter_fasta.py -f $seqs_fna -o $search.fasta -m  $search_map.txt

