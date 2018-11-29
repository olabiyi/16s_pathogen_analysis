#!/usr/bin/env bash

#Author: Olabiyi Aderemi Obayomi
#Created: February 2018

#this script generates a list of whatever you want at anything between 90 - 100% identity to sequences in your given database
#make sure you pass the treatments to be analysed as arguements to the script seperated by spaces
#make sure you add QIIME to your $PATH variable
#the version of QIIME used for the analysis was QIIME 1.9.1
#the version of BLAST used for the analysis was ncbi-blast-2.7.1

#enter the necessary files and save them in variables
echo "what are you searching for: "
read search
echo "Enter the full path to the list of what to search for: "
read pathogens
echo "Enter the full path to the mapping file: "
read mapping_file
echo "Enter the full path to the blast database: "
read database
echo "Enter the full path to the fasta file of what you want to search for: "
read fasta
echo "Enter the column name in your mapping file to be analysed: "
read colName
echo "What will you like to call the file that will contain your combined final results? "
read results_together
echo "Which identity threshold do you want? "
echo "options are 0 to 10 for 90 to 100% identity, respectively"
read ID


#setting the variables
seqs=_$search.fasta
blast_out=_organisms.blast
list=_$search.list
organism_list=_organisms.list
unique=_unique
identity=_identity
final=final_
mappping_file=$(cat $mapping_file)
analyses=_Analysis

#run a blast search on every treatment passed as arguements to the script 
#e.g bash blast_sequences intact bussel potable water will analyse intact bussel potable and water
for treatment; do
   #remove the directory if it already exists and make a new directory for the treatment and then change to the directory
   rm -rf $treatment$analyses
   mkdir $treatment$analyses
   cd $treatment$analyses
    
    #Extract  sequences that belong to the treatment
   extract_seqs_by_sample_id.py -i $fasta -o $treatment$seqs \
   -m $mapping_file -s "$colName:$treatment"

   echo "blasting... $treatment"
   #blast the treatments's sequences against any blast formatted database, retrieving only the best hit
   blastn -db $database -query  $treatment$seqs -outfmt 0 \
   -num_descriptions 1 -num_alignments 1 -out $treatment$blast_out

   #get the list of organisms from the blast output with their identities
    grep ">"  $treatment$blast_out > $treatment$organism_list

   #get the list of the organisms' identity scores
   grep -E "Identities" $treatment$blast_out > $treatment$identity
   #paste the organism list and Identities together
   paste -d" " $treatment$organism_list $treatment$identity > $treatment$organism_list$identity 
   # get the list of organisms with identity scores greater than ID or 100%
        if [ $ID -eq 10 ]; then
             grep -E "(100%)" $treatment$organism_list$identity > $treatment$identity$organism_list
        else 
             grep -E "(9[$ID-9]%)|(100%)" $treatment$organism_list$identity > $treatment$identity$organism_list
        fi
   #create a file that will contain the list of what you are looking for
   touch $treatment$identity$list
   #search for what you are looking for within the organisms list
   grep -wFf $pathogens  $treatment$identity$organism_list > $treatment$identity$list
   ##get the pathogen list and parse it for readability i.e. without the sequence assertion number and removing anything 16S or 18S
   cut -d" " -f2-6  $treatment$identity$list |sed -e 's/1[68]S//g'|sort > $final$treatment$list

   #get the unique list
   sort -u $final$treatment$list > $treatment$unique$list
   #return to the parent directory
   cd ..
done

#####   organize the output files ##################
#remove the folders containing the combined together if they already exist
rm -rf $results_together
rm -rf unique

#create a directory that will contain the final results 
mkdir $results_together
mkdir unique
#put all the final and unique results files in one folder 
find . -type f -name "final_*" -exec mv -vuni '{}' "$results_together" ";" 2> error.log 
find . -type f -name "*unique*" -exec mv -vuni '{}' "unique/" ";" 2>> error.log

#change the extensions of the final and unique lists from .list to .txt so that they may be viewed 
#with a text editor outside unix such as notepad
rename .list .txt $results_together/*
rename .list .txt unique/*

