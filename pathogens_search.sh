#! /bin/bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N soil_pathogens_search
#$ -pe shared 24



###### This script runs the first 3 steps in the 'how_to_run_pathogenAnalysis.txt' file in one go on qsub ######


#source your .bashrc file - it should look like the lines below
#this is just to make available the necessary scripts for the analysis
#pathogen analysis scripts
#export PATH=/path/to/the/folder/containing/the/scripts/for/pathogen_analysis:$PATH
#qiime
#export PATH=/path/to/your/qiime/usr/local/bin/:$PATH
#NCBI blast
#export PATH=/path/to/your/ncbi-blast-2.7.1+/bin/:$PATH

source /path/to/your/.bashrc

#get the list of samples in a QIIME mapping file
#by getting the required field, removing the first line then
#sorting it, getting the unique sample names and saving them in a file called sample_ids.txt
cut -f1 /path/to/your/mapping_file.txt |sed '1d'| sort -u  > sample_ids.txt 

#assign the samples to a variable
samples=$(cat sample_ids.txt)

#convert to samples variable to an array
samples=($samples)



#Enter the following necessary  variables MANUALLY
search=pathogens
pathogens=/path/to/the/human_pathogens_list.txt
database=/path/to/your/database/refseq_16s
fasta=/path/to/the/pathogens.fasta
colName=sample
results_together=results_together
ID=9

#leave these set variables as they are
seqs=_$search.fasta
blast_out=_organisms.blast
list=_$search.list
organism_list=_organisms.list
unique=_unique
identity=_identity
final=final_
analyses=_Analysis

#loop for as many treatments passed as arguements to the script e.g bash blast_sequences intact bussel potable water
#will analyse for intact bussel potable and water
#run a blast search on every treatment passed as arguements to the script 
#e.g bash blast_sequences intact bussel potable water will analyse intact bussel potable and water
for treatment; do
   #remove the directory if it already exists and make a new directory for the treatment and then change to the directory
   rm -rf $treatment$analyses
   mkdir $treatment$analyses
   cd $treatment$analyses
    
    #Extract  sequences that belong to the treatment
   extract_seqs_by_sample_id.py -i $fasta -o $treatment$seqs \
   -m /path/to/your/mapping_file.txt -s "$colName:$treatment"

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