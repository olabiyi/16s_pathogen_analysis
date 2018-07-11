#! /bin/bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N soil_pathogens_search
#$ -pe shared 24

# modify the environment for Qiime
#export PATH="/storage16/app/bioinfo/blast-2.2.22/bin:/fastspace/bioinfo_apps/cdbfasta:/fastspace/bioinfo_apps/microbiomeutil-r20110519/ChimeraSlayer:/fastspace/bioinfo_apps/qiime/usr/local/bin:/fastspace/bioinfo_apps/qiime/local/bin:/bin:/usr/bin" PYTHONPATH="/fastspace/bioinfo_apps/qiime/usr/local" HDF5_DIR="/fastspace/bioinfo_apps/qiime/local/hdf5" RDP_JAR_PATH="/fastspace/bioinfo_apps/qiime/local/rdp_classifier_2.2/rdp_classifier-2.2.jar"QIIME_CONFIG_FP="/fastspace/bioinfo_apps/qiime/qiime_config_biores"


# source my .bashrc file
source /storage/users/gilloro/.bashrc

#get the list of samples in a QIIME mapping file
#by getting the required field, removing the first line then
#sorting it, getting the unique sample sampl names and saving them in a file called sample ids
cut -f1 /gpfs0/biores/users/gilloro/Biyi/community_analysis/soil/00.mapping_file_validation/barrier_exp/protozoa/corrected/nonirr_barrier_protozoa_mapping_corrected.txt |sed '1d'| sort -u  > sample_ids.txt 

#put the samples in a variable
samples=$(cat sample_ids.txt)

#convert to an array
samples=($samples)

#this script generates a list of whatever you want at anything between 90 - 100% identity to sequences in your given database
#make sure you pass the treatments you want to be analysed as arguements to the script seperated by spaces
#make sure to add QIIME to your $PATH variable


#Enter the necessary  variables
search=parasites
pathogens=/gpfs0/biores/users/gilloro/Biyi/pathogen_analysis/human_parasite_list.txt
#/gpfs0/biores/users/gilloro/Biyi/community_analysis/soil/00.mapping_file_validation/soilType_exp/bacteria/corrected/nonirr_soilType_bacteria_mapping_corrected.txt
database=/gpfs0/biores/users/gilloro/Biyi/protozoa_DB/18S_all
fasta=/gpfs0/biores/users/gilloro/Biyi/community_analysis/soil/pathogens99/protozoa/parasites.fasta
colName=sample
results_together=results_together
ID=8
seqs=_$search.fasta
blast_out=_organisms.blast
list=_$search.list
organism_list=_organisms.list
unique=_unique
identity=_identity
final=final_
#mappping_file=$(cat /gpfs0/biores/users/gilloro/Biyi/community_analysis/soil/00.mapping_file_validation/soilType_exp/protozoa/corrected/nonirr_soilType_protozoa_mapping_corrected.txt)
analyses=_Analysis

#loop for as many treatments passed as arguements to the script e.g bash blast_sequences intact bussel potable water
#will analyse for intact bussel potable and water
for treatment in ${samples[*]}; do
   #remove the directory if it already exists and make a new directory for the treatment and then change to the directory
   rm -rf $treatment$analyses
   mkdir $treatment$analyses
   cd $treatment$analyses
    #########################################################################################################
    #Extract  sequences that belong to the treatment
   extract_seqs_by_sample_id.py -i $fasta -o $treatment$seqs \
   -m /gpfs0/biores/users/gilloro/Biyi/community_analysis/soil/00.mapping_file_validation/barrier_exp/protozoa/corrected/nonirr_barrier_protozoa_mapping_corrected.txt -s "$colName:$treatment"

   echo "blasting... $treatment"
   #blast the treatments's sequences against any blast formatted database retrieving on the best hit
   blastn -db $database -query  $treatment$seqs -outfmt 0 \
   -num_descriptions 1 -num_alignments 1 -out $treatment$blast_out -num_threads 24

   #get the list of organisms from the blast output along with their identities
    grep ">"  $treatment$blast_out > $treatment$organism_list

   #get the list of the organisms identity scores
   grep -E "Identities" $treatment$blast_out > $treatment$identity
   #paste the organism list and its Identities together
   paste -d" " $treatment$organism_list $treatment$identity > $treatment$organism_list$identity 
   # get the list of organisms with identity scores greater than ID or 100%
        if [ $ID -eq 10 ]; then
             grep -E "(100%)" $treatment$organism_list$identity > $treatment$identity$organism_list
        else 
             grep -E "(9[$ID-9]%)|(100%)" $treatment$organism_list$identity > $treatment$identity$organism_list
        fi
   #open a file that will contain the list of what you looking for
   touch $treatment$identity$list
   #search for what you looking for within the organisms list
   grep -wFf $pathogens  $treatment$identity$organism_list > $treatment$identity$list
   ##get the pathogen list and sort them in a more readable manner without the sequence assertion number and remove anything 16S or 18S
   cut -d" " -f2-6  $treatment$identity$list |sed -e 's/1[6,8]S//g'|sort > $final$treatment$list

   #get the unique list
   sort -u $final$treatment$list > $treatment$unique$list
   #change back to the previous/parent directory
   cd ..
done

#####   trying to organize the output files ##################
#remove the folders containing the results together if they already exist
rm -rf $results_together
rm -rf unique

#make a directory that will contain all the final results tohether
mkdir $results_together
mkdir unique
#put all the final and unique results files in one folder for each
find . -type f -name "final_*" -exec mv -vuni '{}' "$results_together" ";" 2> error.log 
find . -type f -name "*unique*" -exec mv -vuni '{}' "unique/" ";" 2>> error.log

#change the extensions of  the final and unique lists from .list to .txt that they may be viewd with a text editor out unix like notepad
rename .list .txt $results_together/*
rename .list .txt unique/*
