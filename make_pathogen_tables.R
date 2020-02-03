#!/usr/bin/env Rscript
library(optparse) # Read in package to read in command-line options.



# DESECRIPTION: Script to generate pathogens species table, diversity estimation tables, and pathogen detection / prevalence tables

# Pathogen analysis script
# Author: Olabiyi Aderemi Obayomi
# E-mail: obadbotanist@yahoo.com
# created: February 2018



version <- "1.0"

option_list <- list(
  
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help=paste0("Directory path containing all pathogens files generated from running the analyis on UNIX (required)",
                          "for example results_together/final_*"),
              metavar="filepath"),
  
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Tab-delimited file containing the metadata for the samples(required).",
              metavar="filepath"),
  
  make_option(c("-c", "--categories"), type="character", default=NULL,
              help=paste0("Tab-delimited file with the first column",
                          "corresponding to the categories in your the metadata(required)" ),
              metavar="filepath"),
  make_option(c("-f", "--funcFile"), type="character", default=NULL,
              help="path to the Pathogen_analysis_functions.R file provided (required)",
              metavar="filepath"),
  
  make_option(c("-o", "--outdir"), type="character", default= '.',
              help="Working directory to save output files and .RData file(required). default = current working directory" ,
              metavar="filepath"),
  
  make_option(c("--version"), action = "store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
  
)



opt_parser <- OptionParser(
  option_list=option_list,
  usage = "%prog [options] -i log1.txt,log2.txt",
  description = paste(
    " script to parfor pathogen analysis on the output files produced from running pathogena anlysis on UNIX",
    "\nit outputs a pathogen species table, pathogen detection tables  and the corresponding .RData",
    "USAGE:\n Rscript make_pathogen_tables.R  --indir results_together/ --metadata mapping.txt  
    --categories categories.tsv --funcFile Pathogen_analysis_functions.R --outdir outputdir/",
    sep="")
)

opt <- parse_args(opt_parser)

# Print out version if --version flag set.
if (opt$version) {
  cat("Wrapper version:", version, "\n")
  options_tmp <- options(show.error.messages=FALSE)
  on.exit(options(options_tmp))
  stop()
}

if(is.null(opt$indir)) {
  stop("paths to input files need to be set.")
}

if(is.null(opt$metadata)) {
  stop("paths to metadata table need to be set.")
}

if(is.null(opt$categories)) {
  stop("paths to categories table need to be set.")
}

if(is.null(opt$funcFile)) {
  stop("paths to the functions file must be set.")
}


# Set all necessary variables
working_directory <- opt$outdir

independent_Variables <- opt$categories # 'C:/Users/obayomi/Documents/16s_pathogen_analysis/test/categories.tsv'
independent_Variables <- read.table(file= independent_Variables, header= TRUE, sep="\t", stringsAsFactors = FALSE)
independent_Variables <- independent_Variables[,1] 

mapping_file <- opt$metadata
files_dir <- opt$indir 
functions_file <-  opt$funcFile

setwd(working_directory)

# Source the libraries and functions needed for the analysis
func_env <- new.env()
source(functions_file, local = func_env)
attach(func_env)




# Get the pathogen file names arranged from the lowest to the highest numerically
pathogen_files <- list.files(files_dir) #get all the files in the directory with all the result files
pathogen_files <- gsub(pattern = "__", 
                       replacement = ".",
                       x = mixedsort(gsub(pattern = "\\.", replacement = "__", x = pathogen_files)))

# Get the Treatment or sample names

# Get a list of split file names by substituting "final_" with "" and splitting the file name by "_" 
split_file_names <- strsplit(x = gsub(pattern = "final_", replacement = "", x = pathogen_files),
                             split = "_" )

# A character vector of sample or treatment names to be analysed as they appear in your mapping file

# Loop over the split_file_name list of lists getting the first element of each list
#and setting the respective treatment value to the right index in the vector of
#samples
samples <- map_chr(.x = split_file_names, .f = function(split_file_name) {split_file_name[1] } )
names(samples) <- samples 


#####################################################################################
### Perform pathogen analysis for each treatment or sample
# Turn off annoying warnings that could be ignored
options(warn = -1)
# Perform the analysis for the number of samples
Result_together <- map2(.x = samples, .y = str_c(files_dir, pathogen_files), 
                        .f = pathogen_analysis)


# Change to the parent directory
#setwd('../')
# Turn warnings back on
options(warn = 0)

########## Get the frequency tables for all the taxon levels together 
########## in a list called all_frequency_tables 

# Taxon levels vector
taxons <- c('genera','species', 'strains')
names(taxons) <- taxons

all_frequency_tables <- map(.x = taxons, .f = combine_freq_tables,
                            Treatments = samples,  Result_together = Result_together)



################################## Make taxon tables #####################
# Make raw and normalized taxon/otu tables for all the taxon levels and combine them
combined_taxon_tables <- make_taxon_tables(samples = samples, 
                                           taxons = taxons, 
                                           all_frequency_tables = all_frequency_tables)


# Get the independent variables dataframe
independent_variables.df <- get_independent_variables_df(mapping_file = mapping_file, 
                                                         independent_variables = independent_Variables)


# Get the species taxon table
taxon_table <- combined_taxon_tables$raw_tables$species
names(independent_Variables) <- independent_Variables
# Get the frequency of detection for each pathogen per treatment group
detect <- map(.x = independent_Variables, .f = presence_absence,
              taxon_table = taxon_table, independent_variables.df = independent_variables.df)



#####  Diversisty analysis
# Drop samples with no detects
species_table <-  taxon_table[,-c(which(colSums(taxon_table) == 0))]

# Transpose the species table
species_table_t <- t(species_table)

# Get the common ids between the species table and 
common.ids <- intersect(rownames(independent_variables.df),rownames(species_table_t))

# Get the common ids for the independent variables
independent_variables.df <- independent_variables.df[common.ids, , drop=FALSE]
species_table_t <- species_table_t[common.ids, , drop=FALSE]

# Estimate abundance based richness - observed species, chao1 and ACE
richness <- estimateR(species_table_t)
# Remove the standard error columns
richness <- richness[-c(3,5),]

# Estimate diversity
shanon_index <- diversity(species_table_t) #shannon
simpson_index <- diversity(species_table_t, index = 'simpson') #simpson

# Create a dataframe with the indentification for each sample attached to the species table
raw_table <- cbind(independent_variables.df,species_table_t)

# Create a diversity dataframe
diversity_table <- cbind(independent_variables.df,t(richness),
                         shannon = shanon_index, simpson = simpson_index)


# Get the means and SE of the diversity matrices
diversity_table_metrics <- c("S.obs",  "S.chao1", "S.ACE", "shannon", "simpson" )
names(diversity_table_metrics) <- diversity_table_metrics

diversity_result <- map(.x = independent_Variables, .f = create_mean_table, data= diversity_table )

# write the pathogen detection tables to an excel file for every category
walk2(.x = detect,.y = independent_Variables,
      .f = function(.x,.y){
        write.xlsx(x = .x, 
                   file = 'pathogen_detection.xlsx', 
                   sheetName = .y, 
                   col.names = T, 
                   row.names = F,
                   append = TRUE)
        
      }
)

# Write diversity tables for every category

walk2(.x = diversity_result ,.y = independent_Variables,
      .f = function(.x,.y){
        write.xlsx(x = as.data.frame(.x), 
                   file = 'diversity_tables.xlsx', 
                   sheetName = .y, 
                   col.names = T, 
                   row.names = F,
                   append = TRUE) 
      }
)

# write out the pathogen species table for future analysis either using my workflow or other sofware
write.table(x = taxon_table, file = "pathogen_species_table.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)


save.image("tables.RData")

