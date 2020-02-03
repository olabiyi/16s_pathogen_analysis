#!/usr/bin/env Rscript

#if the package isn't already install it then load it
if (!require(vegan)){ install.packages('vegan'); library(vegan)}

if (!require(tidyverse)){ install.packages("tidyverse"); library(tidyverse)}

if (!require(xlsx)){ install.packages('xlsx'); library(xlsx)}



# DESECRIPTION: Script to generate pathogens species table, diversity estimation tables, and pathogen detection / prevalence tables

# Pathogen analysis script
# Author: Olabiyi Aderemi Obayomi
# E-mail: obadbotanist@yahoo.com
# created: February 2018



version <- "1.0"

option_list <- list(
  
  make_option(c("-i", "--input"), type="character", default=NULL,
              help=paste0("Tab-delimited file containing pathe pathogen counts per table (required)",
                          "for example results_together/final_*"),
              metavar="filepath"),
  
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Tab-delimited file containing the metadata for the samples(required).",
              metavar="filepath"),
  
  make_option(c("-c", "--categories"), type="character", default=NULL,
              help=paste0("Tab-delimited file with the first column",
                          "corresponding to the categories in your the metadata (required)" ),
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
    " script to perform diversity estimation on a count table",
    "\nit outputs diversity eastimate for each category to an excel file  and the corresponding .RData",
    "USAGE:\n Rscript make_pathogen_tables.R  --input results_together/ --metadata mapping.txt  
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

if(is.null(opt$input)) {
  stop("paths to input files need to be set.")
}

if(is.null(opt$metadata)) {
  stop("paths to metadata table need to be set.")
}

if(is.null(opt$categories)) {
  stop("paths to categories table need to be set.")
}


SE <- function(measure){sd(measure,  na.rm=TRUE)/sqrt(n())}

create_mean_table <- function(category, data){ 
  # Mean and Standard error
  mean_table <- data %>%  
    group_by(!!sym(category)) %>% 
    summarise_if(.predicate = is.numeric, .funs = list(mean=mean), na.rm=TRUE)
  SE_table <- data %>%  
    group_by(!!sym(category)) %>% 
    summarise_if(.predicate = is.numeric, .funs = list(SE=SE))                             
  
  inner_join(x=mean_table,y=SE_table,by=category)
  
}


# Set all necessary variables
working_directory <- opt$outdir

setwd(working_directory)

species_table <- opt$input # 'pathogen_species_table.txt' 
mapping_file <- opt$metadata # './mapping.txt'
independent_Variables <- opt$categories # 'C:/Users/obayomi/Documents/16s_pathogen_analysis/test/categories.tsv'
independent_Variables <- read.table(file= independent_Variables, header= TRUE, sep="\t", stringsAsFactors = FALSE)
independent_Variables <- independent_Variables[,1] 

#get the pathogens species table
species_table <- combined_taxon_tables$raw_tables$species

#drop samples with no detects
species_table <-  species_table[,-c(which(colSums(species_table) == 0))]




#transpose the species table
species_table_t <- t(species_table)


# Read in the mapping file
map <- read.table( file = mapping_file,
                   sep='\t', head=T, row.names=1, comment = '')


#get the independent variables to be tested
ind_var <- subset(x = map, select = independent_Variables )

#get the common ids between the species table and 
common.ids <- intersect(rownames(ind_var),rownames(species_table_t))

# Get the common ids for the independent variables
ind_var <- ind_var[common.ids,]
species_table_t <- species_table_t[common.ids,]
abund_species_table_t <- abund_species_table_t[common.ids,]
# Make rarefaction curves
rarecurve(species_table_t)

# Estimate abundance based richness - observed species, chao1 and ACE
richness <- estimateR(species_table_t)
# Remove the standard error columns
richness <- richness[-c(3,5),]


# Estimate diversity
shanon_index <- diversity(species_table_t) #shannon
simpson_index <- diversity(species_table_t, index = 'simpson') #simpson




# Create a diversity dataframe
diversity_table <- cbind(independent_variables.df,t(richness),
                         shannon = shanon_index, simpson = simpson_index)


# Get the means and SE of the diversity matrices
diversity_table_metrics <- c("S.obs",  "S.chao1", "S.ACE", "shannon", "simpson" )
names(diversity_table_metrics) <- diversity_table_metrics

diversity_result <- map(.x = independent_Variables, .f = create_mean_table, data=diversity_table)

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

