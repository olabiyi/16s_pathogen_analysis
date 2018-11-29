#pathogen analysis functions and library script
#Author: Olabiyi Aderemi Obayomi
#E-mail: obadbotanist@yahoo.com
#created: February 2018

#How to install bioconductor packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("packageName")


#loading the required packages
#load ggpubr for plotting and statistics
library('ggpubr')
#package for statistical decriptions
library('DescTools')
#packages for microbiome analysis
library('vegan')
library('phyloseq')
library('DESeq')
library('edgeR')
library('metagenomeSeq')
library('gtools')
library('ComplexHeatmap') # mapping heat maps
library('xlsx') #for writing to excel files
library('data.table')
# make sure that you are using the correct architecture of java
#i.e 64 bit java with 64 bit R version or 32 bit java with 32 bit R version
library("rJava", lib.loc="~/R/win-library/3.4") #load rJava

##############################################################################################
#load in the functions
#main pathogen analysis function
pathogen_analysis <- function(treatment, pathogen_file){
	library('DescTools')
  #Get the frequency of the pathogens detected in a treatment / sample at the genus, specie and starin level. 
  
  # ~ treatment - is a string specifying the name of the treatment / sample being analysed
  # ~ pathogen_file - This is the output file generated from the pathogen search in unix using qiime and bash scripts
  
  #check if the file is empty
  if (file.size(pathogen_file) > 0){
    #get the maximum number of fields in the file
    Ncol <- max(count.fields(pathogen_file, sep = " "))
    #read in the file while removing the excess number of columns using the colClasses arguement
    #row.name = NULL forces row numbering
    pathogens.df <- read.table(pathogen_file, sep = " ", 
                                colClasses = c(rep("character", 5), rep("NULL", (Ncol -5))), 
                                fill = T, row.names = NULL, header = F)

    #convert the columns to character in preparation for pasting
    pathogens.df$V1 <- as.character(pathogens.df$V1)
    pathogens.df$V2 <- as.character(pathogens.df$V2)
    pathogens.df$V3 <- as.character(pathogens.df$V3)
    pathogens.df$V4 <- as.character(pathogens.df$V4)
    pathogens.df$V5 <- as.character(pathogens.df$V5)
    
    #create dataframes for each taxon level by pasting the appropriate columns together
    genera.df <- pathogens.df$V1
    species.df <- paste(pathogens.df$V1,pathogens.df$V2, sep = " ")
    strains.df <- paste(pathogens.df$V1,pathogens.df$V2,pathogens.df$V3,pathogens.df$V4,pathogens.df$V5, sep = " ")
    
    #combine the vectors for each taxon level in one dataframe
    pathogens.df <- data.frame(genera = genera.df, species = species.df, strains = strains.df)
    
    #change the vector columns to factors
    pathogens.df$genera <- as.factor(pathogens.df$genera)
    pathogens.df$species <- as.factor(pathogens.df$species)
    pathogens.df$strains <- as.factor(pathogens.df$strains)
    
    #get the frequency tables 
    genera.df <- Freq(pathogens.df$genera, ord = "desc")
    colnames(genera.df) <- c("genera","freq","perc","cumfreq","cumperc")
    species.df <- Freq(pathogens.df$species, ord = "desc")
    colnames(species.df) <- c("species","freq","perc","cumfreq","cumperc")
    strains.df <- Freq(pathogens.df$strains, ord = "desc")
    colnames(strains.df) <- c("strains","freq","perc","cumfreq","cumperc")
  }else{
    genera.df <- data.frame(genera = "not_detected", freq = 0, perc = 0, cumfreq = 0 , cumperc = 0)
    species.df <- data.frame(species = "not_detected", freq = 0, perc = 0, cumfreq = 0 , cumperc = 0)
    strains.df <- data.frame(strains = "not_detected", freq = 0, perc = 0, cumfreq = 0 , cumperc = 0)
  }
  #tag the the output files with the name of the treatment
  #genera
  genera.df <- data.frame(Treatment = rep(treatment, times = length(genera.df[,1])),genera.df)
  
  #species
  species.df <- data.frame(Treatment = rep(treatment, times = length(species.df[,1])),species.df)
  
  #strains
  strains.df <- data.frame(Treatment = rep(treatment, times = length(strains.df[,1])),strains.df)
  
  
  #combine them in one list to be returned by the function
  result <- list(genera.df, species.df, strains.df)
  
  names(result) <- c('genera','species','strains')
  return(result)
}
#function for combining frequency tables
combine_freq_tables <- function(Treatments,Result_together,frequency_table_name){
	#parameters
  # ~ Treatments -
  # ~ Result_together -
  # ~ frequency_table_name -
  
  #initialize the frequency dataframe
  frequency_table <- data.frame()
  #concatenate the frequency table for all the treatments for each taxonomic level
  for (treat in 1:length(Treatments)){
    treatment_frequency <- Result_together[[Treatments[treat]]][[frequency_table_name]]
    frequency_table <- rbind(frequency_table,treatment_frequency)
  }
  return(frequency_table)
}

#function for drawing pathogen plots
draw_pathogen_plots <- function(Data,main, x, y,Xlab, YLab, plotType = 'box', group) {
	library('ggpubr')
 #function to make pathogen plots
  # ~ Data - is a dataframe for plotting - note that the data for the box plots needs to be stacked
  # ~ main -  is a string specifying the title of the plot
  if (plotType == 'bar'){
    ##make relative abunbance bar charts
    plot.gg <- ggbarplot(data = Data, 
                         x = x, 
                         y = y, 
                         fill = group, 
                         xlab = Xlab,
                         ylab = YLab,
                         title = main)
    plot.gg <- plot.gg + font("xy.text", size = 14, face = "bold",color = "black")
    
    
  }else {
    #making relative abundance boxplots faceted by pathogens
    plot.gg <- ggboxplot(data = Data, 
                         x = x, 
                         y = y,
                         color = x, 
                         add = "mean_se", 
                         facet.by = group,
                         xlab = Xlab,
                         ylab = YLab,
                         title = main)
    #remove tick labels so the plots will look less congested
    plot.gg <- plot.gg + rremove("ticks") 
    plot.gg <- plot.gg + rremove("xy.text")
  }
  #remove the legend title
  plot.gg <- ggpar(plot.gg, legend.title = "")
  print(plot.gg <- plot.gg + 
          font("xlab", size = 16, face = "bold",color = "black")+
          font("ylab", size = 16, face = "bold",color = "black"))#
  
}

# function for making OTU tables
make_OTU_table <- function(taxon_level, samples,samples_column,
                           frequency_column, all_frequency_tables) {
#This funtion generates taxon tables from a set of frequency tables
	###### parameters ##########
  # ~ taxon_level - is a string specifying the taxon level for which the OTU table should
  #be generated e.g genera - only three possiblities c('genera', 'species', 'strains')
  # ~ samples -  is a character vector of the sample names i.e the samples you got their results from blast
  # ~ sample_columns - is a string that specifies the name of the sample column in all_frequency table
  #that contains the sample names usually will be 'Treatment'
  # ~ frequency_column - is a string that specifies the frequency column in all_frequency_tables
  #that contains the frequency for each pathogen usually 'freq'
  # ~ all_frequency_tables -  is a  frequency table for all the taxon levels together
  
  #remove the sample in which no pathogen was detected in, from the frequency table
  not_detected_index <- grep(pattern = "not_detected", x = all_frequency_tables[[taxon_level]][[2]])
  if (length(not_detected_index) > 0){
    all_frequency_tables[[taxon_level]]<- all_frequency_tables[[taxon_level]][-(not_detected_index),]
  }
  #get the unique taxons for all the samples
  unique_taxon <- unique(all_frequency_tables[[taxon_level]][taxon_level])
  #create the taxon matrix filled with zeros
  taxon_table <- matrix(data = 0,nrow = nrow(unique_taxon), ncol = length(samples))
  rownames(taxon_table) <- unique_taxon[,1] #name the rows
  colnames(taxon_table) <- samples #name the columns
  #loop over every sample and within each sample loop over every pathogen
  for (Sample in samples){
    for (pathogen in rownames(taxon_table)){ 
      #check if the present pathogen isn't in the present sample if yes, return the value zero
      check <- all_frequency_tables[[taxon_level]][all_frequency_tables[[taxon_level]][samples_column] == Sample & all_frequency_tables[[taxon_level]][taxon_level] == pathogen,frequency_column]
      if (length(check) == 0){
        taxon_table[pathogen,Sample]<- 0
        
      }else{
        taxon_table[pathogen,Sample]<- all_frequency_tables[[taxon_level]][all_frequency_tables[[taxon_level]][samples_column] == Sample & all_frequency_tables[[taxon_level]][taxon_level] == pathogen,frequency_column]
      }
    }
  }
  return(taxon_table)
}
#normalization functions
#===Normalization functions===#
#Reference: http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003531
#edgeR-TMM,UQ-logFC, "RLE"
edgeRnorm <- function(physeq, ...){
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){
    physeq <- t(physeq)
  }
  x = as(otu_table(physeq), "matrix")
  # See if adding a single observation, 1, 
  # everywhere (so not zeros) prevents errors
  # without needing to borrow and modify 
  # calcNormFactors (and its dependent functions)
  # It did. This fixed all problems. 
  # Can the 1 be reduced to something smaller and still work?
  x <- x + 1
  # Now turn into a DGEList
  y <- edgeR::DGEList(counts=x, remove.zeros=TRUE)
  # Perform edgeR-encoded normalization, using the specified method (...)
  z <- edgeR::calcNormFactors(y, ...)
  # A check that we didn't divide by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
  }
  return(z)
}

#DESeqVS
deseq_varstab <- function(physeq, sampleConditions=rep("A", nsamples(physeq)), ...){
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){
    physeq <- t(physeq)
  }
  x = as(otu_table(physeq), "matrix")
  # The same tweak as for edgeR to avoid NaN problems
  # that cause the workflow to stall/crash.
  x <- x + 1
  cds <- newCountDataSet(x, sampleConditions)
  # First estimate library size factors
  cds <- estimateSizeFactors(cds)
  # Variance estimation, passing along additional options
  cds <- estimateDispersions(cds, ...)
  # Determine which column(s) have the dispersion estimates
  dispcol <- grep("disp\\_", colnames(fData(cds)))
  # Enforce that there are no infinite values in the dispersion estimates
  if( any(!is.finite(fData(cds)[, dispcol])) ){
    fData(cds)[which(!is.finite(fData(cds)[, dispcol])), dispcol] <- 0.0
  }
  vsmat = exprs(varianceStabilizingTransformation(cds))
  otu_table(physeq) <- otu_table(vsmat, taxa_are_rows=TRUE)
  return(physeq)
}

#Proportion
proportion <- function(physeq){
  # Normalize total sequences represented
  normf = function(x, tot=max(sample_sums(physeq))){ tot*x/sum(x) }
  physeq = transform_sample_counts(physeq, normf)
  # Scale by dividing each variable by its standard deviation.
  #physeq = transform_sample_counts(physeq, function(x) x/sd(x))
  # Center by subtracting the median
  #physeq = transform_sample_counts(physeq, function(x) (x-median(x)))
  return(physeq)
}
#Rarefy
randomsubsample <- function(physeq, smalltrim=0.15, replace=TRUE){
  # Set the minimum value as the smallest library quantile, n`smalltrim` 
  samplemin = sort(sample_sums(physeq))[-(1:floor(smalltrim*nsamples(physeq)))][1]
  physeqr = rarefy_even_depth(physeq, samplemin, rngseed=TRUE,
                              replace=replace, trimOTUs=TRUE)
  return(physeqr)
}
#/===Normalization functions===#



normalize_otu_table <- function(OTU_table, method='rarefy') {
  #This function performs normalization on a given count otu table OTU_Table
  #Don't forget to load the libraries and the functions necessary 
  #a.k.a c('phyloseq','metagenomeSeq ','EdegeR', 'vegan')
  # ~ OTU_table - a counts table with samples as columns and rows as features
  # ~ method  - is a string specifying the normalization to be performed on the OTU_table
  #can either be c('proportion', 'RLE', 'logUQ','TMM', 'DESeqVS', 'CSS', 'rarefy')
  #it performs normalization by rarefying by default  
  
  
  #change OTU table into a phyloseq object or MRexperiment obj for metagenomeSeq 
  phylose_otu_table<- phyloseq(otu_table(OTU_table, taxa_are_rows = T))
  meseq_obj<- newMRexperiment(OTU_table + 0.01)
  
  #perform normalizations
  if (method == 'proportion'){
    #normalization by proportions  
    otu_table.norm <- as.data.frame(proportion(physeq = phylose_otu_table))
  }else if (method == 'RLE') {
    #normalization by relative log expression
    otu_table.norm <- otu_table.norm<- as.data.frame(edgeRnorm(physeq = phylose_otu_table, 
                                                               method= "RLE")$counts)
  }else if (method == 'logUQ') {
    #normalization by upper-qauntile log transformation
    otu_table.norm<- as.data.frame(edgeRnorm(physeq = phylose_otu_table, 
                                             method= "upperquartile")$counts)
  }else if (method == 'TMM') {
    #normalization by weighted trimmed mean of M-values
    otu_table.norm<- as.data.frame(edgeRnorm(physeq = phylose_otu_table, 
                                             method= "TMM")$counts)
  }else if (method == 'DESeqVS'){
    # normalization by DESeq variance stabilizing transformations
    otu_table.norm<- as.data.frame(deseq_varstab(phylose_otu_table, 
                                                 method="blind", 
                                                 sharingMode="maximum",
                                                 fitType="local"))
  }else if (method == 'CSS'){
    # normalization by Metagenomeseq cummulative sum 
    cumNorm(meseq_obj)
    otu_table.norm<- as.data.frame(MRcounts(meseq_obj, norm = T, log = T))
    
  }else if (method == 'rarefy'){
    #normalization by rarefying
    otu_table.norm <- as.data.frame(randomsubsample(phylose_otu_table, replace = F))
    
  }else{
    print("incorrect option choosen: options are 'proportion', 'RLE', 'logUQ','TMM', 'DESeqVS', 'CSS', 'rarefy' ")
    
  }
  return(otu_table.norm)
}

#function to remove the samples / treatments that pathogens were not detected in 
remove_zero_samples_from_taxon_table <- function(otu_table) {

# ~ otu_table - is a taxon table with samples as columns and pathogens as rows
  #intialize count variables
  count <- 0 #to count the number of iteration; 
  i <- 0 #to count the number of empty samples
  column_number <- c() #column number vector
  #find the empty samples
  for (sample in colnames(otu_table)){
    count <- count + 1
    sample_sum <- sum(otu_table[,sample])
    if (sample_sum == 0){
      i <- i + 1
      column_number[i] <- count
    }
  }
  #remove the empty samples
  otu_table <- otu_table[,-column_number]
  return(otu_table)
}


#make taxon/otu tables for all the taxon levels
make_taxon_tables <- function(samples,taxons,all_frequency_tables, method = 'logUQ') {
  #This function makes a taxon table from frequency tables supplied which is the output
  #of the combine frequency table function. It outputs raw taxon and normalized taxon tables
  # ~ samples - is a string vector that specifies the sample names in the frequency table
  # ~ taxons - is a string vector that specifies the taxon names which are basically the names
  #of the corresponding frequency table
  # ~ all_frequency_tables - is a list that contains the frequency tables of all the taxons together
  # ~ method - is a string that specifies the normalization method, can be any of 
  #c('proportion', 'RLE', 'logUQ','TMM', 'DESeqVS', 'CSS', 'rarefy')
  
  #intialize the list variables in which the results will be strored    
  combined_taxon_tables <- vector(mode = 'list', length = 2) #both raw and normalized combined
  taxon_tables <- vector(mode = 'list', length = length(taxons)) #raw table
  otu_table.norm <- vector(mode = 'list', length = length(taxons))#normalized table
  
  #loop to make raw and normalized taxon tables for each taxon
  for (i in 1:length(taxons)){
    #otu tables
    taxon_tables[[i]]<- make_OTU_table(taxon_level = taxons[i],
                                       samples = Treatments, 
                                       samples_column = 'Treatment', 
                                       frequency_column = 'freq', 
                                       all_frequency_tables = all_frequency_tables)
    #normalize the OTU table
    otu_table.norm[[i]] <- normalize_otu_table(OTU_table = taxon_tables[[i]], method = method)
  }
  #name the lists of tables
  names(taxon_tables) <- taxons 
  names(otu_table.norm) <- taxons
  combined_taxon_tables[[1]] <- taxon_tables; combined_taxon_tables[[2]] <- otu_table.norm
  names(combined_taxon_tables)<- c('raw_tables', 'normalized_tables' )
  return(combined_taxon_tables)
}

#A function to get the independent variables from a dataframe
get_independent_variables_df <- function(mapping_file, independent_variables){
  #This function extracts the independent variables from a dataframe
  
  # ~ mapping_file - a string speciying that path to your mapping file
  # ~ independent - is a string vector of independent or factor variavles as they appear in the mapping file
  
  #import your qiime formatted mapping file
  map <- read.table( file = mapping_file, sep='\t', head=T, row.names=1, comment = '') 
  #get the independent variables
  independent_variables.df <- subset(x = map, select = independent_variables)
  return(independent_variables.df)
}

#A function to generate table for statistics and plotting for a given independent variable
generate_stats_table <- function(taxon_table, mapping_file, independent_variables){
 #This function genrates a dataframe with the frequencies and samples in the taxon table
 #mapped to the appropriate metadata for a given independent / factor variable
 
 # ~ taxon table - is a matrix or dataframe of the count of taxons per sample
  # ~ mapping_file - is the metadata file containing each sample description
  # ~ independent_variables - is a string vector  specifying the 
  #the independent variables in the supplied mapping file
  
  #import your qiime formatted mapping file
  map <- read.table( file = mapping_file, sep='\t', head=T, row.names=1, comment = '') 
  taxon_table_t <- t(taxon_table) #transpose the OTU table
  #find the overlap between the sample names of the mapping file and the taxon table
  common.ids <- intersect(rownames(map), rownames(taxon_table_t ))
  # get just the overlapping samples
  taxon_table_t <- taxon_table_t[common.ids,]
  map <- map[common.ids,]
  #relative abundance normalization per sample
  taxon_norm <- sweep(taxon_table_t, 1, rowSums(taxon_table_t),'/')
  #get the independent variables
  independent_variables.df <- subset(x = map, select = independent_variables)
  #generate a statistics table for plotting and statistics
  stats_table <- cbind(independent_variables.df,taxon_norm)
  return(stats_table)
}

#function to generate a presence-absence dataframe for the
#pathogens detected per sample or treatment group in a taxon table
presence_absence <- function(taxon_table, independent_variables.df, variable){
# This function generates a presence-absence dataframe for the pathogens 
# detected per sample or treatment group in a taxon table
# ~ taxon_table - is a taxon table dataframe with samples as columns and pathogens as rows
# ~ independent_variables -  is an independent_variables dataframe generated using the 
#get_independent_variables_df function
# ~ variable - a string specifying the independent variable in the independent_variables.df 
#dataframe to be analysed
 
  taxon_table_t <- t(taxon_table) #transpose the OTU table
  #find the overlap between the samples in the independent variables dataframe and the taxon table
  common.ids <- intersect(rownames(independent_variables.df), rownames(taxon_table_t ))
  columns <- colnames(independent_variables.df)
  # get just the overlapping samples
  taxon_table_t <- taxon_table_t[common.ids,]
  independent_variables.df <- independent_variables.df[common.ids,]
  
  if (is.factor(independent_variables.df) == TRUE){
    independent_variables.df <- as.data.frame(independent_variables.df)
    rownames(independent_variables.df)<- common.ids
    colnames(independent_variables.df)<- columns
  }
  
  #find the pathogens detected per sample
  detect_table <- taxon_table_t > 0
  detect_table2 <- data.frame(independent_variables.df,detect_table)
  #get the total number of detects per treatment
  detect <- aggregate(x = detect_table, by = list(detect_table2[,variable]), FUN = sum)
  colNames <- colnames(detect)
  colNames[1] <- variable
  colnames(detect) <- colNames
  
  #get the number of samples per treatment group
  number_of_samples_per_treatment <- aggregate(x = detect_table,
                                               by = list(detect_table2[,variable]), 
                                               FUN = length)[,1:2]
											   
  colnames(number_of_samples_per_treatment) <- c(variable,'number_of_samples')
  result <- data.frame(number_of_samples_per_treatment, detect[,2:ncol(detect)])
  return(result)
}
#function to make PCOA, NMDS and box plots for the pathogens
ordination_plots <- function(taxon_table, mapping_file, independent_variables, label) {
# function to make PCOA, NMDS and box plots for the pathogens from a taxon_table
  # ~ taxon table - is a matrix or dataframe of the count of taxons per sample
  # ~ mapping_file - is the metadata file containing each sample description
  # ~ independent_variables - is a string vector or string specifying the 
  #the independent variables in the supplied mapping file
  # ~ label is a string that specifies the column in the mapping file for the points 
  #on the ordination plots e.g 'sample'
  
  #import your qiime formatted mapping file
  map <- read.table( file = mapping_file, sep='\t', head=T, row.names=1, comment = '')
  
  ##### making PCOA plots #######
  taxon_table_t <- t(taxon_table) #transpose the OTU table
  #find the overlap
  common.ids <- intersect(rownames(map), rownames(taxon_table_t ))
  # get just the overlapping samples
  taxon_table_t <- taxon_table_t[common.ids,]
  map <- map[common.ids,]
  #covert to relative abundance
  taxon_norm <- sweep(taxon_table_t, 1, rowSums(taxon_table_t),'/')
  #get the distance matrix between samples using the bray curtis distance metric 
  d.bray <- vegdist(x = taxon_norm , method = "bray")
  #calculate the principal coordinates
  pc.bray <- cmdscale(d = d.bray, k=2)
  colnames(pc.bray) <- c('PCoa1', 'PCoa2') #name the columns
  
  #preparing to make a biplot
  #First get a normalized version of the OTU table
  # where each OTU's relative abundances sums to 1
  otus.norm <- sweep(taxon_table_t,2,colSums(taxon_table_t),'/')
  # use matrix multiplication to calculate the weighted average of each taxon along each axis
  wa <- t(otus.norm) %*% pc.bray
  colnames(wa) <- c('wa1', 'wa2')
  wa_stat.df <- data.frame(taxon= rownames(wa), wa1=wa[,'wa1'], wa2 = wa[,'wa2'])
  
  # NMDS using Vegan package
  mds.bray <- metaMDS(taxon_norm)$points
  ##name the columns
  colnames(mds.bray) <- c('NMDS1', 'NMDS2')
  
  #get the independent variables
  independent_variables.df <- subset(x = map, select = independent_variables)
  
  #create a PCOA stat dataframe
  pca.stat.df <- data.frame(independent_variables.df,samples = rownames(pc.bray), PC1 = pc.bray[,'PCoa1'], PC2 = pc.bray[,'PCoa2'])
  #create NMDS sta data frame
  NMDS.stat.df <- data.frame(independent_variables.df,samples = rownames(mds.bray), NMDS1 = mds.bray[,'NMDS1'], NMDS2 = mds.bray[,'NMDS2'])
  #generate a statistics table for plotting and statistics
  stats_table <- cbind(independent_variables.df,taxon_norm)
  #plot the PCOA scatter plot
  for (variable in independent_variables) {
    p<- ggscatter(data = pca.stat.df, x = 'PCoa1', y = 'PCoa2', 
                  fill = variable, color = variable, size = 7,
                  label = label)
    
    #making a biplot - simply by adding the weighted average points to
    #the existing scatter plot
    print(ggscatter(data = wa_stat.df, x = 'wa1', y = 'wa2',
                    label = 'taxon', ggp = p))
 
    
    #plot the NMDS scatter plot -ggp adds the plot to an existing ggplot
    print(ggscatter(data = NMDS.stat.df, x = 'NMDS1', y = 'NMDS2', 
                    fill = variable, color = variable, size = 5,
                    label = label, ellipse = T))
    
    
  }
  final_tables <- list(wa_stat.df, pca.stat.df, NMDS.stat.df, stats_table)
  return(final_tables)
}

#function to stack your data from wide to long format
Stack_data <- function(Data,lab_columns,num_of_numericVariables,col_num_of_numericVariables) {
  # ~ data - is the dataframe to stack
  # ~ lab_colums - is a vector of numbers or strings specifying the columns 
  #containing the labels or independent variables
  
  # ~ num_of_numericVariables is a number specifying how many 
  #numeric / denpendent variables you want to stack 
  
  # ~ col_num_of_numericVariables  - is a vector of numbers or strings specifying the columns 
  #containing the numeric or dependent variables 
  
  #subsetting the non-numeric variable columns
  lab <- Data[,lab_columns]
  if (length(lab_columns) == 1){
    tab <- cbind(levels(lab))
    #bind the non-numeric columns row-wise 18 times
    for (i in 1:(num_of_numericVariables-1)){
      tab <- rbind(tab,cbind(levels(lab)))
      
    }
    
  }else{
    
    #initialize an empty dataframe
    tab <- data.frame()
    #bind the non-numeric columns row-wise 18 times
    for (i in 1:num_of_numericVariables){
      tab <- rbind(tab,lab)
    }
  }  
  #stack the numeric variable columns together
  j <- stack(Data[,col_num_of_numericVariables])
  #rename the column names
  colnames(j) <- c("Relative.abundance", "Pathogens")
  #re-arrange the columns
  k <- data.frame(Pathogens = j$Pathogens, Relative.abundance = j$Relative.abundance)
  #combine the stacked variable to form a new stacked dataframe
  stacked_data2 <- cbind(tab,k)
  colnames(stacked_data2) <- c(lab_columns,colnames(k))
  rownames(stacked_data2) <- 1:length(stacked_data2[,1])
  return(stacked_data2)
}

#function to generate stacked frequency and abundance tables
generate_abundance_tables <- function(independent_variables,variable,taxon_table) {
  #function to generate abundance table and frequency table of a specified independent
  #variable from a taxon table
  
  # ~ independent_variables - is a dataframe of independent variables or a factor vector of
  #an indenpendent variable
  # ~ variable  - is a string specifying the name of the independent variable that you
  #want to derive it's frequency and abundance tables
  # ~ taxon_table - is a taxon table generated by make_OTU_tables function
  
  #tanspose the taxon table
  taxon_table <- t(taxon_table)
  #find the overlap
  common.ids <- intersect(rownames(independent_variables), rownames(taxon_table))
  # get just the overlapping samples
  taxon_table <- taxon_table[common.ids,]
  independent_variables <- independent_variables[common.ids,]
  
  
  #check if the given independent_variables is a dataframe so as to decide how to generate the 
  #genreral frequency table 'j'
  if (is.data.frame(independent_variables) == FALSE){
    j<- data.frame(independent_variables, taxon_table)
  }else{
    j<- data.frame(independent_variables[,variable], taxon_table)
  }
  colnames(j)<- c(variable, colnames(taxon_table)) #rename the columns
  #create a list of lists of the same lenght as the levels of the independent variable
  subset_list <- vector(mode = "list", length = length(levels(j[,variable])))
  j_levels <- levels(j[,variable]) #get the levels of the independent variable
  
  #loop over every level of the independent variable and store the subsets in a list
  for (index in 1:length(levels(j[,variable]))){
    subset_list[[index]]<- j[j[,variable]== j_levels[index],]
  }
  names(subset_list) <- j_levels #name the list with the names of the levels of the independent variable
  
  #initialize an empty vector for the pathogen frequencies
  pathogen_freq<- c()
  
  #loop over every level of the independent variable while summing the columns
  #of the subset and row binding to the the previous subset
  for (level in j_levels){
    pathogen_freq <- rbind(pathogen_freq,colSums(subset_list[[level]][,-1]))
  }
  variable.df<- data.frame(j_levels,pathogen_freq) #generate the non-normalized frequency dataframe
  colnames(variable.df) <- c(variable, colnames(pathogen_freq)) #name the columns
  #generate the normalized dataframe of relative abundances
  variable.norm <- data.frame(j_levels,sweep(variable.df[,-1], 1, rowSums(variable.df[,-1]),'/'))
  colnames(variable.norm) <- c(variable, colnames(pathogen_freq))#name the columns
  #creating a vector that contains the  non-normalized  and normalized dataframes
  Data <- list(variable.df, variable.norm)
  #create a stacked_list of lists with its lenght equal to the lenght of data
  stacked_list <- vector(mode = "list", length = length(Data))
  
  #loop over every data, stack them and store them in the stacked_list
  for (i in 1:length(Data)){ 
    stacked_list[[i]] <- Stack_data(Data = Data[[i]], 
                                    lab_columns = variable, 
                                    num_of_numericVariables = length(colnames(pathogen_freq)), 
                                    col_num_of_numericVariables = 2:ncol(Data[[i]]))
    if (i == 1){
      colnames(stacked_list[[i]])<- c(variable, 'Pathogens', 'Frequency')
    }
  }
  names(stacked_list) <- c('frequency','relative.abundance')#name the stacked list
  return(stacked_list) #return the stacked list
}
generate_abundance_tables1 <- function(independent_Variables,variable,taxon_table,stats_table) {
  #function to generate abundance table and frequency table of a specified independent
  #variable from a taxon table
  
  # ~ independent_variables - is a dataframe of independent variables or a factor vector of
  #an indenpendent variable
  # ~ variable  - is a string specifying the name of the independent variable that you
  #want to derive it's frequency and abundance tables
  # ~ taxon_table - is a taxon table generated by make_OTU_tables function
  
  taxon_table <- t(taxon_table)
  #get the independent variables
  independent_variables<- stats_table[,1:length(independent_Variables)]
  
  #test if the independent variables variable is a dataframe if not change it to a dataframe
  if (is.data.frame(independent_variables) == FALSE){
    independent_variables <- data.frame(independent_variables)
    rownames(independent_variables) <- rownames(stats_table)
    common.ids <- intersect(rownames(independent_variables), rownames(taxon_table))
    independent_variables <- independent_variables[common.ids,]
    independent_variables <- data.frame(independent_variables)
    rownames(independent_variables) <- common.ids
    colnames(independent_variables) <- variable
  }else{

 
  #find the overlap
  common.ids <- intersect(rownames(independent_variables), rownames(taxon_table))
  # get just the overlapping samples
  taxon_table <- taxon_table[common.ids,]
  independent_variables <- independent_variables[common.ids,]
  }
  
 
    j<- data.frame(independent_variables[,variable], taxon_table)

  colnames(j)<- c(variable, colnames(taxon_table)) #rename the columns
  #create a list of lists ofthe same lenght as the levels of independent variable
  subset_list <- vector(mode = "list", length = length(levels(j[,variable])))
  j_levels <- levels(j[,variable]) #get the levels of the independent variable
  
  #loop over every level of the independent variable and store the subsets in a list
  for (index in 1:length(levels(j[,variable]))){
    subset_list[[index]]<- j[j[,variable]== j_levels[index],]
  }
  names(subset_list) <- j_levels #name the list with the names of the levels of the independent variable
  
  #initialize an empty vector for the pathogens frequenciess
  pathogen_freq<- c()
  
  #loop over every level of the independent variable while summing the columns
  #of the subset while row binding to the the prior subset
  for (level in j_levels){
    pathogen_freq <- rbind(pathogen_freq,colSums(subset_list[[level]][,-1]))
  }
  variable.df<- data.frame(j_levels,pathogen_freq) #generate the non-normalized frequency dataframe
  colnames(variable.df) <- c(variable, colnames(pathogen_freq)) #name the columns
  #generate the normalized dataframe of relative abundances
  variable.norm <- data.frame(j_levels,sweep(variable.df[,-1], 1, rowSums(variable.df[,-1]),'/'))
  colnames(variable.norm) <- c(variable, colnames(pathogen_freq))#name the columns
  #creating a vector that contains the  non-normalized  and normalized dataframes
  Data <- list(variable.df, variable.norm)
  #create a stacked_list of lists with its lenght equal to the lenght of data
  stacked_list <- vector(mode = "list", length = length(Data))
  
  #loop over every data, stack them and store them in the stacked_list
  for (i in 1:length(Data)){ 
    stacked_list[[i]] <- Stack_data(Data = Data[[i]], 
                                    lab_columns = variable, 
                                    num_of_numericVariables = length(colnames(pathogen_freq)), 
                                    col_num_of_numericVariables = 2:ncol(Data[[i]]))
    if (i == 1){
      colnames(stacked_list[[i]])<- c(variable, 'Pathogens', 'Frequency')
    }
  }
  names(stacked_list) <- c('frequency','relative.abundance')#name the stacked list
  result<- list(stacked_list,variable.df,variable.norm)
  names(result)<- c('stacked_list', 'raw_table','relative_abundance_table')
  return(result) #return the stacked list
}

