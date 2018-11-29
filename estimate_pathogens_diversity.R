
##### Run this script only after you ran pathogen analysis with the Pathogen_analysis.R script ######


####################### how to perform diversity estimation and adonis test ################
#the vegan package
#if the package isn't already install it then load it
if (!require('vegan')){ install.packages('vegan'); library('vegan')}
#library('vegan')

#get the pathogens species table
species_table <- combined_taxon_tables$raw_tables$species

#drop samples with no detects
species_table <-  species_table[,-c(which(colSums(species_table) == 0))]

#get relative abundance species table
abund_species_table <- sweep(x = species_table, 
                             MARGIN = 2, 
                             STATS = colSums(species_table), 
                             FUN = "/")


#transpose the species table
species_table_t <- t(species_table)
abund_species_table_t <- t(abund_species_table)
#read in the mapping file
map <- read.table( file = './mapping.txt',
                   sep='\t', head=T, row.names=1, comment = '')

#get the independent variables to be tested
ind_var <- subset(x = map, select = c('Treatment', 'WaterType', 'Plastic'))

#get the common ids between the species table and 
common.ids <- intersect(rownames(ind_var),rownames(species_table_t))

#get the common ids for the independent variables
ind_var <- ind_var[common.ids,]
species_table_t <- species_table_t[common.ids,]
abund_species_table_t <- abund_species_table_t[common.ids,]
#make rarefaction curves
rarecurve(species_table_t)

#estimate abundance based richness - observed species, chao1 and ACE
richness <- estimateR(species_table_t)
#remove the standard error columns
richness <- richness[-c(3,5),]


#estimate diversity
shanon_index <- diversity(species_table_t) #shannon
simpson_index <- diversity(species_table_t, index = 'simpson') #simpson

#create a dataframe with the samples treatment identification columns
raw_x <- cbind(ind_var,species_table_t)
abund_x <- cbind(ind_var,abund_species_table_t)
#create a diversity dataframe
x_diversity <- cbind(ind_var,t(richness), shannon = shanon_index, simpson = simpson_index)

#test for normality
shapiro.test(x_diversity$S.obs)
shapiro.test(x_diversity$S.chao1)
shapiro.test(x_diversity$S.ACE)
shapiro.test(x_diversity$shannon)
shapiro.test(x_diversity$simpson)

#test if there are signifant differnces between treatments based on the diversity matrices
kruskal.test(formula = S.obs ~ Treatment, data = x_diversity)
kruskal.test(formula = S.chao1 ~ Treatment, data = x_diversity)
kruskal.test(formula = S.ACE ~ Treatment, data = x_diversity)
kruskal.test(formula = shannon ~ Treatment, data = x_diversity)
kruskal.test(formula = simpson ~ Treatment, data = x_diversity)

#load ggpubr package
library('ggpubr')
#find where the differences lie
#compare treatments
compare_means(formula = simpson ~ Treatment, data = x_diversity)

#compare soil types
compare_means(formula = simpson ~ Plastic, data = x_diversity)

#compare water types
compare_means(formula = simpson ~ WaterType, data = x_diversity)


#perform adonis test
#here on treatments
adonis2(formula = species_table_t ~ Treatment, data = ind_var)


#initialize i the count variable for the while loop
i <- 1
treatment <- c()
pathogens <- colnames(abund_x)[4:ncol(abund_x)]

while (i <= length(pathogens)) {

  treatment [i] <- kruskal.test(abund_x[,pathogens[i]] ~ abund_x[,'Treatment'])$p.value
  i = i + 1
}
#correct the p-values for multiple hypothesis testing
treatment_fdr <- p.adjust(treatment, method = 'fdr')


result <- data.frame(pathogens = pathogens,  
                     treatment = treatment, 
                     treatment_fdr = treatment_fdr)


#get the means and SE of the diversity matrices
#mean
mean_diversity <- aggregate(x = x_diversity[,4:ncol(x_diversity)], 
                            by = list(x_diversity$Treatment), 
                            FUN = mean)
#standard error
SE_diversity <- aggregate(x = x_diversity[,4:ncol(x_diversity)], 
                          by = list(x_diversity$Treatment), 
                          FUN = function(x) {sd(x)/sqrt(length(x))})

diversity_result <- data.frame(Treatment = mean_diversity$Group.1, 
                               obseved_species = mean_diversity$S.obs, 
                               SE_observed = SE_diversity$S.obs, 
                               Chao1 = mean_diversity$S.chao1,
                               SE_Chao1 = SE_diversity$S.chao1,
                               shannon = mean_diversity$shannon,
                               SE_shannon = SE_diversity$shannon,
                               simpson = mean_diversity$simpson,
                               SE_simpson = SE_diversity$simpson)

library('xlsx')

write.xlsx(x = diversity_result, 
           file = 'pathogen_diversity.xlsx', 
           sheetName = 'diversity_estimation', 
           col.names = T, 
           row.names = F)

write.xlsx(x = result, 
           file = 'pathogen_diversity.xlsx', 
           sheetName = 'differential_abundance', 
           col.names = T, 
           row.names = F,
           append = T)