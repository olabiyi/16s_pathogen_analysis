# Statistical analyses
#load ggpubr package
library('ggpubr')

# Get relative abundance species table
abund_species_table <- sweep(x = species_table,
                             MARGIN = 2, 
                             STATS = colSums(species_table), 
                             FUN = "/")

abund_species_table_t <- t(abund_species_table)

#get the common ids between the species table and 
common.ids <- intersect(rownames(ind_var),rownames(species_table_t))

# Get the common ids for the independent variables
ind_var <- ind_var[common.ids,, drop=FALSE]

abund_species_table_t <- abund_species_table_t[common.ids, , drop=FALSE]
# Make rarefaction curves
rarecurve(species_table_t)

# Create a dataframe with the indentification for each sample attached to the species table
raw_x <- cbind(ind_var,species_table_t)
abund_x <- cbind(ind_var,abund_species_table_t)

# Test for normality
shapiro.test(x_diversity$S.obs)
shapiro.test(x_diversity$S.chao1)
shapiro.test(x_diversity$S.ACE)
shapiro.test(x_diversity$shannon)
shapiro.test(x_diversity$simpson)

# Test if there are signifant differnces between treatments based on the diversity matrices
kruskal.test(formula = S.obs ~ Treatment, data = x_diversity)
kruskal.test(formula = S.chao1 ~ Treatment, data = x_diversity)
kruskal.test(formula = S.ACE ~ Treatment, data = x_diversity)
kruskal.test(formula = shannon ~ Treatment, data = x_diversity)
kruskal.test(formula = simpson ~ Treatment, data = x_diversity)


#find where the differences lie
#compare treatments
compare_means(formula = simpson ~ Treatment, data = x_diversity)

#compare soil types
compare_means(formula = simpson ~ Plastic, data = x_diversity)

#compare water types
compare_means(formula = simpson ~ WaterType, data = x_diversity)


# Perform adonis test
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


write.xlsx(x = result, 
           file = 'pathogen_diversity.xlsx', 
           sheetName = 'differential_abundance', 
           col.names = T, 
           row.names = F,
           append = T)