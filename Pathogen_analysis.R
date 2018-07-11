#pathogen analysis script
#Author: Olabiyi Aderemi Obayomi
#E-mail: obadbotanist@yahoo.com
#created Feb 2018

#remember to set your working directory to directory containing your results together and your mapping file
######################## set all necessary variables   ####################################

#example of how to set the necessary variables
file_name <- "wash_bac_pathogen_analysis"
independent_Variables <- c('WaterType', 'Plastic', 'Treatment', 'IrrigationType') 
mapping_file <- 'fruit_wash_pathogens_mapping.txt'
# true or false if plots should be made
#I recommend sett it to false when you don't detect alot of pathogens in your samples
#for example when analysing for protists in the water and bacteria in the fruit puree
#but set to true when you analysing bacteria in the soil
make_plots <- FALSE 
Xlab <- c('Water Type','Plastic cover','Treatment', 'Irrigation Type')

#load all the library and functions needed for the anlysis
source('Pathogen_analysis_functions.R')

# sample of character vector of blast files in the current working directory
#change to the directory with all the results together
setwd('results_together/')
#get the pathogen files arranged from the lowest to the highest number
pathogen_files <- list.files() #get all the files in the current directory
pathogen_files<- gsub(pattern = "__", 
                      replacement = ".",
                      x = mixedsort(gsub(pattern = "\\.", replacement = "__", x = pathogen_files)))

#get the Treatment or sample names
#character vector of sample or treatments to be analysed as they appear in your mapping file
Treatments <- c() #intialize an empty treament vector
#get a split file name by substititing "final_" with "" and spliting the file name by "_" 
split_file_name<- strsplit(x = gsub(pattern = "final_", replacement = "", x = pathogen_files),
                           split = "_" )
#loop over the split_file_name list of lists getting the first element of each list
#and setting the respective treatment value to the rightful index in the vector
#Treatments
for (index in 1:length(pathogen_files)) {
  Treatments[index] = split_file_name[[index]][1]
}

#creating a list of list which will contain the outputs for the pathogen analysis
Result_together <- vector(mode = "list",length = length(Treatments))

#####################################################################################

#perform the pathogen analysis for each treatment or sample
#turn off annoying warning that could be ignored
options(warn = -1)
#loop to perform the analysis for the number of treatments required
for (index in 1:length(Treatments)) {
  print(pathogen_files[index])
  Result_together[[index]] <- pathogen_analysis(treatment = Treatments[index], 
                                                pathogen_file = pathogen_files[index])
}
names(Result_together) <- Treatments
#change to the parent directory
setwd('../')
#turn warnings back on
options(warn = 0)
##########Get the frequency tables for all the taxon levels together##################
#in a list called all_frequency_tables
#sample taxonomic level vector
taxons <- c('genera','species', 'strains')
all_frequency_tables <- vector(mode = "list",length = length(taxons))
for (index in 1:length(taxons)) {
  
  taxon_frequency_table <- combine_freq_tables(Treatments = Treatments,
                                               Result_together = Result_together,
                                               frequency_table_name = taxons[index])
  all_frequency_tables[[index]] <- taxon_frequency_table
}
names(all_frequency_tables) <- taxons

##################################make taxon tables #####################
#make raw and normalized taxon/otu tables for all the taxon levels and combine them
combined_taxon_tables <- make_taxon_tables(samples = Treatments, 
                                           taxons = taxons, 
                                           all_frequency_tables = all_frequency_tables)


################################get the number of pathogens detected per treament group #######
#get the independent variables dataframe
independent_variables.df <- get_independent_variables_df(mapping_file = mapping_file, 
                                                         independent_variables = independent_Variables)

#intialize a list to store the results of the pathogens detected per treatment groups
detect <- vector(mode = 'list', length = length(independent_Variables))
#get the taxon table
taxon_table <- combined_taxon_tables$raw_tables$species
#get the frequency of detection for each pathogen per treatment group
for (index in 1:length(independent_Variables)) {
  
detect[[index]]<- presence_absence(taxon_table = taxon_table, 
                          independent_variables.df = independent_variables.df, 
                          variable = independent_Variables[index])
}
names(detect) <- independent_Variables
#write the results to an excel file per independent variavle
for (i in names(detect)){
write.xlsx(x = detect[[i]], 
           file = sprintf(fmt = '%s.xlsx', file_name), 
           sheetName = i, 
           col.names = T, 
           row.names = F,
           append = TRUE)
}


if (make_plots == TRUE){
  

############################### Plotting ##############################################
#open and pdf device and set it's width and height
pdf(file =sprintf('%s.pdf',file_name),width = 14, height = 10.5)

#generate a list of lists that will contain ordination tables for each taxon level
ord_tables <- vector(mode = "list",length = length(taxons))
stats_table <- vector(mode = "list",length = length(taxons))
stacked_stats_table <- vector(mode = "list",length = length(taxons))
abund_table_norm <- vector(mode = "list",length = length(taxons))
for (taxon in 1:length(taxons)){
  #make ordanition plots - PCOA and NMDS and return dataframes with the principal components
  ord_tables[[taxon]] <-ordination_plots(taxon_table = combined_taxon_tables[["normalized_tables"]][[taxon]],
                                         mapping_file = mapping_file, 
                                         independent_variables = independent_Variables,
                                         label =  'samples')
  
  #get the relative abundance table per sample with the independent variables attached
  stats_table[[taxon]] <- ord_tables[[taxon]][[4]]
  num_of_numericVariables <- ncol(stats_table[[taxon]]) - length(independent_Variables)
  col_num_of_numericVariables <- (length(independent_Variables)+1):ncol(stats_table[[taxon]])
  stacked_stats_table[[taxon]]<- Stack_data(Data = stats_table[[taxon]], 
                                            lab_columns = independent_Variables,
                                            num_of_numericVariables = num_of_numericVariables , 
                                            col_num_of_numericVariables = col_num_of_numericVariables)
  
  for (i in 1:length(independent_Variables)){ 
    #making relative abundance boxplots faceted by pathogens with each sample info
    gg <- draw_pathogen_plots(Data = stacked_stats_table[[taxon]], 
                              main = 'Difference between Pathogens',
                              x = independent_Variables[i], 
                              y = 'Relative.abundance', 
                              Xlab = Xlab[i], 
                              YLab = 'Relative abundance', 
                              plotType = 'box', 
                              group = 'Pathogens')
    invisible(print(ggpar(gg, tickslab = FALSE, ticks = FALSE, xlab = FALSE)))
    #make long relative abundance box plots
    long_gg<- ggboxplot(data = stacked_stats_table[[taxon]], x = 'Pathogens', 
                        y = 'Relative.abundance', 
                        #fill = independent_Variables[i],
                        color = independent_Variables[i],
                        add = 'mean_se')
    invisible(print(ggpar(long_gg,x.text.angle = 90, 
          rotate = TRUE, 
          font.x = c(14, "bold", "black"),
          font.y = c(14, "bold", "black"),xlab = taxons[taxon], ylab = 'Relative abundance')))
    
    #generate stacked frequency and abundance tables of an independent variable 
    #first check if the independent variable is not a dataframe
    independent_variables<- stats_table[[1]][,1:length(independent_Variables)]
    if (is.data.frame(independent_variables) == FALSE){
      independent_variables<- as.data.frame(independent_variables)
      rownames(independent_variables) <- rownames(stats_table[[1]])
      colnames(independent_variables) <- independent_Variables
    }
    #generate stacked frequency and abundance tables of an independent variable
    abund_table_norm[[taxon]][[i]] <- generate_abundance_tables(independent_variables = independent_variables, 
                                                                variable = independent_Variables[i], 
                                                                taxon_table = combined_taxon_tables[["raw_tables"]][[taxon]])$relative.abundance
    ##make relative abunbance bar charts
    draw_pathogen_plots(Data = abund_table_norm[[taxon]][[i]], 
                        main = 'Pathogens relative abundance', 
                        x = independent_Variables[i], 
                        y = 'Relative.abundance', 
                        Xlab = Xlab[i], 
                        YLab = 'Relative abundance', 
                        plotType = 'bar', 
                        group = 'Pathogens')
    
  }
  
}
#nameing the lists to correspond to each taxon level
names(ord_tables) <- taxons
names(stats_table) <- taxons
names(stacked_stats_table) <- taxons
names(abund_table_norm) <- taxons

########################## make pathogen Heat map #########################################
#get the raw table that belongs to
x <- generate_abundance_tables1(independent_Variables = independent_Variables,
                                variable = independent_Variables[1],
                                taxon_table = combined_taxon_tables[["raw_tables"]][[2]],
                                stats_table = species_stats_table)$raw_table

#rename the rownames of othe OTU tables 
rownames(x) <- x[,1]
#remove the treatment coolumn
x <- x[,-1]
#normalize the OTU table
x_norm <- normalize_otu_table(OTU_table = x, method = "logUQ")
#calcultate the relative abundance of teh normalized OTU table
x_norm.abund <- as.matrix(sweep(x = x, MARGIN = 1, rowSums(x_norm),'/'))
x_norm.abund <- t(x_norm.abund)
#get the index of the ordered pathogen names
id<- order(rownames(x_norm.abund))
#get their names based on the idex
id<- rownames(x_norm.abund)[id]
#re-arrnge the names alphabetically
x_norm.abund <- x_norm.abund[id, ]

#reorder the columns of a matrix
#col.order <- c("UnIrrigated_clay", "PW+Clay", "TWW+Clay", "UnIrrigated_loam", "PW+Loam", "TWW+Loam", "UnIrrigated_loamSand", "PW+Loamy_sand", "TWW+Loamy_sand")
#col.order <- c("UnIrrigated_clay", "PW-P", "PW+P", "TWW-P", "TWW+P")
col.order <- c("PW", "TWW","non-irrigated", "PW-P", "PW+P", "TWW-P", "TWW+P", "Blank",
               "PW+SD-P", "PW+SD+P", "TWW+SD-P", "TWW+SD+P", "TWW+SSD-P", "TWW+SSD+P")
ordered_x_norm.abund <- x_norm.abund[,col.order]
#rename the colums to represent the treatments the way you want them if need be
#colnames(ordered_x_norm.abund) <- c("non-irrigated","PW-P","PW+P","TWW-P","TWW+P")


#make heatmap of normalized otu table clustered which hclust with bray curtis distance
heatMap <- Heatmap(ordered_x_norm.abund, 
        cluster_columns = FALSE, 
        cluster_rows = FALSE, 
        column_title = 'Treatment', 
        row_names_side = 'left', 
        column_title_side = 'bottom', 
        #clustering_distance_columns = function(x) vegdist(x), 
        row_names_gp = gpar(fontface="bold.italic",cex=0.6), 
        column_names_gp = gpar(fontface="bold"), 
        heatmap_legend_param = list(title= 'Relative abundance'))

#save the heatmap
save(heatMap, file = '../raw_melon_heatmap.RData')

#get the species taxon table
species_otu_table <- combined_taxon_tables[[1]][[2]]

#write the species taxon table to an excel file
write.xlsx(x = t(species_otu_table), file = '../melon_species_table.xlsx', sheetName = 'melon', col.names = T, row.names = T, append = T)

############################################################################################

#turn-off all open devices
dev.off()

# species_relative<- sweep(species_otu_table, 1, rowSums(species_otu_table),'/')
# #re-order factor levels
# mtcars$cyl2 <- factor(mtcars$cyl, levels = c("6","4","8"))
# #re-order decrete scale in ggplot2
# ggplot(mtcars, aes(factor(cyl))) + 
#   geom_bar() + 
#   scale_x_discrete(limits=c(8,4,6))#this is where we re-order
# 
# #get the frequency for each pathogen detected
# species_freq_table<- all_frequency_tables[[2]][1:3]

###################plotting multiple complex heatmas on a page###########################
#1.trying to use the mutipanelFigure package for multiple plots
library('multipanelfigure') #load the package
#create an empty figure
figure <- multi_panel_figure(
width = c(30,40,60),
height = c(40,60,60,60),
panel_label_type = "upper-roman") 

#file the figure with your plots
(figure %<>% fill_panel(
water_heatmap,
row = 1, column = 1))
(figure %<>% fill_panel(
barrier_heatmap,
row = 1, column = 2))
(figure %<>% fill_panel(
wash_heatmap,
row = 1, column = 3))
figure



#2. first convert to grob object the use ggarrage to plot and ararnge the heatmaps
#convert a complex heatmap object to a grob object
water_gb <- grid.grabExpr(draw(water_heatmap))
is.grob(water_gb)#confirm that is a grob object
#test to see it works
grid.draw(water_gb)
#barrier experiment
barrier_gb<- grid.grabExpr(draw(barrier_heatmap))
#soilType experiment
soilType_gb <- grid.grabExpr(draw(soilType_heatmap))
#fruit wash
wash_gb <- grid.grabExpr(draw(wash_heatmap))
#melon puree
melon_gb <- grid.grabExpr(draw(melon_heatmap))

#make multiple plots
#barrier experiment plot
ggarrange(water_gb,barrier_gb,wash_gb, nrow = 1, ncol = 3)

#soiltype experiment plots
ggarrange(water_gb,soilType_gb, nrow = 1, ncol = 2)



#estimating diversity matrics such as chao1 
pool <- poolaccum(t(species_taxon_table))
summary(pool, display = "chao")
qt_model<- estimateR(t(species_taxon_table))
}