
##### Run this script only after you ran pathogen analysis with the Pathogen_analysis.R script ######

#load the heatmap libraries
library(ComplexHeatmap)
library(circlize)

#get the independent variables dataframe
independent_Variables.df <- subset(x = species_stats_table,select = 'Treatment')

#transpose the raw species taxon table such that samples names are now rownames
species_taxon_table_t <- species_taxon_table
#make sure the rownames of the independent variables dataframe are the same as that of the species table
x <- species_taxon_table_t[rownames(independent_Variables.df),]
#column bind the independent variable dataframe with the raw taxon table
raw_table.df <- cbind(independent_Variables.df,x)
#calculate the sum of sequences detected per treatment 
treatment_raw_table.df <- aggregate(formula = .~Treatment, data = raw_table.df, FUN = sum)
#re-assign the rownames such that they represent each treatment and remove the first redundant column
rownames(treatment_raw_table.df) <- treatment_raw_table.df[,1]
#remove the first redundant column
treatment_raw_table.df<- treatment_raw_table.df[,-1]



#get the pathogen names
pathogens<-colnames(treatment_raw_table.df)
#convert to relative abundance table
abuns_treatment.df <- apply(X = treatment_raw_table.df, MARGIN = 1, FUN = function(x) {x/sum(x)})

#replace NAN with zero
library('gtools')
abuns_treatment.df <- na.replace(abuns_treatment.df,0)

id<- order(rownames(abuns_treatment.df))
#get their names based on the idex
id<- rownames(abuns_treatment.df)[id]
#re-arrnge the names alphabetically
abuns_treatment.df <- abuns_treatment.df[id, ]
#get the column order desired

#example of a vector of column order
col.order <- c("PW", "TWW","UnIrrigated_clay", "PW-P", "TWW-P", "PW+P",
               "TWW+P","Blank","PW+SD-P",
               "TWW+SD-P", "PW+SD+P","TWW+SD+P","C_PW",
               "C_TWW", "M_Blank", "M_PW", "M_TWW")



abuns_treatment.df<- abuns_treatment.df[,col.order]
#remove low abundance pathogens
abund_pathogens.df <- remove_low_abund_pathogens(t(abuns_treatment.df))

#get the treatment names
treatments <- colnames(abuns_treatment.df)

#re-arrange the columns
abuns_treatment.df <- abuns_treatment.df[id, ]

f6 <- colorRamp2(seq(0, max(abuns_treatment.df), length = 2), c("white", "red"))

#open a pdf device
pdf(file = 'heatmaps.pdf', onefile = TRUE, width = 14, height = 8 )

#make heatmap (full table)
invisible(print(heatMap5 <- Heatmap(abuns_treatment.df,
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    row_names_side = 'left',
                    column_title_side = 'bottom',
                    row_names_gp = gpar(fontface="bold.italic",cex=0.8),
                    column_names_gp = gpar(fontface="bold"),
                    col = f6,
                    rect_gp = gpar(col="grey",lwd=1.5),
                    heatmap_legend_param = list(title= 'Relative abundance', 
                                                labels_gp=gpar(fontsize=12, fontface="bold"),
                                                title_gp = gpar(fontsize = 14, fontface = "bold")))))

#make heatmap of only the abundant pathogens
invisible(print(heatMap6 <- Heatmap(t(abund_pathogens.df),
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    row_names_side = 'left',
                    column_title_side = 'bottom',
                    row_names_gp = gpar(fontface="bold.italic",cex=1),
                    column_names_gp = gpar(fontface="bold"),
                    col = f6,
                    rect_gp = gpar(col="grey",lwd=1.5),
                    heatmap_legend_param = list(title= 'Relative abundance', 
                                                labels_gp=gpar(fontsize=12, fontface="bold"),
                                                title_gp = gpar(fontsize = 14, fontface = "bold")))))

dev.off()

