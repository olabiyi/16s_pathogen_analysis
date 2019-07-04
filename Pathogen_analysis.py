# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 21:24:15 2018

@author: Olabiyi Aderemi Obayomi
@email: obadbotanist@yahoo.com


"""

import os
import pandas as pd
import seaborn as sns
from skbio.stats.ordination import pcoa
import pathogenAnalysisFunctions as p
from plotnine import *
import warnings
import random


########### Set the required input variables (an example) 
file_name = "barrier_bacteria" 
MAP = "barrier_bacteria_mapping.txt"
independent_variables = ["WaterType","Plastic","Treatment","Season"]
sample_column = "sample" # your mapping file must have a column called sample
working_directory = (r"C:/Users/user/Documents/programming/Python/pathogen_analysis/barrier/results_together/")
group="Treatment"
# PCOA shape and color labels
color_var='WaterType'
shape_var='Plastic'
color_lab= 'Water Type',
shape_lab = 'Plastic'
point_colors=['blue','purple']
# Set color map for the heatmaps
cmap= 'hot' #["white","pink","red"] #


# Define a publication ready theme for making ggplots
publication_format = (theme_bw() +
  theme(panel_grid = element_blank()) +
  theme(axis_ticks_direction='in',
  axis_ticks_length=3.5,
  axis_text_x=element_text(margin={'t':0.1,'r':0,'b':0,'l':0,'units':"in"}),
  axis_text_y=element_text(margin={'t':0,'r':0.1,'b':0,'l':0,'units':"in"}), 
  axis_title = element_text(size = 18, style='italic', weight='bold', color = 'black'), 
  axis_text = element_text(size = 16, style ='italic', weight='bold', color = 'black'),
  legend_position = 'right', 
  legend_text = element_text(size = 14, weight ='bold', color = 'black'),
  legend_title = element_text(size = 15, weight ='bold', color = 'black'),
  strip_text =  element_text(size = 14, weight ='bold', color = 'black')))

# Generate random colors and add them to the custom palette below when you
#don't have enough colors for making your stack bar plot
random.seed(1)
#generate more colors and test if it works 
_,more_colors= p.generate_colors(n= 8)

custom_palette = ["#a6cee3","#1f78b4","#b2df8A","#33a02c","#fb9a99","#e31a1c","#fdbf6f", "#ff7f00",
                  "#6A3D9A","#FF00FFFF","#B15928","#000000","#FFC0CBFF","#8B864EFF","#F0027F",
                  "#666666","#1B9E77", "#E6AB02","#A6761D","#FFFF00FF","#FFFF99","#00FFFFFF",
                  "#B2182B","#FDDBC7","#D1E5F0","#CC0033","#FF00CC","#330033","#999933","#FF9933",
                  "#FFFAFAFF",'#988ea3','#b8aec3','#d8cee3','#382e43','#584e63','#786e83','#5726ac',
                  '#7746cc','#9766ec','#d7a62c', '#f7c64c','#17e66c','#49b926','#69d946','#89f966',
                  '#a91986','#c939a6','#e959c6', '#42f8e3', '#62183','#823823', '#a25843', '#c27863',
                  '#e29883', '#2b8a3', '#22d8c3']





# Turn off annoying warnings
warnings.filterwarnings('ignore')

######## Run pathogen analysis and combine the results for all samples
# then return the combined results and the sample names
os.chdir(working_directory)
files = os.listdir()
result_together,samples = p.combined_results(files)

# Make pathogen taxon table 
taxon_table = p.make_taxon_table(result_together, samples)

# Change to the parent directory
os.chdir("../")

# Write-out the pathogen table to excel 
writer = pd.ExcelWriter("{}_table.xlsx".format(file_name))
taxon_table.to_excel(writer,sheet_name="taxon_table") 


# Get the independent variables dataframe
mapping_file = pd.read_csv(MAP, sep="\t") 
independent_variables.append(sample_column)
independent_variable_df = mapping_file[independent_variables].set_index(sample_column)
independent_variables.remove(sample_column)

# Get the transposed pathogen dataframe and independent variables with common ids
# by joining the independent variable dataframe with the pathogen dataframe
# based on samples common to both dataframes
taxon_table_t = taxon_table.T

independent_variable_df, taxon_table_t, abund_df = p.join_tables(table1=independent_variable_df,
                                                               table2=taxon_table_t)

# Get the number of pathogens detected per treatment and 
# write the resulting table to excel
detect_df = p.presence_abscence(taxon_table_t=taxon_table_t,
                                independent_variables=independent_variable_df)
detect_df.to_excel(writer,sheet_name="detects")

# Remove samples that pathogens were no detected in
taxon_table_t = p.remove_zero_samples(taxon_table_t=taxon_table_t)
# Estimate diversity and write the diversity table to excel
diversity_table = p.estimate_diversity(taxon_table_t= taxon_table_t)
diversity_table.to_excel(writer,sheet_name="diversity")

# Save the excel writer object for the open excel file
writer.save()


# Convert the absolute pathogen taxon table to a relative abundance table
abund_table = taxon_table_t.apply(func=lambda x: x/sum(x),axis=1)
# Drop samples with NA having zero counts
abund_table = abund_table.dropna()



# Calculate bray curtis dissimilarity matrix
brayMat = p.table_to_distances(abund_table.T,p.bray_curtis_distance)
# Perform principal co-ordinate analysis
bc_pc= pcoa(brayMat)

# Extract the first two principle components and combine
# them with the independent variables dataframe
pc_df=bc_pc.samples[["PC1","PC2"]]

independent_variable_df,_,pc_df = p.join_tables(table1=independent_variable_df,
                                                               table2=pc_df)

proportion_explained = bc_pc.proportion_explained[0:2] * 100

xlab = "PCO1({:.2f}%)".format(proportion_explained["PC1"])
ylab = "PCO2({:.2f}%)".format(proportion_explained["PC2"])

# Make the PCoA plot and save it as a png file 
gg = (ggplot(data=pc_df) + 
 geom_point(aes(x='PC1', y='PC2', color=color_var, shape=shape_var), size=6) + 
 labs(x=xlab,y=ylab,color= color_lab, shape= shape_lab) + 
 scale_color_manual(values=point_colors) + publication_format)
gg

ggsave(plot=gg, filename='{}_pcoa.png'.format(file_name),
       width=14, height=8, dpi=600, device='png')



## Prepare data for plotting relative abundances based on treatment category ##########

# Aggregate by group and then stack the data
treat_abund_table = p.get_treatment_abund_table(independent_variable_df, taxon_table_t, group=group)


treat_rare_abund_table=p.remove_low_abundance_pathogens(abun_table=treat_abund_table,
                                               group_low_abund=True,
                                               threshold=0.02)
treat_rare_abund_table[group] = treat_rare_abund_table.index

stacked_treat_abund_table = treat_rare_abund_table.melt(id_vars=[group], 
                                           value_vars=list(treat_rare_abund_table.columns.drop(group)), 
                                           value_name="Relative_abundance", 
                                           var_name="Pathogens")


# Make stacked barplot

gg = (ggplot(data= stacked_treat_abund_table) +
geom_col(aes(x=group,y="Relative_abundance", fill="Pathogens")) + 
 labs(y="Relative abundance(%)") +
 scale_fill_manual(values=custom_palette) +
 publication_format)
gg
ggsave(plot=gg, filename='{}_bar.png'.format(file_name),
       width=14, height=8, dpi=600, device='png')


# Make box plots

gg =(ggplot(data= stacked_treat_abund_table) + 
    geom_boxplot(aes(x="Pathogens",y="Relative_abundance", color=group)) +
    labs(y="Relative abundance(%)") +
    publication_format + theme(axis_text_x=element_text(angle=45)))
gg
ggsave(plot=gg, filename='{}_boxplot.png'.format(file_name),
       width=14, height=8, dpi=600, device='png')

# Prepare to make heatmaps
# Optionally remove pathogens that are not atleast 5% abundant in a sample
abund_pathogens_table = p.remove_low_abundance_pathogens(
                                              abun_table=treat_abund_table,
                                              group_low_abund=False)

# Sort the columns/pathogens in alphabetical order
col_names = abund_pathogens_table.columns.tolist()
col_names.sort()
# Re-arrange the columns/pathogens in alphabetical order
abund_pathogens_table = abund_pathogens_table[col_names]

##### Make heatmaps

# Heatmap for all the pathogens
col_names = treat_abund_table.columns.tolist()
# Sort the columns/pathogens in alphabetical order
col_names.sort()
# Re-arrange the columns/pathogens in alphabetical order
treat_abund_table = treat_abund_table[col_names]
fig,r = plt.subplots(figsize=(14,7))
sns.heatmap(treat_abund_table.T, cmap=cmap,\
              linecolor="black", linewidths=0.01, yticklabels=True, ax=r)
plt.show(r)
# Save the heatmap to a png file
r.figure.savefig("all_pathogens_heatmap.png", dpi=600)

# Heatmap for only the abundant pathogens
fig,a = plt.subplots(figsize=(14,7))
sns.heatmap(abund_pathogens_table.T, cmap=cmap,\
              linecolor="black", linewidths=0.01,
            ax=a, yticklabels=True)
plt.show(a)
# Save the heatmap to a png file
a.figure.savefig("abundant_pathogens_heatmap.png", dpi=600)

# Restore warnings
warnings.filterwarnings('default')