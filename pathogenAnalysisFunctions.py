# -*- coding: utf-8 -*-
# """
# Created on Sat Nov 24 21:12:23 2018

# @author: Biyi

# In order to install and use modules like ecopy on Windows you'll need to 
# download and install visual studio 2015

# 1. From this link https://stackoverflow.com/questions/44290672/
#    how-to-download-visual-studio-community-edition-2015-not-2017

# Downdload the web installer or iso for windows. 
# After this must have installed (it takes a few hours).

# 2. Then run the command 
# a. 'pip install msgpack' in the Anaconda prompt
# b. 'pip install ecopy' in the Anaconda prompt

# To get a work around version of scikit bio for Windows follow the instructions
#  below:

# 1. Install Windows Visual Studio 2015 Community Edition with C++ SDK enabled to
#    meet the requirements of the compiler for installing the package(Scikit-Bio) 
#    from the link https://go.microsoft.com/fwlink/?LinkId=532606&clcid=0x409 
#    (web installer)

# 2. Download the latest source from the Scikit-bio Github repository  
#    https://github.com/biocore/scikit-bio/archive/master.zip.

# 3. Using a tool like 7-zip, unpack it into your python packages directory
#    (C:\Users\user\Anaconda3\pkgs).

# 3. Open a command line client (cmd) and navigate to the source directory
#    i.e. the uzipped folder.
#    cd C:\Users\user\Anaconda3\pkgs\scikit-bio-master


# 4. Using Notepad++ edit the setup.py file
#    "C:\Program Files (x86)\Notepad++\notepad++" setup.py
#   Find the line in which the variable ssw_extra_compile_args 
#   is defined and change it. You can comment on the previous version and 
#   redefine the variable in the new line:

#   #ssw_extra_compile_args = ['-Wno-error=declaration-after-statement']
#   ssw_extra_compile_args = []

# 5. Save the changes, close the editor and start the installation with
#    this command:
#    python setup.py install
#   I hope you do not receive error messages.

# 6. Exit from the directory that you insatlled skbio into
#    cd C:\Users\user\Documents

# 7. Open an Anaconda Python session (using the command python) and check if 
#    Scikit-Bio was installed correctly using print(skbio.art).
#    test insatllation
#   import skbio
#   print(skbio.art) """

import numpy as np
import pandas as pd
import os
import ecopy as ep
import random

def pathogen_analysis(file):
    """
    Function to perform pathogen analysis from files generated
    after running the analysis in bash
    """
    # Check if the file is empty
    if os.stat(file).st_size > 0:
        
        # Read in the file
        sample = pd.read_csv(file,
                             header=None)
        
        # Rename the column
        sample.columns= ['strains']
        # Get the frequency for each strain
        sample['strains'].value_counts()
        
        ####### Get the stain names on the strain column #####
        # Get sthe index of the starins
        index = [x for x in range(len(sample.loc[:,'strains']))]
        
        # Intialize  an empty list that will contain the strain names
        strain_names = []
        # Loop through every row while splitting then appending the strain name 
        # to sb
        for i in index:
            strain_names.append((sample.loc[:,'strains'][i]).split())
            
        # Create a dataframe from sb
        split_names = pd.DataFrame(strain_names)
        
        # Replace none with nan
        split_names = split_names.replace(to_replace='None',
                                          value=np.nan)
        
        # Replace nan with an empty string
        split_names = split_names.fillna(value='')
        # Get the genera dataframe
        genera_df = split_names.iloc[:,0]
        # Get the species dataframe
        species_df = split_names.iloc[:,0] + ' ' \
        + split_names.iloc[:,1]   
        # Get the strains dataframe
        strains_df = split_names.iloc[:,0] + ' ' + \
        split_names.iloc[:,1] + ' ' + split_names.iloc[:,2] + \
         ' ' + split_names.iloc[:,3] #+ ' ' + split_names.iloc[:,4]
        
        taxon_df = [genera_df,species_df,strains_df]
        taxon_names = ['Genus', 'Species', 'Strains']
        taxon = ['genera','species', 'strains']
        for j in range(len(taxon_df)):
    
            #### Make a dataframe of the taxon and the frequency ###
            taxon_freq= pd.DataFrame(taxon_df[j].value_counts())
            # Create an empty dataframe that will contain the taxon
            # and freqency information
            taxon_names[j] = pd.DataFrame()
            # Get the names of the taxon and set it to the taxon
            # column here for the genera
            taxon_names[j][taxon[j]] = taxon_freq.index.values
            # Add a column that will contain the frequency of each genera
            taxon_names[j]['freq'] = taxon_freq.iloc[:,0].values
        # Create a dictionary that will contain
        # the analysis for each taxon level
        final = {taxon[0]: taxon_names[0],
                 taxon[1]: taxon_names[1],
                 taxon[2]: taxon_names[2]}
    else:
        # If the file is empty run this block
        taxon = ['genera','species', 'strains']
        taxon_names = ['Genus', 'Species', 'Strains']
        taxon_names[0] = pd.DataFrame({'genera': 'non_detected','freq': 0},
                                       index=[0],columns= ['genera','freq'])
        taxon_names[1] = pd.DataFrame({'species': 'non_detected','freq': 0},
                                       index=[0],columns= ['species','freq'])
        taxon_names[2] = pd.DataFrame({'strains': 'non_detected','freq': 0},
                                       index=[0],columns= ['strains','freq'])
        
        final = {taxon[0]: taxon_names[0],
                 taxon[1]: taxon_names[1],
                 taxon[2]: taxon_names[2]}
        
    return final


def combined_results(files):
    """
    Function to combined pathogen analyses results
    for all samples analysed
    """
    # Intialize an empty list that will contain 
    # the sample names
    samples = []
    # Intialize a dictionary that will contain
    result_together = {}
    # Append the sample names to the samples list
    for sample in range(len(files)):
        #print(sample)
        #get the sample names
        samples.append(files[sample].split(sep='_')[1])
    
        result_together[samples[sample]] = (pathogen_analysis(file=files[sample]))

    return(result_together,samples)
    
    
def make_taxon_table(result_together, samples):
    """
    Function to make a pathogen taxon table
    """
    ##get a named list
    pathogens = pd.Series()
    for sample in samples:
        pathogens = pathogens.append(result_together[sample]['species']['species'])

    # Get the unique species    
    pathogens = pathogens.unique()
    d = {'pathogens': pathogens}
    taxon_table = pd.DataFrame(d)

    # Remove the non detected pathogens
    taxon_table = taxon_table[taxon_table['pathogens'] != 'non_detected']
    # Create a dataframe with the sample names with values set at zero
    zeros_dataframe = pd.DataFrame(data=0, index=np.arange(len(taxon_table.index)),\
                               columns= samples)
    # Set the index of the zeros dataframe
    zeros_dataframe.index = taxon_table.index
    # Create a frame list
    frame = [taxon_table,zeros_dataframe]
    # Concatenate the dataframes along the columns
    taxon_table = pd.concat(frame, axis=1)
    # Set the index  of the dataframe to the names of the pathogens
    taxon_table = taxon_table.set_index('pathogens')

    # Loop through every sample while getting the frequency for each pathogen
    for sample in samples:
        #print(sample)
        # Get the detect/pathogens for each sample
        detect = result_together[sample]['species']['species']
        # Get the index in a list form
        index = detect.index.tolist()
        # Get all the frequencies for the dtected pathogens
        frequency = result_together[sample]['species']['freq']
        # Loop
        for pathogen in taxon_table.index.tolist():
            for i in index:
                if (pathogen == detect[i]):
                    taxon_table.loc[pathogen,sample] = frequency[i]
                           
    return(taxon_table)
    


# Remember to modify function for any grouping vaiable besides treatment    
def presence_abscence(taxon_table_t, independent_variables):
    """
    Function to calculate the proportion of pathogens detected
    in a pathogen table matrix
    """
    # Find if the pathogens were detected in the samples 
    detect_table = taxon_table_t > 0

    detect_frame = [independent_variables,detect_table]
    detect_table2 = pd.concat(detect_frame,axis=1)

    # Get the number of times each pathogen as detected for each treatment
    detect = detect_table2.groupby(by=['Treatment']).sum()
    number_of_samples_per_treat = detect_table2.groupby(by=['Treatment']).count().iloc[:,0]
    frame2 = [number_of_samples_per_treat,detect]
    detect_df= pd.concat(frame2,axis=1)
    columns = detect_df.columns.tolist()
    columns[0] = 'number of samples'
    detect_df.columns = columns
    
    return detect_df
    

def join_tables(table1,table2,how='inner'):
    """
    Function for joining two pandas dataframes based on their indices 
    """
    table1_col_names = table1.columns
    # Copy the dataframes to avoid inplace assignment
    #to the original dataframe
    table1_copy = table1.copy()
    table2_copy = table2.copy()
    table1_copy['samples'] = table1_copy.index
    table2_copy['samples'] = table2_copy.index
    
    joined_table = table1_copy.join(other=table2_copy,
                              on='samples',
                              how=how,
                              lsuffix='_i',
                              rsuffix='_t')
    joined_table = joined_table.drop(columns=['samples','samples_i','samples_t'])
    table1_copy = joined_table[table1_col_names]
    table2_copy = joined_table.drop(columns=table1_col_names)
    
    return(table1_copy,table2_copy,joined_table)
    

    
def remove_zero_samples(taxon_table_t):
    """
    Function to remove samples that nothing was detected in them 
    to avoid errors with division by zero
    """
    sample_sum = taxon_table_t.sum(axis=1)
    # Get the boolean indices for samples that pathogens were not detected in
    logic_index = sample_sum == 0
    # Get the samples
    zeros_samples =list(sample_sum[logic_index].index)
    # Drop the samples with zero or that pathogens were not detcted in them
    taxon_table_t = taxon_table_t.drop(index=zeros_samples)
    
    return taxon_table_t


def estimate_diversity(taxon_table_t):
    """ Function to estimate common diversity matrics from a taxon table
    with samples as rows and observation/otus/pathogens as column
    the output of this function is a dataframe of the diversity matrices per sample
    """
    # Import the modules for diversity calculations and manipulating dataframes
    # Get the sample names which are the indices of the transposed taxon table
    table_samples = taxon_table_t.index.tolist()

    # Get the number of observed otus / pathogens per sample
    observed_pathogens = [] #initialize an empty list
    # Calculate the observed species per sample
    for sample in table_samples:
        observed_pathogens.append(observed_otus(table= taxon_table_t,
                                                sample_id= sample)) 

    # Estimate diversity
    shannon = ep.diversity(x = taxon_table_t, method = 'shannon')
    simpson = ep.diversity(x = taxon_table_t, method = 'simpson')
    species_richness = ep.diversity(x = taxon_table_t, method = 'spRich')
    eveness = ep.diversity(x = taxon_table_t, method = 'even')

    # Convert the estimates to Pandas series
    shannon = pd.Series(data = shannon, index = table_samples)
    simpson = pd.Series(data = simpson, index = table_samples)
    observed_pathogens = pd.Series(data=observed_pathogens, index= table_samples)
    species_richness = pd.Series(data = species_richness, index = table_samples)
    eveness = pd.Series(data = eveness, index = table_samples)

    diversity_dict = {'observed_species': observed_pathogens, 'species_richness': species_richness,
                  'eveness': eveness, 'shannon':shannon, 'simpson': simpson}
    diversity_table = pd.DataFrame(data = diversity_dict)
    
    return(diversity_table)
    

def stack_data(abund_table,ind_var_df, variable):
    """
    Function to stack an abundance dataframe along with one independent variable
    abund_table = a table normalized to relative abundances with samples on the rows and features on the columns
    ind_var_df = a dataframe with the independent variables in the same order as the samples
    """
    import pandas as pd
    import numpy as np
    
    # Ensure that ind_var_df is a dataframe
    ind_var_df = pd.DataFrame(ind_var_df)
    # Get the independent variable of choice
    ind_var = ind_var_df[variable]
    
    indvar_label = []
    pathogens_label = []
    relative_abundance = []
    
    # Get a list of all the column names
    colNames = list(abund_table)
    # Loop through every column while creating a stack of the independent
    #variable the column name (pathogens)
    # and the values within these columns
    count = 0
    for column in range(0,len(colNames)):
        indvar_label.append(ind_var)
        pathogen = pd.Series((np.repeat(colNames[column],len(ind_var))),name='pathogens')
        pathogens_label.append(pathogen)
        relative_abundance.append(abund_table[colNames[column]])
        count += 1
   
    Ind_var = pd.concat(indvar_label)
    pathogens = pd.concat(pathogens_label)
    pathogens.index = Ind_var.index
    Relative_abund = pd.concat(relative_abundance)
    frame = [Ind_var,pathogens,Relative_abund]
    stacked_table = pd.concat(frame,axis=1)
    stacked_table.columns = [variable,'pathogens','relative_abundance']
    
    return stacked_table


def remove_low_abundance_pathogens(abun_table,threshold = 0.05,group_low_abund=True):
    """
    Function to remove or group the less abundant pathogens
    """
    # Get the maximum value for each column/pathogens
    max_abund = abun_table.max()
    # Get the boolean vector for the less abundant pathogens
    less_abund_column_logic = max_abund < threshold
    # Get the less abundant columns/pathogens
    less_abund_columns = max_abund[less_abund_column_logic]
    # Sum the low abundance pathogens
    rare = abun_table[list(less_abund_columns.index)].sum(axis=1)
    # Drop the less abundant columns/pathogens
    abund_pathogens_table = abun_table.drop(columns=list(less_abund_columns.index))
    # Group the pathogens with low abundance
    if group_low_abund:
        # Create a frame for concatenation
        frame = [abund_pathogens_table,rare]
        # Concatenate the frame
        abund_pathogens_table = pd.concat(frame,axis=1)
        # Rename the last column, which is the column with the rare pathogens
        abund_pathogens_table.columns.values[len(abund_pathogens_table.columns)-1] = 'Rare'
        
    return abund_pathogens_table


def get_treatment_abund_table(independent_variable_df,taxon_table,group="Treatment"):
    # Get the independent variable of interest - here i get the treatment variable
    treatment_df = independent_variable_df[group].to_frame()
    _,_,treatment_abun_table = join_tables(table1=treatment_df,
                          table2=taxon_table)
    #get the sum of counts for each treatment
    treatment_abun_table = treatment_abun_table.groupby([group]).sum()
    #normalize to relative abundance per treatment in percentage
    treatment_abun_table = treatment_abun_table.apply(func=lambda x: x/sum(x),axis=1) * 100
    return treatment_abun_table

##################### these set of function are courtesy of Greg Carapaso  in
#################### the applied bioinformatics tutorial 

def observed_otus(table, sample_id):
    """
    Function to get the number of observed otus in a sample
    when samples are rows and features are columns
    """
    return sum([e > 0 for e in table.loc[sample_id,:]])
 
def get_observed_nodes(tree, table, sample_id, verbose=False):
    """
    Function to get observed nodes
    """
    observed_otus = [obs_id for obs_id in table.index if table[sample_id][obs_id] > 0]
    observed_nodes = set()
     # Iterate over the observed OTUs
    for otu in observed_otus:
        t = tree.find(otu)
        observed_nodes.add(t)
        if verbose:
            print(t.name, t.length, end=' ')
        for internal_node in t.ancestors():
            if internal_node.length is None:
                # We've hit the root
                if verbose:
                    print('')
                else:
                    if verbose and internal_node not in observed_nodes:
                        print(internal_node.length, end=' ')
                    observed_nodes.add(internal_node)
                    
    return observed_nodes


def phylogenetic_diversity(tree, table, sample_id, verbose=False):
    """
    Function to calculate phylogenetic_diversity
    """
    observed_nodes = get_observed_nodes(tree, table, sample_id, verbose=verbose)
    result = sum(o.length for o in observed_nodes)
    
    return result
    
    
def bray_curtis_distance(table, sample1_id, sample2_id):
    """
    Function to calculate bray-curtis distance
    """
    numerator = 0
    denominator = 0
    sample1_counts = table[sample1_id]
    sample2_counts = table[sample2_id]
    for sample1_count, sample2_count in zip(sample1_counts, sample2_counts):
        numerator += abs(sample1_count - sample2_count)
        denominator += sample1_count + sample2_count
        
    return numerator / denominator
    



def table_to_distances(table, pairwise_distance_fn):
    """
    Function to make a distance matrix
    """
    from skbio.stats.distance import DistanceMatrix
    from numpy import zeros
    sample_ids = table.columns
    num_samples = len(sample_ids)
    data = zeros((num_samples, num_samples))
    for i, sample1_id in enumerate(sample_ids):
        for j, sample2_id in enumerate(sample_ids[:i]):
            data[i,j] = data[j,i] = pairwise_distance_fn(table, sample1_id, sample2_id)
            
    return DistanceMatrix(data, sample_ids)
    
    

def unweighted_unifrac(tree, table, sample_id1, sample_id2, verbose=False):
    """
    Function to calculate unweighed unifrac distance
    """
    observed_nodes1 = get_observed_nodes(tree, table, sample_id1, verbose=verbose)
    observed_nodes2 = get_observed_nodes(tree, table, sample_id2, verbose=verbose)
    observed_branch_length = sum(o.length for o in observed_nodes1 | observed_nodes2)
    shared_branch_length = sum(o.length for o in observed_nodes1 & observed_nodes2)
    unique_branch_length = observed_branch_length - shared_branch_length
    unweighted_unifrac = unique_branch_length / observed_branch_length
    
    return unweighted_unifrac

##############################################################################



def generate_colors(n):
    """
    Function to generate n visually distinct RGB colours 
    """
    rgb_values = []
    hex_values = []
    r = int(random.random() * 256)
    g = int(random.random() * 256)
    b = int(random.random() * 256)
    step = 256 / n
    for _ in range(n):
        r += step
        g += step
        b += step
        r = int(r) % 256
        g = int(g) % 256
        b = int(b) % 256
        r_hex = hex(r)[2:]
        g_hex = hex(g)[2:]
        b_hex = hex(b)[2:]
        hex_values.append('#' + r_hex + g_hex + b_hex)
        rgb_values.append((r,g,b))
    return rgb_values, hex_values


