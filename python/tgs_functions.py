# -*- coding: utf-8 -*-
"""
Created on Tue Nov 03 19:18:29 2015
Finished on Tue Dec 15 11:48:00 2015
@author: Erin

Meaning of scores:
    1 endurance , -1 power, 0 intermediate
    1 injury risk, -1 no risk, 0 intermediate 

"""

from math import erf,exp
from itertools import product

import numpy as np
import csv


################################################################################
# Data structure  markers to store information about the multimaarker phenotype
#
markers = \
    {'injury': \
     {'rs12722':  {'gene':'COL5A1','genotype_scores':{'TT':1,'CT':0,'CC':-1},'frequencies':{'TT':0.07,'CT':0.36,'CC':0.57}, 'conf': 1},\
      'rs1800012':{'gene':'COL1A1','genotype_scores':{'CC':0,'AC':0,'AA':-1},'frequencies':{'CC':0.66,'AC':0.30,'AA':0.04}, 'conf': 1},\
      'rs679620': {'gene':'MMP3',  'genotype_scores':{'AA':1,'AG':0,'GG':-1},'frequencies':{'AA':0.31,'AG':0.53,'GG':0.16}, 'conf': 0.75}},\
    'endurance_power': \
     {'rs699':     {'gene':'AGT', 'genotype_scores':{'AA': 1, 'AG': 0, 'GG': -1},'frequencies':{'AA': 0.0676, 'AG': 0.3848, 'GG': 0.5476},'allele_freq':{'A':0.26,'G':0.74}, 'conf': 1},\
      'rs1042713': {'gene':'ADRB2','genotype_scores':{'AA': 1, 'AG': 0, 'GG': -1},'frequencies':{'AA': 0.2304, 'AG': 0.4992, 'GG': 0.2704},'allele_freq':{'A':0.48,'G':0.52}, 'conf': 1},\
      'rs11549465':{'gene':'HIF1A','genotype_scores':{'CC': 1, 'CT': 0, 'TT': -1},'frequencies':{'CC': 0.8575, 'CT': 0.1370, 'TT': 0.0055},'allele_freq':{'C':0.926,'T':0.074}, 'conf': 1}},\
    'vitaminB12': \
     {'rs1801133':    {'gene':'MTHFR', 'genotype_scores':{'GG': 0, 'AG': -0.7, 'AA': -1},'frequencies':{'GG': 0.57, 'AG': 0.33, 'AA': 0.1},'conf': 1},\
      'rs602662': {'gene':'FUT2','genotype_scores':{'AA': 0.7, 'AG': 0, 'GG': -1},'frequencies':{'AA': 0.20, 'AG': 0.35, 'GG': 0.45}, 'conf': 1},\
      'rs526934':{'gene':'TCN1','genotype_scores':{'AA': 0, 'AG': -1, 'GG': -1},'frequencies':{'AA': 0.65, 'AG': 0.30, 'GG': 0.05}, 'conf': 1},\
      'rs1801222':{'gene':'CUBN','genotype_scores':{'CC': 1, 'CT': 0, 'TT': -1},'frequencies':{'CC': 0.65, 'CT': 0.30, 'TT': 0.05}, 'conf': 1}}}



##############################################################################
# Data structure  -markers- extended with additional block of different markers
markers_plus ={ \
# Injury
'injury': \
     {'rs12722':  {'gene':'COL5A1','genotype_scores':{'TT':1,'CT':0,'CC':-1},'frequencies':{'TT':0.07,'CT':0.36,'CC':0.57}, 'conf': 1},\
      'rs1800012':{'gene':'COL1A1','genotype_scores':{'CC':0,'AC':0,'AA':-1},'frequencies':{'CC':0.66,'AC':0.30,'AA':0.04}, 'conf': 1},\
      'rs679620': {'gene':'MMP3',  'genotype_scores':{'AA':1,'AG':0,'GG':-1},'frequencies':{'AA':0.31,'AG':0.53,'GG':0.16}, 'conf': 0.75}},\

# Vilnius endurance_power
     'endurance_power_vilnius': \
                {'rs1799752':     {'gene':'ACE', 'genotype_scores':{'DD': 1, 'ID': 0, 'II': -1},'allele_freq':{'I':0.43,'D':0.57}, 'conf': 1},\
                 'rs1815739':     {'gene':'ACTN3', 'genotype_scores':{'RR': 1, 'RX': 0, 'XX': -1},'allele_freq':{'R':0.65,'X':0.35}, 'conf': 1},\
                 'rs7293': {'gene':'MB','genotype_scores':{'AA': 1, 'AG': 0, 'GG': -1},'allele_freq':{'A':0.5,'G':0.5}, 'conf': 0.75},\
                 'rs17602729':{'gene':'AMPD1','genotype_scores':{'CC': 1, 'CT': -1, 'TT': -1},'allele_freq':{'C':0.85,'T':0.15}, 'conf': 0.75}},\

# Endurance_power
    'endurance_power': \
     {'rs699':     {'gene':'AGT', 'genotype_scores':{'AA': 1, 'AG': 0, 'GG': -1},'frequencies':{'AA': 0.0676, 'AG': 0.3848, 'GG': 0.5476},'allele_freq':{'A':0.26,'G':0.74}, 'conf': 1},\
      'rs1042713': {'gene':'ADRB2','genotype_scores':{'AA': 1, 'AG': 0, 'GG': -1},'frequencies':{'AA': 0.2304, 'AG': 0.4992, 'GG': 0.2704},'allele_freq':{'A':0.48,'G':0.52}, 'conf': 1},\
      'rs11549465':{'gene':'HIF1A','genotype_scores':{'CC': 1, 'CT': 0, 'TT': -1},'frequencies':{'CC': 0.8575, 'CT': 0.1370, 'TT': 0.0055},'allele_freq':{'C':0.926,'T':0.074}, 'conf': 1}}}

#################################################################
# Labels to generate synthetic data
labels={'injury':{'risk':0.50,'no-risk':0.5},'endurance_power':{'endurance':0.35,'mix':0.35,'power':0.30},'vitaminB12':{'high':0.20,'higher':0.20,'normal':0.20,'lower':0.20,'low':0.20}}
##################################################################
# TGS computation
def total_genotype_score(markers_list,markers_dict):
    """
    total_genotype_score(markers_list,markers_dict)
    Function: computes a total genotype score TGS

    Input: 
        markers_list: 
          string of comma separated pairs of a marker rs id
          followed by whitespace followed by a genotype. Example:
          markers_list='rs12722 TT, rs1800012 AA, rs679620 AG'
        markers_dict: 
          dictionary of markers. Example: markers['injury']  

    Return value: 
        float TGS:
          total genotype score computed as a dot product of two vectors 
          a vector of individual genotype scores and 
          a vector of marker confidence scores 
    """

    id_gt_pairs = markers_list.split(',')
    n_gt_pairs  = len(id_gt_pairs)
#    if n_gt_pairs <> 3 : # number of markers may be arbitrary, here 3  
    if n_gt_pairs < 2 : # test if it is multimarker, STOP if it is one marker
            return -101

    gt_score    = list() # vector of marker scores
    marker_conf = list() # vector of confidence scores
    for id_gt in id_gt_pairs:
        if len(id_gt.split()) <> 2:
            print "<p> part is missing</p>"
            return -104

        (identity,genotype) = id_gt.split()
        #print identity, genotype
        if identity in markers_dict:
            marker_conf.append(float(markers_dict[identity]['conf']))
        else:
            print "<p> Invalid identity </p>"
            return -102 # Invalid identity

        if genotype in markers_dict[identity]['genotype_scores']:
            gt_score.append(float(markers_dict[identity]['genotype_scores'][genotype]))
        else:
            print "<p> Invalid genotype </p>"
            return -103 # Invalid genotype    
      
    # apply error function to map the score into [-1 , 1] interval
    score = sum([gt_score[i]*marker_conf[i] for i in range(0,n_gt_pairs)])

    # normalization constant is a sum of confidences
    normalization_constant = sum(markers_dict[i]['conf'] for i in markers_dict.keys())

    # normalized score
    normalized_score = score/normalization_constant

    #erf_score = erf(normalized_score)
    # apply sigmoid function and constant shift into interval [-1,1]
    #sigmoid_score = 1/(1+exp(-normalized_score))-0.5
    # return original score and a normalized_score
    #return [score,normalized_score,erf_score,sigmoid_score]
    return [score,normalized_score]

#####################################################################################  
# test total_genotype_score
def genotype_combinations_tgs(markers_dict):

    """
    Testing TGS function on all possible genotype combinations 
    It works with any type of markers and any number of marker states
    """
    genotypes_list    = list() # list of genotypes of each marker
    marker_names_list = list() # list of marker names (rs ids)

    for marker_id in markers_dict:
        genotypes_list.append(markers_dict[marker_id]['genotype_scores'].keys())
        marker_names_list.append(marker_id)
    
    # Enumerate explicitly all possible combinations of genotypes.  
    # Generate string for a command: list(itertools.product(['TT','CT','CC'],...,['AA','AG','GG']))

    arg_i = list() # arguments for itertools.product()
    for marker_genotypes in genotypes_list:
        arg_i.append(str(marker_genotypes))

    #Evaluate itertools.product command and get all combinations of genotype indexes
    genotype_combinations = eval("list(product("+",".join(arg_i)+"))")

    tgs_vec               = list()
    table_gt_combinations = list()
    header = "Index"+","+",".join(marker_names_list)+","+",".join(["TGS","TGS_normalized"])
    table_gt_combinations.append(header)

    for gt_index,gt_combination in enumerate(genotype_combinations):
        arg_i=list() # arg_i parameter string to call TGS function
        for index,genotype in enumerate(gt_combination):
            arg_i.append(marker_names_list[index]+" "+genotype)
        
        tgs = total_genotype_score(",".join(arg_i),markers_dict) # Compute TGS
        tgs_str = ["%.4f" % tgs[i] for i in range(len(tgs))]
        line=str(gt_index)+","+",".join(gt_combination)+","+",".join(tgs_str) #Add row to the table
        table_gt_combinations.append(line) #Add row to the table
        tgs_vec.append(tgs[0]) # TGS non-normalized 

    tgs_frequency_dist=list()     # Create dictionary for TGS frequency distribution 
    tgs_frequency_dist.append(['Value','Count'])
    tgs_values=list(set(tgs_vec)) # Count unique score occurrences
    for x in tgs_values:
        tgs_frequency_dist.append([x,tgs_vec.count(x)])
    # sort by value
    #print tgs_frequency_dist
    return [table_gt_combinations,tgs_frequency_dist]


#####################################################################################  
# generate sample of size n for multiple markers 
# add random phenotypic information
def generate_marker_table(markers,n,labels):
    """
    Input:
        markers - markers[phenotype] Ex: markers['injury']
        n -  number of rows in the sample
        labels - dictionary of phenotype labels Ex: {'endurance':0.35,'mix':0.35,'power':0.30} 
    Output:
        table in which each column represents a marker
        last column  class represents phenotypic information 
    """
    table = {}
    for rsid in markers.keys():
        m = np.random.choice(markers[rsid]['frequencies'].keys(),n,markers[rsid]['frequencies'].values())
        table[rsid] = m
    table['class'] = np.random.choice(labels.keys(),n,labels.values())
    return table


#####################################################################################  
# write synthetic multi marker table to file
def marker_table_to_csv(marker_table, output_csv_file_name):

    # create object to write to csv file
    csv_writer = csv.writer(open(output_csv_file_name,'w'),delimiter=',')
    header = marker_table.keys()
    n = len(marker_table[header[0]])
    csv_writer.writerow(header)
    for i in range(n):
        csv_writer.writerow([marker_table[k][i] for k in header])

####### Data with labels prepare for testing ##############################
# Read markers from csv , compute TGS, 
# write markers and TGS to csv for testing in another envirinment
def class_table_csv_for_testing(input_csv_file_name, phenotype_markers, output_csv_file_name):
    class_table = list()
    with open(input_csv_file_name) as csvfile:
        reader = csv.DictReader(csvfile)
        #only extract keys that begin with rs
        for row in reader:
            markers = list()
            rsids = [s for s in row.keys() if 'rs' in s] 
            #print rsids
            # create markers_list input to tgs_function
            for rsid in rsids:
                markers.append(' '.join([rsid,row[rsid]]))
            markers_list = ','.join(markers)
            # compute total genotype score
            rez =total_genotype_score(markers_list, phenotype_markers) 
            rezs = [str(r) for r in rez]
            s=','.join(rezs) 

            line = ','.join([markers_list,row['class'],','.join(rezs)])
            class_table.append(line.split(','))
    #print output_csv_file_name
    with open(output_csv_file_name,'w') as outf:
        for row in class_table:
            row_str=",".join(row)
            #print row_str
            outf.write("%s\n" % row_str)

    return class_table
####### END OF FUNCTIONS ##############################
  



