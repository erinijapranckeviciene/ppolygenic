#!/usr/bin/python
"""
Call: python this_code.py 

Code to test total_genotype_score functions 

"""
import sys
import itertools
import csv

from tgs_functions import total_genotype_score, markers, markers_plus, genotype_combinations_tgs, generate_marker_table, marker_table_to_csv,class_table_csv_for_testing,labels

##############################################################################
# Testing markers computing all genotype combinations and TGS  frequency table
##############################################################################

####### Injury markers

injury_gt_list='rs12722 TT, rs1800012 AA, rs679620 GG'

tgs_injury = total_genotype_score(injury_gt_list,markers['injury'])
print "Injury markers: ",injury_gt_list 
print "TGS computed: ", tgs_injury

co = genotype_combinations_tgs(markers['injury'])
tg=co[0]
ft=co[1]

print "Injury markers: all genotype combinations and their TGS:" 
for x in tg:
    print x

print "Injury markers: table of TGS frequency distribution:" 
for x in ft:
    print x

####### Endurance and power markers

endurance_power_gt_list='rs699 AA, rs1042713 AG, rs11549465 CT'

tgs_endurance_power = total_genotype_score(endurance_power_gt_list,markers['endurance_power'])
print "Endurance_power markers: ",endurance_power_gt_list 
print "TGS computed: ", tgs_endurance_power

co = genotype_combinations_tgs(markers['endurance_power'])
tg=co[0]
ft=co[1]

print "Endurance and power markers: all genotype combinations and their TGS:" 
for x in tg:
    print x

print "Endurance and power markers: table of TGS frequency distribution:" 
for x in ft:
    print x


####### Vitamin B12 markers

vitaminB12_gt_list='rs1801133 GG, rs602662 AA, rs526934 AG, rs1801222 CC'

tgs_vitaminB12 = total_genotype_score(vitaminB12_gt_list,markers['vitaminB12'])
print "Endurance_power markers: ",vitaminB12_gt_list 
print "TGS computed: ", tgs_vitaminB12

co = genotype_combinations_tgs(markers['vitaminB12'])
tg=co[0]
ft=co[1]

print "Vitamin B12 status markers: all genotype combinations and their TGS:" 
for x in tg:
    print x

print "Vitamin B12 status markers: table of TGS frequency distribution:" 
for x in ft:
    print x

##############################################################################
# Testing generation of synthetic data, csv format , output file names
#      synthetic_injury , synthetic_endurance_power, synthetic_vitaminB12 
##############################################################################

# Phenotype labels to assign at random import labels from tgs_functions.py
#labels={'injury':{'risk':0.50,'no-risk':0.5},'endurance_power':{'endurance':0.35,'mix':0.35,'power':0.30}}

####### Injury markers

# file name to write: chage accordingly
output_name_injury =  '/var/www/html/browse/sessions/injury_synthetic_data.txt'
print 'Output file name :' ,output_name_injury
# phenotype in markers
phenotype_injury = 'injury'
print 'Phenotype : ', phenotype_injury
#sample size = 100, change accordingly 
n = int(100) 
print 'Sample size : ', n

# generate sample for phenotypic markers
rez=generate_marker_table(markers_plus[phenotype_injury],n,labels[phenotype_injury])

#write csv for injury
marker_table_to_csv(rez,output_name_injury)
print 'Wrote synthetic data for injury' 

####### Endurance vs. Power  markers

# file name to write: chage accordingly
output_name_endurance_power =  '/var/www/html/browse/sessions/endurance_power_synthetic_data.txt'
print 'Output file name :' ,output_name_endurance_power
# phenotype in markers
phenotype_endurance_power = 'endurance_power'
print 'Phenotype : ', phenotype_endurance_power
#sample size = 100, change accordingly 
n = int(100) 
print 'Sample size : ', n

# generate sample for phenotypic markers
rez=generate_marker_table(markers_plus[phenotype_endurance_power],n,labels[phenotype_endurance_power])

#write csv for endurance power
marker_table_to_csv(rez,output_name_endurance_power)
print 'Wrote synthetic data for endurance power'

####### Endurance vs. Power  markers

# file name to write: chage accordingly
output_name_vitaminB12 =  '/var/www/html/browse/sessions/vitaminB12_synthetic_data.txt'
print 'Output file name :' ,output_name_vitaminB12
# phenotype in markers
phenotype_vitaminB12 = 'vitaminB12'
print 'Phenotype : ', phenotype_vitaminB12
#sample size = 100, change accordingly 
n = int(100) 
print 'Sample size : ', n

# generate sample for phenotypic markers
rez=generate_marker_table(markers[phenotype_vitaminB12],n,labels[phenotype_vitaminB12])

#write csv for endurance power
marker_table_to_csv(rez,output_name_vitaminB12)
print 'Wrote synthetic data for vitaminB12'
 
##############################################################################
# Testing reading of synthetic data
#      synthetic_injury.csv , synthetic_endurance_power.csv 
##############################################################################
print 'Test reading synthetic csv data for injury'
print 'Process the labeled data' 
##
#Test TGS for data reading from file
#

####### Injury markers
inputf = '/var/www/html/browse/sessions/injury_synthetic_data.txt'
outputf = '/var/www/html/browse/sessions/injury_testing.txt'
class_table = class_table_csv_for_testing(inputf,markers_plus['injury'],outputf)
print "injury data"
print class_table

####### Endurance power  markers
inputf = '/var/www/html/browse/sessions/endurance_power_synthetic_data.txt'
outputf = '/var/www/html/browse/sessions/endurance_testing.txt'
class_table = class_table_csv_for_testing(inputf,markers['endurance_power'],outputf)
print "endurance power data"
print class_table


####### Vilnius markers
inputf = '/var/www/html/browse/sessions/athletes.csv'
outputf = '/var/www/html/browse/sessions/vilnius_testing.txt'
class_table = class_table_csv_for_testing(inputf,markers_plus['endurance_power_vilnius'],outputf)
print "Vilnius data"
print class_table


####### Vitamin B12 markers
inputf = '/var/www/html/browse/sessions/vitaminB12_synthetic_data.txt'
outputf = '/var/www/html/browse/sessions/vitaminB12_testing.txt'
class_table = class_table_csv_for_testing(inputf,markers['vitaminB12'],outputf)
print class_table
print "Vitamin B12 data"
