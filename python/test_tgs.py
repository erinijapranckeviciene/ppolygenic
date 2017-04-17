#!/usr/bin/python

import cgi, cgitb, csv
from tgs_functions import markers,markers_plus,labels,total_genotype_score,genotype_combinations_tgs,generate_marker_table,marker_table_to_csv,class_table_csv_for_testing

form = cgi.FieldStorage()
# which button?

# injury group
injury_button         = form.getvalue('injury_markers_submit')
injury_gt_comb_button = form.getvalue('injury_genotype_combinations')
injury_markers        = form.getvalue('injury_markers')
injury_synthetic_button = form.getvalue('injury_synthetic')
injury_analyze_button   = form.getvalue('injury_analyze')

# endurance_power group
endurance_power_button         = form.getvalue('endurance_power_submit')
endurance_power_gt_comb_button = form.getvalue('endurance_power_gt_combinations')
endurance_power_markers        = form.getvalue('endurance_power_markers')
endurance_power_synthetic_button = form.getvalue('endurance_power_synthetic')
endurance_power_analyze_button   = form.getvalue('endurance_power_analyze')

# vitaminB12 group 
vitaminb12_button         = form.getvalue('vitaminb12_submit')
vitaminb12_gt_comb_button = form.getvalue('vitaminb12_genotype_combinations')
vitaminb12_markers        = form.getvalue('vitaminb12_markers')
vitaminb12_synthetic_button = form.getvalue('vitaminb12_synthetic')
vitaminb12_analyze_button   = form.getvalue('vitaminb12_analyze')


# vilnius
vilnius_button = form.getvalue('vilnius_analyze')
vilnius_gt_comb_button = form.getvalue('vilnius_genotype_combinations')

print "Content-type:text/html\r\n\r\n"
print "<html>"
print "<head>"
print "<title>Testing functions</title>"
print "</head>"
print "<body>"
# injury group
#print "<p>injury_button:  %s <p>" % (injury_button)
#print "<p>injury_gt_comb_button:  %s <p>" % (injury_gt_comb_button)


# endurance_power group
#print "<p>endurance_power_button:  %s <p>" % (endurance_power_button)
#print "<p>endurance_power_gt_comb_button:  %s <p>" % (endurance_power_gt_comb_button)

if injury_button:
    print "<p>injury_markers:  %s <p>" % (injury_markers)
    # submit input to the function
    tgs1 = total_genotype_score(injury_markers,markers_plus['injury'])
    print "<p>TGS and TGS normalized: %s </p>" % (tgs1)

##########################################################
if injury_gt_comb_button:
    # submit input to the function
    tables = genotype_combinations_tgs(markers_plus['injury'])
    gt_comb = tables[0]
    freq_dist = tables[1]
    print "<p>"
    for row in gt_comb:
         print " %s <br>" % (row)
    print "<br>"
    for row in freq_dist:
         print " %s <br>" % (row)

    print "</p>"

##########################################################
if injury_synthetic_button:
    table = generate_marker_table(markers_plus['injury'],300,labels['injury'])
    print "<p> Generated table is in  synthetic_data.txt</p>"
    print "<p>"
    header = table.keys()
    for name in header:
        print"%s" %(name)
    print "<br>"

    for row in range(300):
        line = [table[column][row] for column in header]
        row_str = ",".join(line)
        print "%s " % (row_str)
        print "<br>"
    print "</p>"

    fout = "/var/www/html/browse/sessions/injury_synthetic_data.txt"
    marker_table_to_csv(table,fout)

#########################################################
if injury_analyze_button:
    print "<p> Generated injury data  to test the model<br>"
    inputf = "/var/www/html/browse/sessions/injury_synthetic_data.txt"
    outputf = "/var/www/html/browse/sessions/injury_testing.txt"
    class_table = class_table_csv_for_testing(inputf,markers_plus['injury'],outputf)
    for row in class_table:
       row_str = ",".join(row)
       print "%s " % (row_str)
       print "<br>"
       
    print "</p>"

#################################################################
if endurance_power_button:
    print "<p>endurance_power_markers:  %s <p>" % (endurance_power_markers)
    # submit input to the function
    tgs1 = total_genotype_score(endurance_power_markers,markers_plus['endurance_power'])
    print "<p>TGS and TGS normalized: %s </p>" % (tgs1)

##################################################################
if endurance_power_gt_comb_button:
    # submit input to the function
    tables = genotype_combinations_tgs(markers_plus['endurance_power'])
    gt_comb = tables[0]
    freq_dist = tables[1]
    print "<p>"
    for row in gt_comb:
         print " %s <br>" % (row)
    print "<br>"
    for row in freq_dist:
         print " %s <br>" % (row)

    print "</p>"

##################################################################
if endurance_power_synthetic_button:
    table = generate_marker_table(markers_plus['endurance_power'],300,labels['endurance_power'])
    print "<p> Generated table is in markers table data</p>"
    print "<p>"
    header = table.keys()
    for name in header:
        print"%s" %(name)
    print "<br>"

    for row in range(300):
        line = [table[column][row] for column in header]
        row_str = ",".join(line)
        print "%s " % (row_str)
        print "<br>"
    print "</p>"

    fout = "/var/www/html/browse/sessions/endurance_power_synthetic_data.txt"
    marker_table_to_csv(table,fout)

#########################################################
if endurance_power_analyze_button:
    print "<p> Endurance power generated data to test the model<br>"
    inputf = "/var/www/html/browse/sessions/endurance_power_synthetic_data.txt"
    outputf = "/var/www/html/browse/sessions/endurance_power_testing.txt"
    class_table = class_table_csv_for_testing(inputf,markers_plus['endurance_power'],outputf)
    for row in class_table:
       row_str = ",".join(row)
       print "%s " % (row_str)
       print "<br>"
       
    print "</p>"


##########################################################################
if vitaminb12_button:
    print "<p>vitaminb12_markers:  %s <p>" % (vitaminb12_markers)
    # submit input to the function
    tgs1 = total_genotype_score(vitaminb12_markers,markers['vitaminB12'])
    print "<p>TGS and TGS normalized: %s </p>" % (tgs1)

##########################################################
if vitaminb12_gt_comb_button:
    #print "<p> button resp</p>"
    # submit input to the function
    tables = genotype_combinations_tgs(markers['vitaminB12'])
    gt_comb = tables[0]
    freq_dist = tables[1]
    print "<p>"
    for row in gt_comb:
         print " %s <br>" % (row)
    print "<br>"
    for row in freq_dist:
         print " %s <br>" % (row)

    print "</p>"

##########################################################
if vitaminb12_synthetic_button:
    table = generate_marker_table(markers['vitaminB12'],300,labels['vitaminB12'])
    print "<p> Generated table is in  ...synthetic_data.txt</p>"
    print "<p>"
    header = table.keys()
    for name in header:
        print"%s" %(name)
    print "<br>"

    for row in range(300):
        line = [table[column][row] for column in header]
        row_str = ",".join(line)
        print "%s " % (row_str)
        print "<br>"
    print "</p>"

    fout = "/var/www/html/browse/sessions/vitaminB12_synthetic_data.txt"
    marker_table_to_csv(table,fout)

#########################################################
if vitaminb12_analyze_button: 
    print "<p> Generated vitamin B12 status data  to test the model<br>"
    inputf = "/var/www/html/browse/sessions/vitaminB12_synthetic_data.txt"
    outputf = "/var/www/html/browse/sessions/vitaminB12_testing.txt"
    class_table = class_table_csv_for_testing(inputf,markers['vitaminB12'],outputf)
    for row in class_table:
       row_str = ",".join(row)
       print "%s " % (row_str)
       print "<br>"
       
    print "</p>"


#########################################################
if vilnius_button:
    print "<p> LT endurance power athletes <br>"
    inputf = "/var/www/html/browse/sessions/athletes.csv"
    outputf = "/var/www/html/browse/sessions/vilnius_testing.txt"
    class_table = class_table_csv_for_testing(inputf,markers_plus['endurance_power_vilnius'],outputf)
    for row in class_table:
       row_str = ",".join(row)
       print "%s " % (row_str)
       print "<br>"
       
    print "</p>"

##################################################################
if vilnius_gt_comb_button:
    # submit input to the function
    tables = genotype_combinations_tgs(markers_plus['endurance_power_vilnius'])
    gt_comb = tables[0]
    freq_dist = tables[1]
    print "<p>"
    for row in gt_comb:
         print " %s <br>" % (row)
    print "<br>"
    for row in freq_dist:
         print " %s <br>" % (row)

    print "</p>"


print "</body>"
print "</html>"
