###############################################################################

from __future__ import division
import csv
import time
import argparse
from screenFun import *
import sys
import numpy as np


###############################################################################    
# Version number

current_version = '0.7'


###############################################################################    

# Parses input using argparse module

# Initiates input parser
parser = argparse.ArgumentParser(description='Combines data')

# Non-optional arguments:
parser.add_argument('res_file1', help='Result file one', type=str)

parser.add_argument('res_file2', help='Result file two', type=str)

parser.add_argument('name', help='Name for output files', type=str)

# Optional arguments:
parser.add_argument('-p', '--processors', dest='nums',
                help='Number of processors to use', type=int,
                default=20)

parser.add_argument('-r', '--reference',
                help='Location of reference files', type=str,
                default='GenRef')

parser.add_argument('-of', '--override_file', action='store_true',
                help='Override Result file output location')

# Saves all input to object args
args = parser.parse_args()


###############################################################################
# Processes and checks input arguments

# Checks if parallel requested
if args.nums == 1:
    single = True
else:
    single = False

# Determines output file
if args.override_file:
    file_out = args.name
else:
    file_out = os.path.join('Results', args.name)

# Checks if can write to output
try:
    with open(file_out + '_record.txt', 'w') as out_open:
        pass
except:
    sys.exit('Cannot write to output file:\n' + file_out + '\n'
                + 'Use -of or --override_file to change')


###############################################################################
# Finds record files

print('Retrieving records')

name1 = args.res_file1[: -4]
rec_file1 = name1 + '_record.txt'

name2 = args.res_file2[: -4]
rec_file2 = name2 + '_record.txt'

try:
    # Parses record file 1
    with open(rec_file1, 'r') as rec_open:
        rec_csv = csv.reader(rec_open, delimiter='\t')
        version1 = rec_csv.next()[1]
        last_time1 = rec_csv.next()[1]
        unt_file1 = rec_csv.next()[1]
        trt_file1 = rec_csv.next()[1]
        file_out1 = rec_csv.next()[1]
        thresh1 = int(rec_csv.next()[1])
        screen_type1 = rec_csv.next()[1]
        K1 = float(rec_csv.next()[1])
        neg_name1 = rec_csv.next()[1]
        split_mark1 = rec_csv.next()[1]
        zero_files1 = rec_csv.next()[1]
        off_rate1 = float(rec_csv.next()[1])
        like_fun1 = eval(rec_csv.next()[1])
        I_step1 =  float(rec_csv.next()[1])

    # Parses record file 2
    with open(rec_file2, 'r') as rec_open:
        rec_csv = csv.reader(rec_open, delimiter='\t')
        version2 = rec_csv.next()[1]
        last_time2 = rec_csv.next()[1]
        unt_file2 = rec_csv.next()[1]
        trt_file2 = rec_csv.next()[1]
        file_out2 = rec_csv.next()[1]
        thresh2 = int(rec_csv.next()[1])
        screen_type2 = rec_csv.next()[1]
        K2 = float(rec_csv.next()[1])
        neg_name2 = rec_csv.next()[1]
        split_mark2 = rec_csv.next()[1]
        zero_files2 = rec_csv.next()[1]
        off_rate2 = float(rec_csv.next()[1])
        like_fun2 = eval(rec_csv.next()[1])
        I_step2 =  float(rec_csv.next()[1])

except IOError:
    sys.exit('Record of result file not found\n'
                + 'Change file name or rerun analysis')

if version1 != version2:
    sys.exit('Versions not consistent\n'
                + 'Rerun analysis')

if version1 != current_version:
    sys.exit('Version number not current\n'
                + 'Rerun analysis')


###############################################################################
# Pulls in gene info, IDs, as well as GO terms

print('Retrieving gene information')

# Retrieves gene info
geneID2Name, geneID2Info, geneName2ID, geneEns2Name = retrieveInfo()

# Retrieves GO information
geneID2Comp, geneID2Proc, geneID2Fun = retrieveGO()


###############################################################################
# Pulls in untreated and treated counts, and filters by defined threshold

print('Filtering reads')

# Retreives filtered counts for auxilary function
untreated1, treated1, stats1, time_zero1 = filterCounts(unt_file1,
                                                trt_file1, thresh1,
                                                zero_files1)

untreated2, treated2, stats2, time_zero2 = filterCounts(unt_file2,
                                                trt_file2, thresh2,
                                                zero_files2)


###############################################################################
# Calculates enrichment values

print('Calculating enrichment values')

# Retrieves enrichment values from auxilary function
element_rhos1, gene_rhos1, neg_rhos1, tar_rhos1, gene_ref1 = enrich_all(untreated1,
		treated1, neg_name1, split_mark1, K1, time_zero1)

element_rhos2, gene_rhos2, neg_rhos2, tar_rhos2, gene_ref2 = enrich_all(untreated2,
		treated2, neg_name2, split_mark2, K2, time_zero2)


###############################################################################
# Determines grid search

print('Determining search grid')

I_step = min([I_step1, I_step2])

# Auxiliary function determines single grid and unifies gene IDs
add_gene_rhos1, add_gene_rhos2, gene_span = comboSpan(gene_rhos1, gene_rhos2)


###############################################################################
# Retrieves and combines likelihoods using auxiliary functions

print('Calculating first likelihoods')

data1 = retrieveLikelihoods(add_gene_rhos1, neg_rhos1, tar_rhos1, off_rate1,
                                like_fun1, args.nums, gene_span, I_step)

print('Calculating second likelihoods')

data2 = retrieveLikelihoods(add_gene_rhos2, neg_rhos2, tar_rhos2, off_rate2,
                                like_fun2, args.nums, gene_span, I_step)

print('Combining likelihoods')

geneI, geneL, geneInterval = retrieveCombo(data1, data2, gene_span, I_step)


###############################################################################
# Writes output in human readable format

print('Outputing file')

with open(file_out + '.csv','w') as out_open:
    out_csv = csv.writer(out_open, delimiter=',',lineterminator = '\n')

    # Writes a header
    out_csv.writerow(['#GeneID','Symbol','GeneInfo','Localization','Process',
                        'Function','# elements 1','Effect 1','L stat 1',
                        '# elements 2','Effect 2','L stat 2',
                        'Combo effect', 'Combo L', 'L p-value',
			'Min I', 'Max I', 'Elements'])

    for gene in gene_span:

        # Finds appropriate IDs
	if gene in geneID2Name:
	    geneID = gene
	    name = geneID2Name[geneID]

	elif gene in geneName2ID:
	    name = gene
	    geneID = geneName2ID[gene]

        elif gene in geneEns2Name:
            name = geneEns2Name[gene]
            geneID = geneName2ID[name]

	else:
	    geneID = gene
	    name = gene

        # Retrieves information about gene
        name = geneID2Name[geneID]
        info = geneID2Info[geneID]
        comp = geneID2Comp[geneID]
        proc = geneID2Proc[geneID]
        fun = geneID2Fun[geneID]

        # Retrieves analysis of gene; rounding where appropriate
        num1 = len(add_gene_rhos1[gene])
        num2 = len(add_gene_rhos2[gene])

        # Retreives likelihood states
        effect1 = sigDig(data1[0][gene])
        effect2 = sigDig(data2[0][gene])
        rat1 = sigDig(data1[1][gene])
        rat2 = sigDig(data2[1][gene])

        # Retrieves combo scores
        effect = sigDig(geneI[gene])
        rat = sigDig(geneL[gene])
        pRat = 'N/A'
        min_CI = geneInterval[gene][0]
        max_CI = geneInterval[gene][1]

        # Writes to file
        out_csv.writerow([gene, name, info, comp, proc, fun,
                                num1, effect1, rat1,
                                num2, effect2, rat2,
                                effect, rat, pRat, min_CI, max_CI])


###############################################################################
# Saves record files

# Saves record for downstream analysis
with open(file_out + '_record.txt', 'w') as rec_open:
    rec_csv = csv.writer(rec_open, delimiter='\t')
    rec_csv.writerow(['analyzeCombo.py version', current_version])
    rec_csv.writerow(['Data/Time', time.strftime("%d:%H:%M:%eS")])
    rec_csv.writerow(['Result file 1', args.res_file1])
    rec_csv.writerow(['Result file 2', args.res_file2])
    rec_csv.writerow(['Untreated file 1', unt_file1])
    rec_csv.writerow(['Untreated file 2', unt_file2])
    rec_csv.writerow(['Treated file 1', trt_file1])
    rec_csv.writerow(['Treated file 2', trt_file2])
    rec_csv.writerow(['Output file', file_out])
    rec_csv.writerow(['Threshold 1', thresh1])
    rec_csv.writerow(['Threshold 2', thresh2])
    rec_csv.writerow(['Screen Type 1', screen_type1])
    rec_csv.writerow(['Screen Type 2', screen_type2])
    rec_csv.writerow(['Normalization 1', K1])
    rec_csv.writerow(['Normalization 2', K2])
    rec_csv.writerow(['Negative name 1', neg_name1])
    rec_csv.writerow(['Negative name 2', neg_name2])
    rec_csv.writerow(['Split mark 1', split_mark1])
    rec_csv.writerow(['Split mark 2', split_mark2])
    rec_csv.writerow(['Time zero files 1', zero_files1])
    rec_csv.writerow(['Time zero files 2', zero_files2])
    rec_csv.writerow(['Off-target rate 1', off_rate1])
    rec_csv.writerow(['Off-target rate 2', off_rate2])
    rec_csv.writerow(['Likelihood function 1', like_fun1.__name__])
    rec_csv.writerow(['Likelihood function 2', like_fun2.__name__])
    rec_csv.writerow(['Draw Number 1', int(np.median(map(len, add_gene_rhos1.values())))])
    rec_csv.writerow(['Draw Number 2', int(np.median(map(len, add_gene_rhos2.values())))])
    rec_csv.writerow(['I step', I_step])
    rec_csv.writerow(['Number of processers', args.nums])
    
# Saves permanent record
with open(os.path.join('Records', 'analyzeCombo' + time.strftime("%d%H%M%S")), 'w') as back_open:
    rec_csv = csv.writer(back_open, delimiter='\t')
    rec_csv.writerow(['analyzeCombo.py version', current_version])
    rec_csv.writerow(['Data/Time', time.strftime("%d:%H:%M:%eS")])
    rec_csv.writerow(['Result file 1', args.res_file1])
    rec_csv.writerow(['Result file 2', args.res_file2])
    rec_csv.writerow(['Untreated file 1', unt_file1])
    rec_csv.writerow(['Untreated file 2', unt_file2])
    rec_csv.writerow(['Treated file 1', trt_file1])
    rec_csv.writerow(['Treated file 2', trt_file2])
    rec_csv.writerow(['Output file', file_out])
    rec_csv.writerow(['Threshold 1', thresh1])
    rec_csv.writerow(['Threshold 2', thresh2])
    rec_csv.writerow(['Screen Type 1', screen_type1])
    rec_csv.writerow(['Screen Type 2', screen_type2])
    rec_csv.writerow(['Normalization 1', K1])
    rec_csv.writerow(['Normalization 2', K2])
    rec_csv.writerow(['Negative name 1', neg_name1])
    rec_csv.writerow(['Negative name 2', neg_name2])
    rec_csv.writerow(['Split mark 1', split_mark1])
    rec_csv.writerow(['Split mark 2', split_mark2])
    rec_csv.writerow(['Time zero files 1', zero_files1])
    rec_csv.writerow(['Time zero files 2', zero_files2])
    rec_csv.writerow(['Off-target rate 1', off_rate1])
    rec_csv.writerow(['Off-target rate 2', off_rate2])
    rec_csv.writerow(['Likelihood function 1', like_fun1.__name__])
    rec_csv.writerow(['Likelihood function 2', like_fun2.__name__])
    rec_csv.writerow(['Draw Number 1', int(np.median(map(len, add_gene_rhos1.values())))])
    rec_csv.writerow(['Draw Number 2', int(np.median(map(len, add_gene_rhos2.values())))])
    rec_csv.writerow(['I step', I_step])
    rec_csv.writerow(['Number of processers', args.nums])


###############################################################################
