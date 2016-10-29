###############################################################################
# David Morgens
# 04/06/2016
###############################################################################
#

from __future__ import division
import csv
import time
import argparse
from screenFun import *
import sys
import numpy as np


###############################################################################    
# Version number

current_version = '1.0'


###############################################################################    

# Parses input using argparse module

# Initiates input parser
parser = argparse.ArgumentParser(description='Combines data for two screens')

# Non-optional arguments:
parser.add_argument('res_file1', help='Result file one', type=str)

parser.add_argument('res_file2', help='Result file two', type=str)

parser.add_argument('name', help='Name for output files', type=str)

# Optional arguments:
parser.add_argument('-p', '--processors', dest='nums',
                help='Number of processors to use; default is 20', type=int,
                default=20)

parser.add_argument('-r', '--reference',
                help='Location of reference files. Default is GenRef', type=str,
                default='GenRef')

parser.add_argument('-of', '--override_file', action='store_true',
                help='Override Result file output location.')

parser.add_argument('-m', '--mouse', action='store_true')

parser.add_argument('-s', '--span', type=int, default=2)

# Saves all input to object args
args = parser.parse_args()


###############################################################################
# Processes and checks input arguments

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

stats1, files1, info1, param1 = retrieveRecord(args.res_file1, current_version)
stats2, files2, info2, param2 = retrieveRecord(args.res_file2, current_version)

script1, version1, last_time1 = stats1
unt_file1, trt_file1, zero_files1, file_out1 = files1
screen_type1, neg_name1, split_mark1, exclude1 = info1
thresh1, K1, back1, I_step1, scale1, draw_num1 = param1

script2, version2, last_time2 = stats2
unt_file2, trt_file2, zero_files2, file_out2 = files2
screen_type2, neg_name2, split_mark2, exclude2 = info2
thresh2, K2, back2, I_step2, scale2, draw_num2 = param2

if version1 != version2:
    sys.exit('Error: Versions not consistent\n'
                + 'Rerun analysis')


###############################################################################
# Pulls in gene info, IDs, as well as GO terms

print('Retrieving gene information')

# Retrieves gene info
geneID2Name, geneID2Info, geneName2ID, geneEns2Name = retrieveInfo(mouse=args.mouse)

# Retrieves GO information
geneID2Comp, geneID2Proc, geneID2Fun = retrieveGO()


###############################################################################
# Pulls in untreated and treated counts, and filters by defined threshold

print('Filtering reads')

# Retreives filtered counts for auxilary function
untreated1, treated1, stats1, time_zero1 = filterCounts(unt_file1,
                                                trt_file1, thresh1,
                                                zero_files1, exclude1)

untreated2, treated2, stats2, time_zero2 = filterCounts(unt_file2,
                                                trt_file2, thresh2,
                                                zero_files2, exclude2)


###############################################################################
# Calculates enrichment values

print('Calculating enrichment values')

# Retrieves enrichment values from auxilary function
element_rhos1, gene_rhos1, neg_rhos1, tar_rhos1, gene_ref1 = enrich_all(untreated1,
		treated1, neg_name1, split_mark1, K1, time_zero1, back1)

element_rhos2, gene_rhos2, neg_rhos2, tar_rhos2, gene_ref2 = enrich_all(untreated2,
		treated2, neg_name2, split_mark2, K2, time_zero2, back2)


###############################################################################
#

# Selects background
if back1 == 'all':
    back_rhos1 = neg_rhos1 + tar_rhos1

elif back1 == 'neg':
    back_rhos1 = neg_rhos1

elif back1 == 'tar':
    back_rhos1 = tar_rhos1

else:
    sys.exit('Error: Unrecognized background option: ' + back1)

if back2 == 'all':
    back_rhos2 = neg_rhos2 + tar_rhos2

elif back2 == 'neg':
    back_rhos2 = neg_rhos2

elif back2 == 'tar':
    back_rhos2 = tar_rhos2

else:
    sys.exit('Error: Unrecognized background option: ' + back2)


###############################################################################
# Determines grid search

print('Determining search grid')

I_step = min([I_step1, I_step2])
scale = max([scale1, scale2])

# Auxiliary function determines single grid and unifies gene IDs
add_gene_rhos1, add_gene_rhos2, gene_span = comboSpan(gene_rhos1, gene_rhos2, args.span)


###############################################################################
# Retrieves and combines likelihoods using auxiliary functions

print('Calculating first likelihoods')

data1 = retrieveLikelihoods(add_gene_rhos1, back_rhos1, args.nums,
                                gene_span, scale, I_step)

print('Calculating second likelihoods')

data2 = retrieveLikelihoods(add_gene_rhos2, back_rhos2, args.nums,
                                gene_span, scale, I_step)

print('Combining likelihoods')

geneI, geneL, geneInterval = retrieveCombo(data1, data2, gene_span, I_step)


###############################################################################
# Writes output in human readable format

print('Outputing file')

with open(file_out + '.csv','w') as out_open:
    out_csv = csv.writer(out_open, delimiter=',',lineterminator = '\n')

    # Writes a header
    out_csv.writerow(['#GeneID','Symbol','GeneInfo','Localization','Process',
                        'Function','# elements 1','casTLE Effect 1','casTLE Score 1',
                        '# elements 2','casTLE Effect 2','casTLE Score 2',
                        'Combo casTLE Effect', 'Combo casTLE Score', 'casTLE p-value',
			'Minimum Effect Estimate', 'Maximum Effect Estimate'])

    for gene in gene_span:

        # Finds appropriate IDs
        geneID, name, entrez = retrieveIDs(gene, geneID2Name, geneName2ID, geneEns2Name)

        # Retrieves information about gene
        info = geneID2Info[entrez]
        comp = geneID2Comp[entrez]
        proc = geneID2Proc[entrez]
        fun = geneID2Fun[entrez]

        # Retrieves analysis of gene; rounding where appropriate
        num1 = len(add_gene_rhos1[gene])
        num2 = len(add_gene_rhos2[gene])

        # Retreives likelihood states
        effect1 = sigDig(data1[0][gene] / (10 ** scale))
        effect2 = sigDig(data2[0][gene]  / (10 ** scale))
        rat1 = sigDig(data1[1][gene])
        rat2 = sigDig(data2[1][gene])

        # Retrieves combo scores
        effect = sigDig(geneI[gene]  / (10 ** scale))
        rat = sigDig(geneL[gene])
        pRat = 'N/A'
        min_CI = sigDig(geneInterval[gene][0] / (10 ** scale))
        max_CI = sigDig(geneInterval[gene][1] / (10 ** scale))

        # Writes to file
        out_csv.writerow([geneID, name, info, comp, proc, fun,
                                num1, effect1, rat1,
                                num2, effect2, rat2,
                                effect, rat, pRat, min_CI, max_CI])


###############################################################################
# Saves record files

# Saves record for downstream analysis
with open(file_out + '_record.txt', 'w') as rec_open:
    rec_csv = csv.writer(rec_open, delimiter='\t')
    rec_csv.writerow(['analyzeCombo.py', current_version])
    rec_csv.writerow(['Date', time.strftime("%d:%m:%Y")])
    rec_csv.writerow(['Result file 1', args.res_file1])
    rec_csv.writerow(['Result file 2', args.res_file2])
    rec_csv.writerow(['Output file', file_out])
    rec_csv.writerow(['I step', I_step])
    rec_csv.writerow(['Scale', scale])
    rec_csv.writerow(['Number of processers', args.nums])
    
# Saves permanent record
with open(os.path.join('Records', 'analyzeCombo' + time.strftime("%Y%m%d%H%M%eS")) + '_record.txt', 'w') as back_open:
    rec_csv = csv.writer(back_open, delimiter='\t')
    rec_csv.writerow(['analyzeCombo.py', current_version])
    rec_csv.writerow(['Date', time.strftime("%d:%m:%Y")])
    rec_csv.writerow(['Result file 1', args.res_file1])
    rec_csv.writerow(['Result file 2', args.res_file2])
    rec_csv.writerow(['Output file', file_out])
    rec_csv.writerow(['I step', I_step])
    rec_csv.writerow(['Scale', scale])
    rec_csv.writerow(['Number of processers', args.nums])


###############################################################################
