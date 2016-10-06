###############################################################################
# David Morgens
# 04/06/2016
###############################################################################

from __future__ import division
import csv
import time
import argparse
from screenFun import *
import sys


###############################################################################    
# Version number

current_version = '1.0'


###############################################################################    

# Parses input using argparse module

# Initiates input parser
parser = argparse.ArgumentParser(description='Adds permutations to calculate pvalues')

# Non-optional arguments: The results file and the number of permutations
parser.add_argument('res_file', help='Results file location', type=str)

parser.add_argument('perm_num', help='Number of permutations', type=int)

# Optional arguments
parser.add_argument('-p', '--proccessors', dest='nums',
                help='Number of processors to use. Default is 20.', type=int,
                default=20)

parser.add_argument('-e', '--erase', action='store_true',
		help='Don\t keep previous permutations.')

parser.add_argument('-r', '--ratio_col', default=13, type=int,
                        help='Column containing ratio scores. Default is 13.')

parser.add_argument('-m', '--mouse', action='store_true',
                        help='Flag for mouse screen')

# Saves all input to object args
args = parser.parse_args()


###############################################################################
# Finds records

print('Retrieving records')

name = args.res_file[: -4]
rec_file = name + '_record.txt'

try:
    # Parses record file
    with open(rec_file, 'r') as rec_open:
        rec_csv = csv.reader(rec_open, delimiter='\t')
        script, version = rec_csv.next()

        if version != current_version:
            sys.exit('Error: Version number not current\n'
                        + 'Rerun analysis')

        if script != 'analyzeCombo.py':
            sys.exit('Error: File is not a combo file')

        prev_time = rec_csv.next()[1]
        res_file1 = rec_csv.next()[1]
        res_file2 = rec_csv.next()[1]
        file_out = rec_csv.next()[1]
        I_step = float(rec_csv.next()[1])
        scale = int(rec_csv.next()[1])

except IOError:
    sys.exit('Record of result file not found\n'
                + 'Change file name or rerun analysis')

except StopIteration:
    sys.exit('Empty record file found\n'
                + 'Rerun analysis')

stats1, files1, info1, param1 = retrieveRecord(res_file1, current_version)
stats2, files2, info2, param2 = retrieveRecord(res_file2, current_version)

script1, version1, last_time1 = stats1
unt_file1, trt_file1, zero_files1, file_out1 = files1
screen_type1, neg_name1, split_mark1, exclude1 = info1
thresh1, K1, back1, I_step1, scale1, draw_num1 = param1

script2, version2, last_time2 = stats2
unt_file2, trt_file2, zero_files2, file_out2 = files2
screen_type2, neg_name2, split_mark2, exclude2 = info2
thresh2, K2, back2, I_step2, scale2, draw_num2 = param2


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
                                                trt_file1, thresh1, zero_files1,
                                                exclude1)

untreated2, treated2, stats2, time_zero2 = filterCounts(unt_file2,
                                                trt_file2, thresh2, zero_files2,
                                                exclude2)


###############################################################################
# Calculates enrichment values

print('Calculating enrichment values')

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
# Calculates null distribution of ratio statistic

permI, permL, permInterval = comboPerm(draw_num1, draw_num2, args.perm_num,
                                        back_rhos1, back_rhos2,
                                        tar_rhos1, tar_rhos2,
                                        args.nums, I_step, scale)


###############################################################################
# Access previously calculate ratios and permutations

print('Calculating p-values')

gene2line, geneP, header, all_perm_num = calculatePval(args.res_file, permI, permL,
                                        args.erase, args.ratio_col)

print('Permutations used: ' + all_perm_num)


###############################################################################
# Writes output in human readable format

print('Outputing file')

with open(args.res_file, 'w') as out_open:
    out_csv = csv.writer(out_open, delimiter=',', lineterminator='\n')

    # Writes a header
    out_csv.writerow(header)

    for geneID in gene2line:

        line = gene2line[geneID]
        line[args.ratio_col + 1] = sigDig(geneP[geneID])

        # Writes to file
        out_csv.writerow(line)


###############################################################################
# Appends note to existing record

with open(rec_file, 'a') as rec_open:
    rec_csv = csv.writer(rec_open, delimiter='\t')
    rec_csv.writerow(['addCombo.py', current_version])
    rec_csv.writerow(['Date', time.strftime("%d:%m:%Y")])
    rec_csv.writerow(['Number of permutations', all_perm_num])


###############################################################################
