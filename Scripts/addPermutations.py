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
parser = argparse.ArgumentParser(description='Adds permutations for pvalues')

# Non-optional arguments: The results file and the number of permutations
parser.add_argument('res_file', help='Results file', type=str)

parser.add_argument('perm_num', help='Number of permutations', type=int)

# Optional arguments
parser.add_argument('-p', '--proccessors', dest='nums',
                help='Number of prcessors to use', type=int,
                default=20)

parser.add_argument('-e', '--erase', action='store_true',
		help='Don\t keep previous permutations')

parser.add_argument('-t', '--out_time', action='store_true',
		help='Output timestamps.')

parser.add_argument('-r', '--ratio_col', default=8, type=int,
                        help='Column containing ratio scores')

parser.add_argument('-m', '--mouse', action='store_true',
                        help='Uses mouse gene information.')

# Saves all input to object args
args = parser.parse_args()


###############################################################################
# Finds records

print('Retrieving records')

stats, files, info, param = retrieveRecord(args.res_file, current_version)

unt_file, trt_file, zero_files, file_out = files
screen_type, neg_name, split_mark, exclude = info
thresh, K, back, I_step, scale, draw_num = param

print('Drawing ' + str(draw_num) + ' random elements')


###############################################################################
# Pulls in gene info, IDs, as well as GO terms

print('Retrieving gene information')

# Uses different genetic information depending whether a human or mouse screen
geneID2Name, geneID2Info, geneName2ID, geneEns2Name = retrieveInfo(mouse=args.mouse)

# Retrieves GO information
geneID2Comp, geneID2Proc, geneID2Fun = retrieveGO()


###############################################################################
# Pulls in untreated and treated counts, and filters by defined threshold

print('Filtering reads')

# Retreives filtered counts for auxilary function
untreated, treated, stats, time_zero = filterCounts(unt_file,
                                                trt_file, thresh,
                                                zero_files, exclude)


###############################################################################
# Calculates enrichment values

print('Calculating enrichment values')

element_rhos, gene_rhos, neg_rhos, tar_rhos, gene_ref = enrich_all(untreated,
		treated, neg_name, split_mark, K, time_zero, back)

# Chooses appropriate background from record file
if back == 'all':
    back_rhos = neg_rhos + tar_rhos
elif back == 'neg':
    back_rhos = neg_rhos
elif back == 'tar':
    back_rhos = tar_rhos
else:
    sys.exit('Unrecognized background option')


###############################################################################
# Calculates null distribution of ratio statistic

print('Running permutations')

if args.out_time:
    print time.strftime("%d:%H:%M:%S")

# Retrieves pvalues from auxilary function
permI, permL, permInterval = retrievePerm(draw_num, args.perm_num, back_rhos,
                                                tar_rhos, args.nums, scale, I_step)

if args.out_time:
    print time.strftime("%d:%H:%M:%S")


###############################################################################
# Access previously calculate ratios and permutations

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

name = args.res_file[: -4]
rec_file = name + '_record.txt'

with open(rec_file, 'a') as rec_open:
    rec_csv = csv.writer(rec_open, delimiter='\t')
    rec_csv.writerow(['addPermutations.py', current_version])
    rec_csv.writerow(['Date', time.strftime("%d:%m:%Y")])
    rec_csv.writerow(['Number of permutations', all_perm_num])


###############################################################################
