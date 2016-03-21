###############################################################################

from __future__ import division
import csv
import time
import argparse
from screenFun import *
import sys


###############################################################################    
# Version number

current_version = '0.7'


###############################################################################    

# Parses input using argparse module

# Initiates input parser
parser = argparse.ArgumentParser(description='Adds permutations for pvalues')

# Non-optional arguments: The results file and the number of permutations
parser.add_argument('res_file', help='Results file', type=str)

parser.add_argument('perm_num', help='Number of permutations', type=int)

# Optional arguments
parser.add_argument('-p', '--proccessors', dest='nums',
                help='Number of processors to use', type=int,
                default=20)

parser.add_argument('-e', '--erase', action='store_true',
		help='Don\t keep previous permutations')

parser.add_argument('-t', '--out_time', action='store_true',
		help='Output timestamps')

parser.add_argument('-r', '--ratio_col', default=13, type=int,
                        help='Column containing ratio scores')

# Saves all input to object args
args = parser.parse_args()


###############################################################################
# Finds records

print('Retrieving records')

# Finds record file
name = args.res_file[: -4]
rec_file = name + '_record.txt'

try:
    # Parses record file
    with open(rec_file, 'r') as rec_open:
        rec_csv = csv.reader(rec_open, delimiter='\t')
        version = rec_csv.next()[1]

        if version != current_version:
            sys.exit('Version number not current\n'
                        + 'Rerun analysis')

        last_time = rec_csv.next()[1]
        res_file1 = rec_csv.next()[1]
        res_file2 = rec_csv.next()[1]
        unt_file1 = rec_csv.next()[1]
        unt_file2 = rec_csv.next()[1]
        trt_file1 = rec_csv.next()[1]
        trt_file2 = rec_csv.next()[1]
        file_out = rec_csv.next()[1]
        thresh1 = int(rec_csv.next()[1])
        thresh2 = int(rec_csv.next()[1])
        screen_type1 = rec_csv.next()[1]
        screen_type2 = rec_csv.next()[1]
        K1 = float(rec_csv.next()[1])
        K2 = float(rec_csv.next()[1])
        neg_name1 = rec_csv.next()[1]
        neg_name2 = rec_csv.next()[1]
        split_mark1 = rec_csv.next()[1]
        split_mark2 = rec_csv.next()[1]
        zero_files1 = rec_csv.next()[1]
        zero_files2 = rec_csv.next()[1]
        off_rate1 = float(rec_csv.next()[1])
        off_rate2 = float(rec_csv.next()[1])
        like_fun1 = eval(rec_csv.next()[1])
        like_fun2 = eval(rec_csv.next()[1])
        draw_num1 = int(rec_csv.next()[1])
        draw_num2 = int(rec_csv.next()[1])
        I_step = float(rec_csv.next()[1])

except IOError:
    sys.exit('Record of result file not found\n'
                + 'Change file name or rerun analysis')


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
                                                trt_file1, thresh1, zero_files1)

untreated2, treated2, stats2, time_zero2 = filterCounts(unt_file2,
                                                trt_file2, thresh2, zero_files2)


###############################################################################
# Calculates enrichment values

print('Calculating enrichment values')

element_rhos1, gene_rhos1, neg_rhos1, tar_rhos1, gene_ref1 = enrich_all(untreated1,
                                treated1, neg_name1, split_mark1, K1, time_zero1)

element_rhos2, gene_rhos2, neg_rhos2, tar_rhos2, gene_ref2 = enrich_all(untreated2,
                                treated2, neg_name2, split_mark2, K2, time_zero2)


###############################################################################
# Calculates null distribution of ratio statistic

print('Running permutations')

gene2line = {}
gene2rat = {}

# Reads in previously calculated ratios
with open(args.res_file, 'r') as res_open:

    res_csv = csv.reader(res_open, delimiter=',', lineterminator='\n')
    header = res_csv.next()

    for line in res_csv:
        gene2line[line[0]] = line
        gene2rat[line[0]] = float(line[args.ratio_col])

ref_file = name + '_ref.csv'

if args.out_time:
    print time.strftime("%d:%H:%M:%S")

# Retrieves pvalues from auxilary function
geneP, all_perm_num = comboPerm(draw_num1, draw_num2, args.perm_num,
                                neg_rhos1, neg_rhos2,
                                tar_rhos1, tar_rhos2,
                                off_rate1, off_rate2,
                                like_fun1, like_fun2,
                                args.nums, ref_file,
                                gene2rat, I_step, args.erase)

if args.out_time:
    print time.strftime("%d:%H:%M:%S")

print('Permutations used: ' + all_perm_num)


###############################################################################
# Writes output in human readable format

print('Outputing file')

with open(args.res_file,'w') as out_open:
    out_csv = csv.writer(out_open, delimiter=',', lineterminator = '\n')

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
    rec_csv.writerow(['addCombo.py version', '0.7'])
    rec_csv.writerow(['Data/Time', time.strftime("%d:%H:%M:%S")])
    rec_csv.writerow(['Number of permutations', all_perm_num])


###############################################################################
