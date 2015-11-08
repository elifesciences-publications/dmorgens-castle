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
                help='Number of prcessors to use', type=int,
                default=20)

parser.add_argument('-e', '--erase', action='store_true',
		help='Don\t keep previous permutations')

parser.add_argument('-t', '--out_time', action='store_true',
		help='Output timestamps')

parser.add_argument('-r', '--ratio_col', default=10, type=int,
                        help='Column containing ratio scores')

# Saves all input to object args
args = parser.parse_args()

if args.nums == 1:
    single = True
else:
    single = False


###############################################################################
# Finds records

print('Retrieving records')

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
        unt_file = rec_csv.next()[1]
        trt_file = rec_csv.next()[1]
        file_out = rec_csv.next()[1]
        thresh = int(rec_csv.next()[1])
        screen_type = rec_csv.next()[1]
        K = float(rec_csv.next()[1])
        neg_name = rec_csv.next()[1]
        split_mark = rec_csv.next()[1]
        zero_files = rec_csv.next()[1]
        off_rate = float(rec_csv.next()[1])
        like_fun = eval(rec_csv.next()[1])
        I_step = float(rec_csv.next()[1])
        draw_num = int(rec_csv.next()[1])
        back = rec_csv.next()[1]

except IOError:
    sys.exit('Record of result file not found\n'
                + 'Change file name or rerun analysis')

print('Drawing ' + str(draw_num) + ' random elements')


###############################################################################
# Pulls in gene info, IDs, as well as GO terms

print('Retrieving gene information')

# Uses different genetic information depending whether a human or mouse screen
geneID2Name, geneID2Info, geneName2ID, geneEns2Name = retrieveInfo()

# Retrieves GO information
geneID2Comp, geneID2Proc, geneID2Fun = retrieveGO()


###############################################################################
# Pulls in untreated and treated counts, and filters by defined threshold

print('Filtering reads')

# Retreives filtered counts for auxilary function
untreated, treated, stats, time_zero = filterCounts(unt_file,
                                                trt_file, thresh, zero_files)


###############################################################################
# Calculates enrichment values

print('Calculating enrichment values')

shRNA_rhos, gene_rhos, neg_rhos, tar_rhos, gene_ref = enrich_all(untreated,
                                treated, neg_name, split_mark, K, time_zero)

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

gene2line = {}
gene2rat = {}

# Reads in previously calculated ratios
with open(args.res_file, 'r') as res_open:
    res_csv = csv.reader(res_open, delimiter=',', lineterminator = '\n')
    header = res_csv.next()
    for line in res_csv:
        gene2line[line[1]] = line
        gene2rat[line[1]] = float(line[args.ratio_col])

ref_file = name + '_ref.csv'

if args.out_time:
    print time.strftime("%d:%H:%M:%S")

# Retrieves pvalues from auxilary function
geneP, all_perm_num = retrievePerm(draw_num, args.perm_num, back_rhos, tar_rhos,
                off_rate, like_fun, args.nums, ref_file, gene2rat, I_step, args.erase)

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
    rec_csv.writerow(['addPermutations.py version', '0.7'])
    rec_csv.writerow(['Data/Time', time.strftime("%d:%H:%M:%S")])
    rec_csv.writerow(['Number of permutations', all_perm_num])


###############################################################################
