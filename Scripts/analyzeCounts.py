###############################################################################
# David Morgens
# 04/06/2016
###############################################################################
# Import neccessary modules

'''
Function for comparing count files.
'''

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
parser = argparse.ArgumentParser(description='Compares count files using casTLE')

# Non-optional arguments:
parser.add_argument('unt_file', help='File for untreated counts', type=str)

parser.add_argument('trt_file', help='File for treated counts', type=str)

parser.add_argument('name', help='Name for output files', type=str)

# Options for element IDs
parser.add_argument('-n', '--neg_name',
                help='Symbol used to denote negative controls. Default is 0.',
                type=str, default='0')

parser.add_argument('-s', '--split', dest='split_mark',
                        help='Delimiter for element name. Default is _',
                type=str, default='_')

parser.add_argument('-x', '--exclude', type=str, nargs='+',
                        help='Only include elements containing substrings.')

# Optional arguments for how analysis is performed:
parser.add_argument('-t', '--threshhold', dest='thresh',
                help='Read cutoff for small count numbers. Default is 10.',
                        type=int, default=10)

parser.add_argument('-k', '--strength', dest='K',
                        help='Normalizing constant. Default is 1.',
                        type=float, default=1.0)

parser.add_argument('-b', '--back',
                help='Background population for noise estimation. Default is neg.',
                default='neg', choices=['all', 'neg', 'tar'])

parser.add_argument('-z', '--zero_files',
                help='Time zero count files. Optional.',
                nargs=2, type=str, default='')

# Options for how computation is performed
parser.add_argument('-c', '--scale', type=int, default=3,
                help='Scale of calculations; default is 3')

parser.add_argument('-I', '--I_step', type=float, default=0.1,
                help='Step size in grid search; default is 0.1')

parser.add_argument('-p', '--proccessors', dest='nums',
                help='Number of proccessors to use; default is 20', type=int,
                default=20)

# Options for overriding behavior
parser.add_argument('-r', '--reference',
                help='Location of reference files; default is GeneRef', type=str,
                default='GenRef')

parser.add_argument('-of', '--override_file', action='store_true',
                help='Overrides restriction of output to Results folder')

parser.add_argument('-m', '--mouse', action='store_true',
                help='Uses mouse gene information.')

parser.add_argument('-ro', '--record', action='store_false',
                help='Allows script to run without record of count files.')

# Saves all input to object args
args = parser.parse_args()


###############################################################################
# Processes and checks input arguments

# Creates output file name
if args.override_file:
    file_out = args.name
else:
    file_out = os.path.join('Results', args.name)

# Checks if can write to output
try:
    with open(file_out + '_record.txt', 'w') as out_open:
        pass

    os.remove(file_out + '_record.txt')

except:
    sys.exit('Cannot write to output file:\n' + file_out + '\n'
                + 'Use -of or --override_file to change')


###############################################################################
# Finds record files

if args.record:

    print('Retrieving records')

    # Retrieves record files for count files
    unt_rec_name = args.unt_file[: -11]
    unt_rec_file = unt_rec_name + '_record.txt'

    trt_rec_name = args.trt_file[: -11]
    trt_rec_file = trt_rec_name + '_record.txt'

    try:
        # Parses record files
        with open(unt_rec_file, 'r') as rec_open:
            rec_csv = csv.reader(rec_open, delimiter='\t')
            unt_version = rec_csv.next()[1]
            unt_time1 = rec_csv.next()[1]
            unt_seq_file = rec_csv.next()[1]
            unt_seq_add_file = rec_csv.next()[1]
            unt_out = rec_csv.next()[1]
            unt_screen_type = rec_csv.next()[1]

        with open(trt_rec_file, 'r') as rec_open:
            rec_csv = csv.reader(rec_open, delimiter='\t')
            trt_version = rec_csv.next()[1]
            trt_time1 = rec_csv.next()[1]
            trt_seq_file = rec_csv.next()[1]
            trt_seq_add_file = rec_csv.next()[1]
            trt_out = rec_csv.next()[1]
            trt_screen_type = rec_csv.next()[1]

    except IOError:

        print(unt_rec_file)
        print(trt_rec_file)

        sys.exit('Record of count file not found\n'
                + 'Change file name or rerun makeCounts.py')

    # Checks comparison makes sense
    if unt_screen_type != trt_screen_type:
        sys.exit('Screen types do not match.')

else:

    print('Warning: Record file overriden by -ro flag')

    unt_version = 'Overridden'
    unt_time1 = 'Overridden'
    unt_seq_file = 'Overridden'
    unt_seq_add_file = 'Overridden'
    unt_out = 'Overridden'
    unt_screen_type = 'Overridden'

    trt_version = 'Overridden'
    trt_time1 = 'Overridden'
    trt_seq_file = 'Overridden'
    trt_seq_add_file = 'Overridden'
    trt_out = 'Overridden'
    trt_screen_type = 'Overridden'


###############################################################################
# Pulls in gene info, IDs, as well as GO terms

print('Retrieving gene information')

# Retrieves gene information
geneID2Name, geneID2Info, geneName2ID, geneEns2Name = retrieveInfo(mouse=args.mouse)

# Retrieves GO information
geneID2Comp, geneID2Proc, geneID2Fun = retrieveGO()


###############################################################################
# Pulls in untreated and treated counts and filters by defined threshold

print('Filtering reads')

# Retrieves filtered counts for auxilary function
untreated, treated, stats, time_zero = filterCounts(args.unt_file,
                                                args.trt_file, args.thresh,
                                                args.zero_files, args.exclude)

# Outputs statistics
belowTrt, belowUnt, removed = stats

print('Total untreated counts: ' + str(sum(untreated.values())))
print('Total treated counts: ' + str(sum(treated.values())))
print('Missing/low untreated elements: ' + str(belowUnt))
print('Missing/low treated elements: ' + str(belowTrt))
print('Elements removed: ' + str(removed))
print('Number of distinct elements: ' + str(len(untreated)))


###############################################################################
# Calculates enrichment values

print('Calculating enrichment values')

# Retrieves enrichment values from auxilary function
element_rhos, gene_rhos, neg_rhos, tar_rhos, gene_ref = enrich_all(untreated,
		treated, args.neg_name, args.split_mark, args.K, time_zero, args.back)

print('Number of negative controls = ' + str(len(neg_rhos)))

if args.back == 'neg' and len(neg_rhos) == 0:
    sys.exit('No negative contols found.\n' + 
                'Change negative indicator with -n or --negative')


###############################################################################
# Selects background based on input options

if args.back == 'all':
    back_rhos = neg_rhos + tar_rhos

elif args.back == 'neg':
    back_rhos = neg_rhos

elif args.back == 'tar':
    back_rhos = tar_rhos

else:
    sys.exit('Unrecognized background option')


###############################################################################
# Determines grid search.

print('Determining search grid')

gene_span = {}

for gene, rhos in gene_rhos.items():

    # Maximum and minimum effect checked is twice the max and min of elements
    max_I = int(2 * max(rhos + [-1, 0, 1]))
    min_I = int(2 * min(rhos + [-1, 0, 1]))

    gene_span[gene] = (min_I, max_I)


###############################################################################
# Finds likelihoods for each gene

print('Calculating likelihoods')

# Retrieves effect and log-likelihood ratio from auxilary function
geneI, geneL, geneInterval, geneDist = retrieveLikelihoods(gene_rhos, back_rhos,
                                                                args.nums, gene_span,
                                                                args.scale, args.I_step)


###############################################################################
# Writes output in human readable format

print('Outputing file')

with open(file_out + '.csv', 'w') as out_open:
    out_csv = csv.writer(out_open, delimiter=',', lineterminator='\n')

    # Writes a header
    out_csv.writerow(['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process',
                        'Function', 'Element #',
                        'casTLE Effect', 'casTLE Score', 'casTLE p-value',
			'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements'])

    for gene in gene_rhos:

	# Converts IDs as neccessary
        geneID, name, entrez = retrieveIDs(gene, geneID2Name, geneName2ID, geneEns2Name)

        # Retrieves information about gene
        info = geneID2Info[entrez]
        comp = geneID2Comp[entrez]
        proc = geneID2Proc[entrez]
        fun = geneID2Fun[entrez]

        # Retrieves analysis of gene; rounding where appropriate
        num = len(gene_rhos[gene])

        # Retreives effect, log-likelihood ratio, and 95% likelihood interval
        max_e = sigDig(float(geneI[gene]) / (10 ** args.scale))
        rat = sigDig(geneL[gene])
        pRat = 'N/A'

        min_I = sigDig(float(geneInterval[gene][0]) / (10 ** args.scale))
        max_I = sigDig(float(geneInterval[gene][1]) / (10 ** args.scale))

        # Reformats guide reference as human readable
        elements = [str(sigDig(val)) + ' : ' + element.split(args.split_mark)[-1] for val,
                    element in sorted(gene_ref[gene], key=lambda x : x[0])]

        # Writes to file
        out_csv.writerow([geneID, name, info, comp, proc, fun, num,
                                max_e, rat, pRat, min_I, max_I] + elements)


###############################################################################
# Saves record files

# Saves record for downstream analysis
with open(file_out + '_record.txt', 'w') as rec_open:
    rec_csv = csv.writer(rec_open, delimiter='\t')
    rec_csv.writerow(['analyzeCounts.py', current_version])
    rec_csv.writerow(['Date', time.strftime("%d:%m:%Y")])
    rec_csv.writerow(['Untreated count file', args.unt_file])
    rec_csv.writerow(['Treated count file', args.trt_file])
    rec_csv.writerow(['Time zero files', args.zero_files])
    rec_csv.writerow(['Output file', file_out])
    rec_csv.writerow(['Screen Type', trt_screen_type])
    rec_csv.writerow(['Negative symbol', args.neg_name])
    rec_csv.writerow(['Split mark', args.split_mark])
    rec_csv.writerow(['Excluded strings', args.exclude])
    rec_csv.writerow(['Count threshold', args.thresh])
    rec_csv.writerow(['Selection strength', args.K])
    rec_csv.writerow(['Background', args.back])
    rec_csv.writerow(['Effect step', args.I_step])
    rec_csv.writerow(['Scale', args.scale])
    rec_csv.writerow(['Draw number', int(np.median(map(len, gene_rhos.values())))])
    rec_csv.writerow(['Number of processers', args.nums])
    rec_csv.writerow(['Total untreated counts', sum(untreated.values())])
    rec_csv.writerow(['Total treated counts', sum(treated.values())])
    rec_csv.writerow(['Missing/low untreated elements', belowUnt])
    rec_csv.writerow(['Missing/low treated elements', belowTrt])
    rec_csv.writerow(['Elements removed', removed])
    rec_csv.writerow(['Number of distinct elements', len(untreated)])

# Saves permanent record
with open(os.path.join('Records', 'analyzeCounts' + time.strftime("%Y%m%d%H%M%eS") + '_record.txt'), 'w') as back_open:
    rec_csv = csv.writer(back_open, delimiter='\t')
    rec_csv.writerow(['analyzeCounts.py version', current_version])
    rec_csv.writerow(['Date', time.strftime("%d:%m:%Y")])
    rec_csv.writerow(['Untreated count file', args.unt_file])
    rec_csv.writerow(['Treated count file', args.trt_file])
    rec_csv.writerow(['Time zero files', args.zero_files])
    rec_csv.writerow(['Output file', file_out])
    rec_csv.writerow(['Screen Type', trt_screen_type])
    rec_csv.writerow(['Negative symbol', args.neg_name])
    rec_csv.writerow(['Split mark', args.split_mark])
    rec_csv.writerow(['Excluded strings', args.exclude])
    rec_csv.writerow(['Count threshold', args.thresh])
    rec_csv.writerow(['Selection strength', args.K])
    rec_csv.writerow(['Background', args.back])
    rec_csv.writerow(['Effect step', args.I_step])
    rec_csv.writerow(['Scale', args.scale])
    rec_csv.writerow(['Draw number', int(np.median(map(len, gene_rhos.values())))])
    rec_csv.writerow(['Number of processers', args.nums])
    rec_csv.writerow(['Total untreated counts', sum(untreated.values())])
    rec_csv.writerow(['Total treated counts', sum(treated.values())])
    rec_csv.writerow(['Missing/low untreated elements', belowUnt])
    rec_csv.writerow(['Missing/low treated elements', belowTrt])
    rec_csv.writerow(['Elements removed', removed])
    rec_csv.writerow(['Number of distinct elements', len(untreated)])


###############################################################################
