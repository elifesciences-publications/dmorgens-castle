###############################################################################
# David Morgens
# 11/01/15
###############################################################################
# Import neccessary modules

'''
Function for aligning fastq files.
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

current_version = '0.7'


###############################################################################    

# Parses input using argparse module

# Initiates input parser
parser = argparse.ArgumentParser(description='Compares count files using casTLE')

# Non-optional arguments:
parser.add_argument('unt_file', help='File for untreated counts', type=str)

parser.add_argument('trt_file', help='File for treated counts', type=str)

parser.add_argument('name', help='Name for output files', type=str)

# Optional arguments:
parser.add_argument('-n', '--negative', dest='neg_name',
                help='Symbol used to denote negative controls',
                type=str)

parser.add_argument('-k', '--strength', dest='K', help='Normalizing constant',
                type=float, default=1.0)

parser.add_argument('-I', '--I_step', help='Step size in grid search',
                type=float, default=0.025)

parser.add_argument('-s', '--split', dest='split_mark', help='Delimiter for element name',
                type=str, default='_')

parser.add_argument('-t', '--threshhold', dest='thresh',
                help='Read cutoff for small count numbers', type=int, default=10)

parser.add_argument('-p', '--proccessors', dest='nums',
                help='Number of proccessors to use', type=int,
                default=20)

parser.add_argument('-r', '--reference',
                help='Location of reference files', type=str,
                default='GenRef')

parser.add_argument('-l', '--like', help='Likelihood function',
                default='casTLE',
                choices=['uniform', 'shifted', 'casTLE'])

parser.add_argument('-b', '--back', help='Background population for noise estimation',
                default='neg', choices=['all', 'neg', 'tar'])

parser.add_argument('-o', '--off_rate', help='Off-target rate',
                type=float, choices=[0.0, 0.1, 0.2, 0.3], default=0.0)

parser.add_argument('-z', '--zero_files', help='Time zero count files',
                nargs=2, type=str, default='')

parser.add_argument('-of', '--override_file', action='store_true',
                help='Overrides output to Results folder')

parser.add_argument('-x', '--exclude', type=str, help='Overrides output to Results folder')

# Saves all input to object args
args = parser.parse_args()


###############################################################################
# Processes and checks input arguments

# Short cut
split_mark = args.split_mark

# Checks if using single processor
if args.nums == 1:
    single = True
else:
    single = False

# Creates output file name
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

# Assigns default negative symbol for Cas9 or shRNA screens
if unt_screen_type in ['Cas9-4']:
    if not args.neg_name:
	args.neg_name = 'neg'

elif unt_screen_type in ['shRNA']:
    if not args.neg_name:
	args.neg_name = '0'

elif not args.neg_name:
    sys.exit('Screen type not recognized.\n' + 
		'Indicate negative elements -n or --negative')

# Assigns likelihood function based on input
if args.like == 'uniform':
    like_fun = likeshRNA
elif args.like == 'shifted':
    like_fun = casLike
elif args.like == 'casTLE':
    like_fun = likeEB
else:
    sys.exit('Likelihood function not recognized, use -l to change')


###############################################################################
# Pulls in gene info, IDs, as well as GO terms

print('Retrieving gene information')

# Retrieves gene information
geneID2Name, geneID2Info, geneName2ID, geneEns2Name = retrieveInfo()

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
		treated, args.neg_name, args.split_mark, args.K, time_zero)

print('Number of negative controls = ' + str(len(neg_rhos)))

if len(neg_rhos) == 0:
    sys.exit('No negative contols found.\n' + 
                'Change negative indicator with -n or --negative')


###############################################################################
# Calculates KS and MW p-values for each gene versus the negative controls

print('Calculating p-values')

# Selects background based on input options
if args.back == 'all':
    back_rhos = neg_rhos + tar_rhos
elif args.back == 'neg':
    back_rhos = neg_rhos
elif args.back == 'tar':
    back_rhos = tar_rhos
else:
    sys.exit('Unrecognized background option')

# Retrieves pvalues from auxilary function
rhoMW, rhoKS = retrievePvals(gene_rhos, neg_rhos, args.nums)


###############################################################################
# Determines grid search.

print('Determining search grid')

gene_span = {}

for gene, rhos in gene_rhos.items():

    # Maximum and minimum effect checked is twice the max and min of elements
    max_I = 2 * max(rhos + [0])
    min_I = 2 * min(rhos + [0])

    gene_span[gene] = (min_I, max_I)


###############################################################################
# Finds likelihoods for each gene

print('Calculating likelihoods')  

# Retrieves effect and log-likelihood ratio from auxilary function
geneI, geneL, geneInterval, geneDist = retrieveLikelihoods(gene_rhos, back_rhos,
                                        tar_rhos, args.off_rate, like_fun,
                                        args.nums, gene_span, args.I_step)


###############################################################################
# Writes output in human readable format

print('Outputing file')

with open(file_out + '.csv', 'w') as out_open:
    out_csv = csv.writer(out_open, delimiter=',', lineterminator='\n')

    # Writes a header
    out_csv.writerow(['#GeneID','Symbol','GeneInfo','Localization','Process',
                        'Function','Element #','Rho pMW','Rho pKS',
                        'Effect', 'L statistic', 'L p-value',
			'Min I', 'Max I', 'Elements'])

    for gene in gene_rhos:

	# Converts IDs as neccessary
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
        info = geneID2Info[geneID]
        comp = geneID2Comp[geneID]
        proc = geneID2Proc[geneID]
        fun = geneID2Fun[geneID]

        # Retrieves analysis of gene; rounding where appropriate
        num = len(gene_rhos[gene])
        pMW = sigDig(rhoMW[gene])
        pKS = sigDig(rhoKS[gene])

        # Retreives effect, log-likelihood ratio, and 95% likelihood interval
        max_e = sigDig(geneI[gene])
        rat = sigDig(geneL[gene])
        pRat = 'N/A'
        min_I = sigDig(geneInterval[gene][0])
        max_I = sigDig(geneInterval[gene][1])

        # Reformats guide reference as human readable
        elements = [str(sigDig(val)) + ' : ' + element for val,
                element in sorted(gene_ref[gene], key=lambda x : x[0])]

        # Writes to file
        out_csv.writerow([geneID, name, info, comp, proc, fun, num, pMW, pKS,
                                max_e, rat, pRat, min_I, max_I] + elements)


###############################################################################
# Saves record files

# Saves record for downstream analysis
with open(file_out + '_record.txt', 'w') as rec_open:
    rec_csv = csv.writer(rec_open, delimiter='\t')
    rec_csv.writerow(['analyzeCounts.py version', current_version])
    rec_csv.writerow(['Data/Time', time.strftime("%d:%H:%M:%eS")])
    rec_csv.writerow(['Untreated count file', args.unt_file])
    rec_csv.writerow(['Treated count file', args.trt_file])
    rec_csv.writerow(['Output file', file_out])
    rec_csv.writerow(['Count threshold', args.thresh])
    rec_csv.writerow(['Screen Type', trt_screen_type])
    rec_csv.writerow(['Selection strength', args.K])
    rec_csv.writerow(['Negative symbol', args.neg_name])
    rec_csv.writerow(['Split mark', args.split_mark])
    rec_csv.writerow(['Time zero files', args.zero_files])
    rec_csv.writerow(['Off-target rate', args.off_rate])
    rec_csv.writerow(['Likelihood function', like_fun.__name__])
    rec_csv.writerow(['I step', args.I_step])
    rec_csv.writerow(['Draw number', int(np.median(map(len, gene_rhos.values())))])
    rec_csv.writerow(['Background', args.back])
    rec_csv.writerow(['Number of processers', args.nums])
    rec_csv.writerow(['Total untreated counts', sum(untreated.values())])
    rec_csv.writerow(['Total treated counts', sum(treated.values())])
    rec_csv.writerow(['Missing/low untreated elements', belowUnt])
    rec_csv.writerow(['Missing/low treated elements', belowTrt])
    rec_csv.writerow(['Elements removed', removed])
    rec_csv.writerow(['Number of distinct elements', len(untreated)])

# Saves permanent record
with open(os.path.join('Records', 'analyzeCounts' + time.strftime("%d%H%M%S")), 'w') as back_open:
    rec_csv = csv.writer(back_open, delimiter='\t')
    rec_csv.writerow(['analyzeCounts.py version', current_version])
    rec_csv.writerow(['Data/Time', time.strftime("%d:%H:%M:%eS")])
    rec_csv.writerow(['Untreated count file', args.unt_file])
    rec_csv.writerow(['Treated count file', args.trt_file])
    rec_csv.writerow(['Output file', file_out])
    rec_csv.writerow(['Count threshold', args.thresh])
    rec_csv.writerow(['Screen Type', trt_screen_type])
    rec_csv.writerow(['Selection strength', args.K])
    rec_csv.writerow(['Negative symbol', args.neg_name])
    rec_csv.writerow(['Split mark', args.split_mark])
    rec_csv.writerow(['Time zero files', args.zero_files])
    rec_csv.writerow(['Off-target rate', args.off_rate])
    rec_csv.writerow(['Likelihood function', like_fun.__name__])
    rec_csv.writerow(['I step', args.I_step])
    rec_csv.writerow(['Draw number', int(np.median(map(len, gene_rhos.values())))])
    rec_csv.writerow(['Background', args.back])
    rec_csv.writerow(['Number of processers', args.nums])
    rec_csv.writerow(['Total untreated counts', sum(untreated.values())])
    rec_csv.writerow(['Total treated counts', sum(treated.values())])
    rec_csv.writerow(['Missing/low untreated elements', belowUnt])
    rec_csv.writerow(['Missing/low treated elements', belowTrt])
    rec_csv.writerow(['Elements removed', removed])
    rec_csv.writerow(['Number of distinct elements', len(untreated)])


###############################################################################
