###############################################################################
# David Morgens
# 04/06/2016
###############################################################################
# Imports neccessary modules

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import scipy.stats as st
import scipy.stats.mstats as ms
import csv
import sys
from collections import defaultdict
import argparse
from screenFun import *
import warnings
import math


###############################################################################    
# Version number

current_version = '1.0'


###############################################################################
# Parses input using argparse module

# Initiates input parser
parser = argparse.ArgumentParser(description='Visualizes replicate data')

# Non-optional arguments: The files containing results, as well as an output
parser.add_argument('res_file1', help='File for untreated results', type=str)

parser.add_argument('res_file2', help='File for treated results', type=str)

parser.add_argument('name', help='Name for output files', type=str)

# Optional arguments:
parser.add_argument('-hi', '--hist', action='store_false',
                        help='Flag to show data instead of 2D histogram')

parser.add_argument('-f', '--file_type', default='png',
                        help='File ending/type. Default is "png"')

# Arguments for statistics
parser.add_argument('-l', '--line', action='store_false',
                        help='Don\'t include linear regression line')

parser.add_argument('-t', '--title', action='store_false',
                        help='Don\'t include title')

# Arguments for labeling points
parser.add_argument('-n', '--names', help='List of genes to label', nargs='+', default=[])

parser.add_argument('-m', '--mouse', action='store_true',
                        help='Uses mouse gene information.')

# Arguments for axes
parser.add_argument('-xl', '--xlim', nargs=2, type=float,
                        help='x axis range')

parser.add_argument('-yl', '--ylim', nargs=2, type=float,
                        help='y axis range')

parser.add_argument('-x', '--x_axis', default='replicate 1',
                        help='Label for x axis. Default is "replicate 1"')

parser.add_argument('-y', '--y_axis', default='replicate 2',
                        help='Label for y axis. Default is "replicate 2"')

# Arguments for changing columns
parser.add_argument('-cr1', '--rat_col1', type=int,
                        help='Manual selection of significance column.')

parser.add_argument('-ce1', '--effect_col1', type=int,
                        help='Manual selection of effect column.')

parser.add_argument('-cr2', '--rat_col2', type=int,
                        help='Manual selection of significance column.')

parser.add_argument('-ce2', '--effect_col2', type=int,
                        help='Manual selection of effect column.')

parser.add_argument('-g', '--gene_col', type=int, default=1,
                        help='Manual selection of gene name column.')

# Override options
parser.add_argument('-o', '--record', action='store_true')

parser.add_argument('-of', '--override_file', action='store_true',
                        help='Override automatic targeting to Results folder')

# Saves all input to object args
args = parser.parse_args()


###############################################################################
# Checks arguments

if args.override_file:
    file_out = args.name
else:
    file_out = os.path.join('Results', args.name)

try:
    with open(file_out + '_effect.' + args.file_type, 'w') as out_open:
        pass

    os.remove(file_out + '_effect.' + args.file_type)

except IOError:
    sys.exit('Cannot write to output file:\n' + file_out + '\n'
                + 'Use -of or --override_file to change')

if os.path.isfile(args.res_file1) and os.path.isfile(args.res_file2):
    pass

else:
    sys.exit('Result files not found')


###############################################################################
# Finds count files

print('Retrieving records')

if not args.record:

    file_name1 = args.res_file1.split('.')[0]
    rec_file1 = file_name1 + '_record.txt'

    file_name2 = args.res_file2.split('.')[0]
    rec_file2 = file_name2 + '_record.txt'

    try:

        with open(rec_file1, 'r') as rec_open:
            rec_csv = csv.reader(rec_open, delimiter='\t')
            script1, version1 = rec_csv.next()

            if version1 != current_version:
                sys.exit('Version number not current\n'
                            + 'Rerun analysis')

        with open(rec_file2, 'r') as rec_open:
            rec_csv = csv.reader(rec_open, delimiter='\t')
            script2, version2 = rec_csv.next()

            if version2 != current_version:
                sys.exit('Version number not current\n'
                            + 'Rerun analysis')

    except IOError:
        sys.exit('Records not found')

    if version1 != version2:
        sys.exit('Versions not consistent\n'
                    + 'Rerun analysis')

else:
    script1 = 'NONE'
    script2 = 'NONE'
    print('Warning: record overridden')


###############################################################################
#

if args.rat_col1 and args.effect_col1:
    pass

elif script1 == 'analyzeCounts.py':

    if not args.rat_col1:
        args.rat_col1 = 8

    if not args.effect_col1:
        args.effect_col1 = 7

elif script1 == 'analyzeCombo.py':

    if not args.rat_col1:
        args.rat_col1 = 13

    if not args.effect_col1:
        args.effect_col1 = 12

else:
    sys.exit('Error: Result format not recognized')

if args.rat_col2 and args.effect_col2:
    pass

elif script2 == 'analyzeCounts.py':

    if not args.rat_col2:
        args.rat_col2 = 8

    if not args.effect_col2:
        args.effect_col2 = 7

elif script2 == 'analyzeCombo.py':

    if not args.rat_col2:
        args.rat_col2 = 13

    if not args.effect_col2:
        args.effect_col2 = 12

else:
    sys.exit('Error: Result format not recognized')


###############################################################################
# Pulls in gene info, IDs, as well as GO terms

print('Retrieving gene information')

# Uses different genetic information depending whether a human or mouse screen
geneID2Name, geneID2Info, geneName2ID, geneEns2Name = retrieveInfo(mouse=args.mouse)

# Retrieves GO information
geneID2Comp, geneID2Proc, geneID2Fun = retrieveGO()


###############################################################################
# Parses results file

print('Parsing results')

gene2effect1 = {}
gene2rat1 = {}

gene2effect2 = {}
gene2rat2 = {}

with open(args.res_file1, 'r') as res_open:

    res_csv = csv.reader(res_open, delimiter=',', lineterminator='\n')
    res_csv.next()

    for line in res_csv:
        gene = line[args.gene_col].upper()

        gene2effect1[gene] = float(line[args.effect_col1])
        gene2rat1[gene] = float(line[args.rat_col1])

with open(args.res_file2, 'r') as res_open:

    res_csv = csv.reader(res_open, delimiter=',', lineterminator='\n')
    res_csv.next()

    for line in res_csv:
        gene = line[args.gene_col].upper()

        gene2effect2[gene] = float(line[args.effect_col2])
        gene2rat2[gene] = float(line[args.rat_col2])


###############################################################################
# Pairs genes

gene2effects = {}
gene2rats = {}

for gene in gene2effect1:

    effect1 = gene2effect1[gene]
    rat1 = gene2rat1[gene]

    geneID, name, entrez = retrieveIDs(gene, geneID2Name, geneName2ID, geneEns2Name)

    if gene in gene2effect2:
        gene2 = gene

    elif geneID in gene2effect2:
        gene2 = geneID

    elif name in gene2effect2:
        gene2 = name

    elif entrez in gene2effect2:
        gene2 = entrez

    else:
        continue

    effect2 = gene2effect2[gene2]
    rat2 = gene2rat2[gene2]

    gene2effects[gene] = (effect1, effect2)
    gene2rats[gene] = (math.copysign(rat1, effect1), math.copysign(rat2, effect2))


###############################################################################
# Finds labels

gene_labels = []

if args.names:
    for gene in args.names:

        geneID, name, entrez = retrieveIDs(gene, geneID2Name, geneName2ID, geneEns2Name)

        if gene.upper() in gene2effects:
            gene_labels.append((gene, gene.upper()))

        elif geneID in gene2effects:
            gene_labels.append((gene, geneID))

        elif name in gene2effects:
            gene_labels.append((gene, name))

        elif entrez in gene2effects:
            gene_labels.append((gene, entrez))

        else:
            print('Warning: ' + gene + ' not found')


###############################################################################
# Plots gene against gene using signed ratio statistic

print('Plotting figures')

x = []
y = []

for gene in gene2rats:

    rat1, rat2 = gene2rats[gene]

    x.append(rat1)
    y.append(rat2)

plt.figure(dpi=400)

if args.hist:
    plt.hist2d(x, y, (100, 100), cmap=cm.jet, norm=LogNorm())

else:
    plt.plot(x, y, '.', color='0', alpha=0.25, markersize=10, markeredgewidth=0)

plt.xlabel('casTLE score for ' + args.x_axis)
plt.ylabel('casTLE score for ' + args.y_axis)

slope, intercept, r_value, p_value, std_err = st.linregress(x, y)

if args.ylim:
    miny, maxy = args.ylim
    plt.ylim([miny, maxy])

if args.xlim:
    minx, maxx = args.xlim
    plt.xlim([minx, maxx])

if args.line:
    xgrid = np.linspace(min(x), max(x), 100)
    plt.plot(xgrid, slope * xgrid + intercept,'-', color='0')

if args.names:
    for name, gene in gene_labels:
	rat1, rat2 = gene2rats[gene]
	plt.text(rat1, rat2, name)

if args.title:
    plt.title('Reproducibility of casTLE score by gene; R2 = '
                + str(round(r_value ** 2, 2)))

plt.savefig(file_out + '_score.' + args.file_type)
plt.close()


###############################################################################
# Plots gene against gene using estimated effect size

x = []
y = []

for gene in gene2effects:

    effect1, effect2 = gene2effects[gene]

    x.append(effect1)
    y.append(effect2)

plt.figure(dpi=400)

plt.plot(x, y, '.', color='0', alpha=0.25, markersize=10, markeredgewidth=0)
plt.xlabel('casTLE effect estimate for ' + args.x_axis)
plt.ylabel('casTLE effect estimate for ' + args.y_axis)

slope, intercept, r_value, p_value, std_err = st.linregress(x, y)

if args.line:
    xgrid = np.linspace(min(x), max(x), 100)
    plt.plot(xgrid, slope * xgrid + intercept,'-', color='0')

if args.names:
    for name, gene in gene_labels:
	effect1, effect2 = gene2effects[gene]
	plt.text(effect1, effect2, name)

if args.title:
    plt.title('Reproducibility of casTLE effect by gene; R2 = '
                + str(round(r_value ** 2, 2)))

plt.savefig(file_out + '_effect.' + args.file_type)
plt.close()

print('Finished')


###############################################################################
