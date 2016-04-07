###############################################################################
# David Morgens
# 04/06/2016
###############################################################################
# Imports neccessary modules

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import scipy.stats as st
import scipy.stats.mstats as ms
import csv
import sys
from collections import defaultdict
import time
import argparse
from screenFun import *


###############################################################################    
# Version number

current_version = '1.0'


###############################################################################
# Parses input using argparse module

# Initiates input parser
parser = argparse.ArgumentParser(description='Make gene plots')

# Non-optional arguments:
parser.add_argument('res_file', help='Result file', type=str)

parser.add_argument('gene_list', help='List of genes', type=str, nargs='+')

# Optional arguments to change images
parser.add_argument('-c', '--cloud', action='store_false',
                        help='Don\'t generate cloud plot')

parser.add_argument('-hi', '--hist', action='store_false',
                        help='Don\'t generate histogram plot')

parser.add_argument('-sl', '--legend', help='Include legend',
                        action='store_true')

# Optional arguments for output
parser.add_argument('-of', '--override_file', type=str,
                        help='Override automatic targeting to Results folder')

parser.add_argument('-f', '--file_type', default='png',
                        help='File ending/type. Default is "png"')

parser.add_argument('-m', '--mouse', action='store_true',
                        help='Uses mouse gene information.')

# Optional arguments for axes
parser.add_argument('-y', '--ylim', type=float,
                        help='y axis maximum.')

parser.add_argument('-x', '--xlim', type=float,
                        help='x axis range.')

# Optional arguments for analysis

parser.add_argument('-e', '--effect_col', type=int, default=7,
                        help='Manual selection of effect size column')

parser.add_argument('-t', '--thresh', type=int,
                        help='Override threshold for counts')

# Saves all input to object args
args = parser.parse_args()


###############################################################################
# Finds records file

print('Retrieving records')

stats, files, info, param = retrieveRecord(args.res_file, current_version)

unt_file, trt_file, zero_files, file_out = files
screen_type, neg_name, split_mark, exclude = info
thresh, K, back, I_step, scale, draw_num = param

file_name = args.res_file.split('.')[0]
rec_file = file_name + '_record.txt'

if args.override_file:
    file_out = args.override_file

else:
    file_out = file_name

if args.thresh:
    thresh = args.thresh


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
# Calculates where zero is by considering the median of the negative control
# guides.

print('Calculating enrichment values')
    
element_rhos, gene_rhos, neg_rhos, tar_rhos, gene_ref = enrich_all(untreated,
		treated, neg_name, split_mark, K, time_zero, back)

for gene, ref in gene_ref.items():

    geneID, name, entrez = retrieveIDs(gene, geneID2Name, geneName2ID, geneEns2Name)

    gene_ref[geneID] = ref
    gene_ref[name] = ref
    gene_ref[entrez] = ref


###############################################################################
# Plots genes given by genelist

print('Retrieving results')

# Collects the count values for all negative guides
neg_x = []
neg_y = []

all_x = []
all_y = []

# Collects the negative control counts   
for element in untreated:

    count_x = math.log(untreated[element], 2)
    count_y = math.log(treated[element], 2)

    all_x.append(count_x)
    all_y.append(count_y)

    if element.split(split_mark)[0].startswith(neg_name):
        neg_x.append(count_x)
        neg_y.append(count_y)

if back == 'neg':

   back_x = neg_x
   back_y = neg_y

   negDist = st.gaussian_kde(neg_rhos)

elif back == 'tar' or back == 'all':

   back_x = all_x
   back_y = all_y

   negDist = st.gaussian_kde(tar_rhos)


###############################################################################
#

gene2effect = {}

# Reads in previously calculated ratios
with open(args.res_file, 'r') as res_open:

    res_csv = csv.reader(res_open, delimiter=',', lineterminator = '\n')
    header = res_csv.next()

    for line in res_csv:

        effect = float(line[args.effect_col])
        gene2effect[line[0]] = effect

for gene, effect in gene2effect.items():

    geneID, name, entrez = retrieveIDs(gene, geneID2Name, geneName2ID, geneEns2Name)

    gene2effect[geneID] = effect
    gene2effect[name] = effect
    gene2effect[entrez] = effect


###############################################################################
#

ref_list = []

# Finds the elements enrichments for each gene, looking by symbol and id
for gene in args.gene_list:

    geneID, name, entrez = retrieveIDs(gene, geneID2Name, geneName2ID, geneEns2Name)

    if gene.upper() in gene_ref:

        elements = gene_ref[gene.upper()]
        effect = gene2effect[gene.upper()]

        ref_list.append((gene, elements, effect))

    elif geneID in gene_ref:

        elements = gene_ref[geneID]
        effect = gene2effect[geneID]

        ref_list.append((gene, elements, effect))

    elif name in gene_ref:

        elements = gene_ref[name]
        effect = gene2effect[name]

        ref_list.append((gene, elements, effect))

    elif entrez in gene_ref:

        elements = gene_ref[entrez]
        effect = gene2effect[entrez]

        ref_list.append((gene, elements, effect))

    else:
        print('Warning: ' + gene + ' not found')


###############################################################################
# 

print('Making graphs')

# A custom color map for easy visualization
cdict = {'red': [(0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)],
        'green': [(0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)],
        'blue': [(0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)]}
custom_map = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 256)

# Determines the appropriate span
span = max(abs(max(neg_rhos)), abs(min(neg_rhos)),
           abs(max(tar_rhos)), abs(min(tar_rhos)))

# Override that appropriate range
if args.xlim:
    span = args.xlim

# Creates an xgrid so kernel distributions can be graphed
xgrid = np.linspace(-span, span, 100)


###############################################################################
#

# Moves through each gene entered, making a graph for each
for gene, ref, effect in ref_list:

    # Information is stored as four lists, each tracking together
    x = []
    y = []
    enrich = []

    # Uses reference to fill in lists
    for rho, element in ref:
        x.append(math.log(untreated[element], 2))
        y.append(math.log(treated[element], 2))
        enrich.append(rho)

    if args.cloud:
    
        # Starts figure as a scatter plot of mildly transparent negative guides
        plt.figure(dpi=400)
        plt.scatter(back_x, back_y, 40, 'lightgrey', alpha=0.25,
                    linewidth=0, label='Negative Control Guides')

        # Places over negative guides a colorful set of guides for the GOI
        plt.scatter(x, y, 40, c=enrich, linewidth=1, cmap=custom_map,
                    vmax=span, vmin=-span, label='Guides targeting ' + gene)

        # Places colorbar, axes labels, title, legend etc
        plt.xlim(min(all_x), max(all_x))
        plt.ylim(min(all_y), max(all_y))

        plt.colorbar()
        plt.xlabel('Log2 counts untreated')
        plt.ylabel('Log2 counts treated')
        plt.title(gene + ' counts')

        if args.legend:
            plt.legend(loc=2)

        plt.axis('image')

        # Saves image by gene name
        plt.savefig(file_out + '_' + gene + '_cloud.' + args.file_type)
        plt.close()

    if args.hist:

        plt.figure(dpi=400)

        hitDist = st.gaussian_kde(enrich)

        # Plots histogram and kernel estimation for negative guides
        plt.plot(xgrid, negDist(xgrid), '-r', linewidth=4)
        plt.plot(xgrid, hitDist(xgrid), '-b', linewidth=4)

        if args.ylim:
            plt.ylim([0, args.ylim])
        plt.ylim(ymin=0)

        ymin, ymax = plt.ylim()
        plt.xlim([-span,span])

        plt.vlines(enrich, 0, ymax / 5.0, colors='blue')
        plt.axvline(x=effect, color='b', ls='dashed', linewidth=2)

        plt.legend(['Null distribution', gene + ' distribution'])
        plt.xlabel('Enrichment')
        plt.ylabel('Frequency')
    
        plt.savefig(file_out + '_' + gene + '_hist.' + args.file_type)
        plt.close()

    print(gene + ' figures created')


###############################################################################
