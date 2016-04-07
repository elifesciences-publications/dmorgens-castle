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
parser = argparse.ArgumentParser(description='Visualizes replicate data element-wise')

# Non-optional arguments: The files containing results, as well as an output
parser.add_argument('res_file1', help='File for untreated results', type=str)

parser.add_argument('res_file2', help='File for treated results', type=str)

parser.add_argument('name', help='Name for output file', type=str)

# Optional arguments:
parser.add_argument('-of', '--override_file', action='store_true',
                        help='Override automatic targeting to Results folder')

parser.add_argument('-x', '--x_axis', default='replicate 1',
                        help='Label for x axis. Default is "replicate 1"')

parser.add_argument('-y', '--y_axis', default='replicate 2',
                        help='Label for y axis. Default is "replicate 2"')

parser.add_argument('-f', '--file_type', default='png',
                        help='File ending/type. Default is "png"')

# Arguments for statistics
parser.add_argument('-l', '--line', action='store_false',
                        help='Don\'t include linear regression line')

parser.add_argument('-t', '--title', action='store_false',
                        help='Don\'t include title')

# Saves all input to object args
args = parser.parse_args()

bins = 50

###############################################################################
# Checks arguments

if args.override_file:
    file_out = args.name
else:
    file_out = os.path.join('Results', args.name)

try:
    with open(file_out + '_rho.' + args.file_type, 'w') as out_open:
        pass

    os.remove(file_out + '_rho.' + args.file_type)

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

if screen_type1 != screen_type2:
    print('Warning: Screen_types do not match.')
    args.rhos = False


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
# Plots guides against guides and calculates linear regression

print('Plotting figures')

x = []
y = []

# Collects all entries
for element in element_rhos1:
    if element in element_rhos2:
        x.append(element_rhos1[element])
        y.append(element_rhos2[element])

# Finds linear fit
slope, intercept, r_value, p_value, std_err = st.linregress(x,y)

plt.figure()
plt.hist2d(x, y, (100, 100), cmap=cm.jet, norm=LogNorm())
plt.xlabel('Enrichment for ' + args.x_axis)
plt.ylabel('Enrichment for ' + args.y_axis)

# Plots fit line
if args.line:
    xgrid = np.linspace(min(x), max(x), 100)
    plt.plot(xgrid, slope * xgrid + intercept,'-', color='0')

# Adds title
if args.title:
    plt.title('Reproducibility of enrichment by element; R2 = ' + str(round(r_value ** 2, 2)))

plt.savefig(file_out + '_rho.' + args.file_type)
plt.close()

print('Finished')


###############################################################################
