##############################################################################
#

import argparse
import sys
import csv
import subprocess
import shlex
import os


##############################################################################
# Initiates argument parser
parser = argparse.ArgumentParser(description='Screen type, oligo file, name')

parser.add_argument('screen_type', help='The screen type', type=str)
parser.add_argument('oligo_file', help='Input oligo file', type=str)
parser.add_argument('name', help='Name output files', type=str)

# Optional arguments: base trimming of fasta
parser.add_argument('-s','--strim', help='Number of bases to'
                    'be trimmed from the start', default=0, type=int)

parser.add_argument('-e', '--etrim', help='Number of bases to'
                    'be trimmed from the end', default=0, type=int)

parser.add_argument('-o', '--override', help='Rename indices', action='store_true')


##############################################################################
# Saves input to args object
args = parser.parse_args()

Index_file = os.path.join('Indices','screen_type_index.txt')

# Check whether screen type and name already exist
with open(Index_file, 'r') as index_open:

    Index_csv = csv.reader(index_open, delimiter='\t')

    for rows in Index_csv:

        prev_name, prev_location = rows

        if prev_name == args.screen_type:

            if not args.override:
                sys.exit('Screen name taken')

            else:
                print('Previous screen overwritten')
                #try:
                #    subprocess.check_call('rm ' + prev_location + '.*',
                #                                                    shell=True)
                #except:
                #    sys.exit()


##############################################################################
# Convert oligo file into bowtie-compatible fasta

with open(args.oligo_file, 'r') as oligo_file:
    oligos = csv.reader(oligo_file, delimiter=',')
    oligo_list = []
    
    if args.etrim > 0:
        for row in oligos:
            oligo_list.append(['>' + row[0], row[1][args.strim: -args.etrim]])
    else:
        for row in oligos:
            try:
                oligo_list.append(['>' + row[0], row[1][args.strim: ]])

            except:
                print row
                sys.exit()

fasta_location = os.path.join('Indices', 'temp_fasta.fna')

with open(fasta_location, 'w') as fasta_open:
    fasta_csv = csv.writer(fasta_open, delimiter='\n')
    for i in oligo_list:
        fasta_csv.writerow(i)


##############################################################################
# Call bowtie-build to build new index

try:
    subprocess.check_call('bowtie-build ' + fasta_location + ' ' + 
                                            os.path.join('Indices', args.name),
                                                                    shell=True)
except:
    sys.exit()

# delete fasta files
try:
    subprocess.check_call('rm ' + fasta_location, shell=True)
except:
    sys.exit()


##############################################################################
# Add new index to screentype_index file
Index_name = os.path.join('Indices', args.name)

with open(Index_file,'a') as index_open:
    index_csv = csv.writer(index_open, delimiter='\t')
    index_csv.writerow([args.screen_type, Index_name])


##############################################################################
