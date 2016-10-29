###############################################################################
# David Morgens
# 03/18/2016
###############################################################################
# Imports neccessary modules

'''
Module containing important auxilary screening functions
'''

from __future__ import division
import numpy
import numpy as np
import math
import csv
from collections import defaultdict
import os
import scipy.misc
import scipy.stats as st
import scipy.stats
import sys
import random
import re


###############################################################################
# Function to round off significant digits

def sigDig(x, num=3):
    '''
    Function rounds number in a reasonable manner. Default is to three
    significant digits.
    '''
    # Don't attempt to calculate the log of zero
    if x == 0:
        return 0

    else:
        order_of_magnitude = int(math.floor(math.log10(abs(x))))
        digits = (num - 1) - order_of_magnitude

    return round(x, digits)


#########################################################################
# Calculates the GC content of an inputed string

def getGC(guide):

    numGC = 0
    total = len(guide)

    for nuc in list(guide):
        if nuc == 'C' or nuc == 'G' or nuc == 'c' or nuc == 'g':
            numGC += 1

    return float(numGC) / total


###############################################################################
# Function which retrieves info about genes

def retrieveInfo(ref_base='GenRef', mouse=False):

    '''
    Retrieves gene info for the screen type. Location of reference
    files can be changed, defaults to nearby GenRef folder.
    '''

    # Finds info files downloaded from NCBI
    org_file_human = os.path.join(ref_base, 'Homo_sapiens.gene_info')
    org_file_mouse = os.path.join(ref_base, 'Mus_musculus.gene_info')

    # Custom Ensemble ID to gene name file
    ens_file = os.path.join(ref_base, 'ensRef.csv')

    geneID2Name = defaultdict(lambda: 'N/A')
    geneID2Info = defaultdict(lambda: 'N/A')
    geneName2ID = defaultdict(lambda: 'N/A')
    geneEns2Name = defaultdict(lambda: 'N/A')

    # Reads in Ensemble data
    try:
        with open(ens_file, 'r') as ens_open:

            ens_csv = csv.reader(ens_open, delimiter=',')

            for line in ens_csv:
                geneEns2Name[line[1]] = line[0].upper()

    except IOError:
        print('Ensembl information file not found.\n'
                + 'Use -r to change file location')

    # Reads in Mouse data
    try:
        with open(org_file_mouse, 'r') as org_open:

            org_csv = csv.reader(org_open, delimiter='\t')
            org_csv.next()  # Skips header

            for line in org_csv:
                # Entrez
                geneID2Name[line[1]] = line[2].upper()
                geneID2Info[line[1]] = line[8]
                geneName2ID[line[2].upper()] = line[1]

    except IOError:
        print('Mouse information file not found.\n'
                + 'Use -r to change file location')

    if not mouse:
        # Reads in Human data
        try:
            with open(org_file_human, 'r') as org_open:

                org_csv = csv.reader(org_open, delimiter='\t')
                org_csv.next()  # Skips header

                # For each line in file, save that gene information
                for line in org_csv:

                    geneID2Name[line[1]] = line[2].upper()
                    geneID2Info[line[1]] = line[8]
                    geneName2ID[line[2].upper()] = line[1]


        except IOError:
            print('Human information file not found.\n'
                    + 'Use -r to change file location')

    return geneID2Name, geneID2Info, geneName2ID, geneEns2Name


###############################################################################
# Retreives GO information

def retrieveGO(ref_base='GenRef'):
    '''
    Returns GO component, process, and function data by geneID.
    '''

    go_file = os.path.join(ref_base, 'gene2go')

    # Stores as dictionary of strings. 
    geneID2Comp = defaultdict(str)
    geneID2Proc = defaultdict(str)
    geneID2Fun = defaultdict(str)

    # Checks that file exists, if not returns empty dictionaries
    if not os.path.isfile(go_file):

        print('GO reference file not found; use -r to change file location')
        return geneID2Comp, geneID2Proc, geneID2Fun

    # Reads in GO data
    with open(go_file, 'r') as go_open:
        go_csv = csv.reader(go_open, delimiter='\t')

        for line in go_csv:

            # Checks that line is correct length
            if len(line) == 8:

                # Skips NOT GO terms
                if line[4] == 'NOT':
                    continue

                if line[7] == 'Component':
                    geneID2Comp[line[1]] += line[5] + '|'

                elif line[7] == 'Process':
                    geneID2Proc[line[1]] += line[5] + '|'

                elif line[7] == 'Function':
                    geneID2Fun[line[1]] += line[5] + '|'

    geneID2Comp.default_factory = lambda: 'None'
    geneID2Proc.default_factory = lambda: 'None'
    geneID2Fun.default_factory = lambda: 'None'

    return geneID2Comp, geneID2Proc, geneID2Fun


###############################################################################
# Processes time zero count files

def timeZero(zero_files, thresh):

    # Defaults values to count threshold
    zero_unt = defaultdict(lambda: thresh)
    zero_trt = defaultdict(lambda: thresh)

    # If no time zero file provided, returns defaults only
    if not zero_files:
        return zero_unt, zero_trt
    else:
        zero_unt_file, zero_trt_file = zero_files

    # Reads and filters in time zero untreated file
    with open(zero_unt_file, 'rU') as zero_unt_open:

        dialect = csv.Sniffer().sniff(zero_unt_open.read(1024), delimiters='\t ,')
        zero_unt_csv = csv.reader(zero_unt_open, dialect)
        #zero_unt_csv = csv.reader(zero_unt_open, delimiter=',')

        for line in zero_unt_csv:

            if int(line[1]) > thresh:
                zero_unt[line[0]] = int(line[1])

            else:
                zero_unt[line[0]] = thresh

    # Reads and filters in time zero treated file
    with open(zero_trt_file, 'rU') as zero_trt_open:

        dialect = csv.Sniffer().sniff(zero_trt_open.read(1024), delimiters='\t ,')
        zero_trt_csv = csv.reader(zero_trt_open, dialect)
        #zero_trt_csv = csv.reader(zero_trt_open, delimiter=',')

        for line in zero_trt_csv:

            if int(line[1]) > thresh:
                zero_trt[line[0]] = int(line[1])

            else:
                zero_trt[line[0]] = thresh

    return zero_unt, zero_trt


###############################################################################
# Filters count file by a threshold. If counts are below threshold, redefines
# them to equal threshold.  If counts in both samples are below threshold,
# throws them out.

def filterCounts(unt_file, trt_file, thresh, zero_files, exclude=False):
    '''
    Takes untreated and treated count files and filters them according
    to threshold.
    '''

    # Processes time zero files in auxilary function
    zero_unt_raw, zero_trt_raw = timeZero(zero_files, thresh)

    # Stores untreated counts as dictionary of name to count
    untreated_raw = {}
    treated_raw = {}

    with open(unt_file, 'rU') as unt_open:

        dialect = csv.Sniffer().sniff(unt_open.read(1024), delimiters='\t ,')
        unt_open.seek(0)
        unt_csv = csv.reader(unt_open, dialect)
        #unt_csv = csv.reader(unt_open, delimiter=',')

        for line in unt_csv:

            # Skips blank lines
            if not line or not line[0]:
                continue

            # If no exclusion characters, save line
            if not exclude:
                    untreated_raw[line[0]] = int(float(line[1]))

            # If exclusion character if it does not contain substring
            else:
                for ex in exclude:
                    if ex in line[0]:
                        untreated_raw[line[0]] = int(line[1])
                        break

    # Stores treated counts as dictionary of name to count
    with open(trt_file, 'rU') as trt_open:

        dialect = csv.Sniffer().sniff(trt_open.read(1024), delimiters='\t ,')
        trt_open.seek(0)
        trt_csv = csv.reader(trt_open, dialect)
        #trt_csv = csv.reader(trt_open, delimiter=',')

        for line in trt_csv:

            # Skips blank lines
            if not line or not line[0]:
                continue

            # If no exclusion characters, save line
            if not exclude:
                treated_raw[line[0]] = int(float(line[1]))

            # If exclusion character if it does not contain substring
            else:
                for ex in exclude:
                    if ex in line[0]:
                        treated_raw[line[0]] = int(line[1])
                        break         

    # Tracks some filtering statistics
    belowUnt, belowTrt = 0, 0
    removed = 0

    # Stores filtered counts as dictionary of name to count
    treated = {}
    untreated = {}
    zero_unt = {}
    zero_trt = {}

    # Loops over untreated counts, looks for that entry in the the treated
    # counts. Nonpresence indicates zero counts.  If both counts are less
    # than the threshold, the entry is filtered out.  If one is less, then
    # it is assigned to the threshold value. Elsewise saves the values.
    for entry in untreated_raw:

        # Indicator variable of meeting the threshold in each count file
        un = 0
        tr = 0

        # Checks if over threshold in untreated sample
        if untreated_raw[entry] < thresh:
            un = 1
	    belowUnt += 1

        # Checks if over threshold in treated sample
        if entry not in treated_raw or treated_raw[entry] < thresh:
            tr = 1
            belowTrt += 1

        # If under in both, don't save the entry
        if un and tr:
            removed += 1
            continue

        # If under threshold in untreated, save as threshold value
        if un:
            untreated[entry] = thresh
            
        else:
            untreated[entry] = untreated_raw[entry]

        # If under threshold in treated, save as threshold value
        if tr:
            treated[entry] = thresh
        else:
            treated[entry] = treated_raw[entry]

        # Looks up time zero counts
        zero_unt[entry] = zero_unt_raw[entry]
        zero_trt[entry] = zero_trt_raw[entry]

    # Loops over treated, looking for entries missed in the untreated counts
    for entry in treated_raw:
        if entry not in untreated_raw:

            # If too small in both, do not save.
            if treated_raw[entry] < thresh:
                removed += 1

            # Else save with untreated value equal to threshold
            else:
                treated[entry] = treated_raw[entry]
                untreated[entry] = thresh
                belowUnt += 1

                # Looks up time zero counts
                zero_unt[entry] = zero_unt_raw[entry]
                zero_trt[entry] = zero_trt_raw[entry]

    # Saves stats and time zero files
    stats = (belowTrt, belowUnt, removed)
    time_zero = (zero_unt, zero_trt)

    return untreated, treated, stats, time_zero
 

###############################################################################
# Function to calculate enrichment values

def enrich(count1, sum1, zero1, sum_zero1,
           count2, sum2, zero2, sum_zero2,
           shift, norm):
    '''
    Function calculates enrichment values
    '''

    # Calculates proportions
    prop1 = float(count1) / sum1
    prop2 = float(count2) / sum2
    prop_zero1 = float(zero1) / sum_zero1
    prop_zero2 = float(zero2) / sum_zero2

    # Calculates ratio and log ratio 
    log_enrich = math.log(prop1 / prop2, 2) - math.log(prop_zero1 / prop_zero2, 2)

    # Normalizes log ratio by shifting around 'zero' and stretching
    # appropriately
    shift_enrich = log_enrich - shift
    norm_enrich = shift_enrich / norm

    return norm_enrich


###############################################################################
# Function to calculate enrichments

def enrich_all(untreated, treated, neg_name, split_mark, K, time_zero, back):
    '''
    Auxilary function to calculate enrichment values
    '''

    # Finds total counts in time zero files
    zero_unt, zero_trt = time_zero
    total_zero_unt = sum(zero_unt.values())
    total_zero_trt = sum(zero_trt.values())

    # Finds total counts in count files
    total_unt = sum(untreated.values())
    total_trt = sum(treated.values())

    # Stores enrichments of negative controls
    neg_raw = []
    tar_raw = []
    back_raw = []

    # Moves over each element, and, if it is a control element, calculates its enrichment
    for entry in untreated:

        # Select only negative controls
        if entry.split(split_mark)[0].startswith(neg_name):

            # Calls enrichment function with a 0 shift
            neg_raw.append(enrich(treated[entry], total_trt,
                                    zero_trt[entry], total_zero_trt,
                                    untreated[entry], total_unt,
                                    zero_unt[entry], total_zero_unt,
                                    0, 1))
        else:
            tar_raw.append(enrich(treated[entry], total_trt,
                                    zero_trt[entry], total_zero_trt,
                                    untreated[entry], total_unt,
                                    zero_unt[entry], total_zero_unt,
                                    0, 1))

    if back == 'neg':
        shift = np.median(neg_raw)  # Calculates the shift as a median

    elif back == 'tar':
        shift = np.median(tar_raw)

    elif back == 'all':
        shift = np.median(neg_raw + tar_raw)

    else:
        sys.exit('Unrecognized option for background choice: ' + back)

    entry_rhos = {}

    # With the shift in hand, calculates the enrichment for each guide
    for entry in treated:

        # Note the calculated neg_shift is used now
        entry_rhos[entry] = enrich(treated[entry], total_trt,
                                    zero_trt[entry], total_zero_trt,
                                    untreated[entry], total_unt,
                                    zero_unt[entry], total_zero_unt,
                                    shift, K)

    # Gatheros rho values for each gene and seperates out negative controls
    gene_rhos = defaultdict(list)
    gene_rhos_int = defaultdict(list)
    gene_ref = defaultdict(list)

    # Stores all negative element rhos and targeting rhos
    neg_rhos = []
    tar_rhos = []

    for entry in entry_rhos:

        # Checks if entry is a negative control
        if entry.split(split_mark)[0].startswith(neg_name):
            neg_rhos.append(entry_rhos[entry])

        else:

            # Gathers rhos of elements targeting each gene
            gene = entry.split(split_mark)[0].upper()
            gene_rhos[gene] += [entry_rhos[entry]]

            # Saves element name and enrichment for output
            entry_split = entry.split(split_mark)
            gene_ref[gene] += [(entry_rhos[entry], entry)]
            tar_rhos.append(entry_rhos[entry])

    return entry_rhos, gene_rhos, neg_rhos, tar_rhos, gene_ref


###############################################################################
# The likelihood function for casTLE

def likeGrid(rhos, I, hit_rate, back_dist):

    '''
    Takes precomputed likelihoods and free parameter values and returns the
    log likelihood.
    '''

    like = 0  # Initiates log likelihood

    # Calculates negative rate by elimination
    back_rate = 1 - hit_rate

    # Checks that rates make sense
    if back_rate < -0.01 or back_rate > 1.01:
        sys.exit('Error: Impossible background rate')

    # Normalization constant
    hit_norm = back_dist[0] * abs(I) + 1

    # For each entry, determines whether it falls within the hit region, then
    # calculate the appropriate likelihood
    if I < 0:
        for rho in rhos:

            # Indicates enrichment is below the effect estimate, meaning most likely
            # true effect is I
            if rho < I:
                like += math.log(hit_rate * back_dist[rho - I] / hit_norm +
                                        back_rate * back_dist[rho])

            # Indicates enrichment is the opposite sign, meaning most likely
            # true effect is 0
            elif rho > 0:
                like += math.log(hit_rate * back_dist[rho] / hit_norm +
                                        back_rate * back_dist[rho])

            # Indicates enrichment is within the bounded region, meaning most likely
            # true effect is rho
            else:
                like += math.log(hit_rate * back_dist[0] / hit_norm +
                                        back_rate * back_dist[rho])

    elif I > 0:
        for rho in rhos:

            # Indicates enrichment is above the effect estimate, meaning most likely
            # true effect is I
            if rho > I:
                like += math.log(hit_rate * back_dist[rho - I] / hit_norm +
                                        back_rate * back_dist[rho])

            # Indicates enrichment is the opposite sign, meaning most likely
            # true effect is 0
            elif rho < 0:
                like += math.log(hit_rate * back_dist[rho] / hit_norm +
                                        back_rate * back_dist[rho])

            # Indicates enrichment is within the bounded region, meaning most likely
            # true effect is rho
            else:
                like += math.log(hit_rate * back_dist[0] / hit_norm +
                                        back_rate * back_dist[rho])

    else:
        for rho in rhos:

            # If I is zero, then true effect is most likely 0
            like += math.log(back_rate * back_dist[rho])

    return like


###############################################################################
# Finds credible intervals

def findInterval(log_data, target, step):

    # Calculates likelihood fractions for interval calculation
    norm_log_data = log_data - scipy.misc.logsumexp(log_data[~numpy.isnan(log_data)])
    data = numpy.exp(norm_log_data)

    # Sanity check
    if numpy.nansum(data) <= 0.95 or numpy.nansum(data) >= 1.05:
        print(str(numpy.nansum(data)))
        sys.exit('Error: Unstable computation')

    # Initiates search at the starting value
    start = numpy.nanargmax(data)
    total_weight = data[start] + 0
    min_ind, max_ind = start, start

    steps = 0

    left_edge, right_edge = False, False

    # Continues adding to the interval until the target weight is reached
    while 1:

        steps += 1

        # Ends if total weight exceeded
        if total_weight >= target:
            break

	# Checks it has not reached the edge
        if min_ind == 0:
            left_edge = True
        else:
            left = data[min_ind - step]

            if numpy.isnan(left):
                left_edge = True

        if max_ind >= (len(data) - step):
            right_edge = True
        else:
            right = data[max_ind + step]

            if numpy.isnan(right):
                right_edge = True

	# If it gets stuck, error out
        if right_edge and left_edge:
            print(str(steps))
            sys.exit('Error: Interval calculation failure')

        elif right_edge:
            total_weight += left
            min_ind -= step

        elif left_edge:
            total_weight += right
            max_ind += step

	# Picks the larger of the neighboring likelihoods
        elif left >= right:
            total_weight += left
            min_ind -= step

        else:
            total_weight += right
            max_ind += step

    return total_weight, min_ind, max_ind


###############################################################################
# The subprocess code, which runs the parallel processes

def trial(gene_rhos, back_dist, gene_span, step):
    '''
    Function that runs subprocesses
    '''

    # Output stored in dictionaries
    geneI = {}
    geneL = {}
    geneInterval = {}
    geneDist = {}

    # For each gene given, find the MLE
    for gene in gene_rhos:

        # Retrieves the rho values and precomputes their likelihood to be drawn
        # from the negative distributions, note lower bound to prevent underflows
        rhos = gene_rhos[gene]

        # Sets up the grid search,
        min_I, max_I = gene_span[gene]
        pos_hit_rate = numpy.linspace(0.1, 0.9, 9)

        # Finds the likelihood of the null model, where no elements work slash
        # the gene has no effect
        I = 0
        hit_rate = 0

        like0 = likeGrid(rhos, I, hit_rate, back_dist)

        # Initiates grid of likelihoods
        all_dist = {}
        all_dist[0] = [like0]
        I2marglog = {}

        # Set up grid calculation
        buff_zero = 10 ** 3

        marg_logs = numpy.zeros(max_I - min_I + buff_zero * 2) - 1000
        marg_logs[:] = numpy.NAN
        marg_logs[-1 * min_I + buff_zero] = like0

	# Grid search across the positive values
        for I in range(step, max_I + step, step):

            # Determines effect size
            dist = []

            # For each fraction of effective reagents, calculate and save likelihood
            for hit_rate in pos_hit_rate:

                # Finds likelihood of parameter combination
                like = likeGrid(rhos, I, hit_rate, back_dist)
                dist.append(like)

            # Saves raw likelihoods for additional analysis
            all_dist[I] = dist

	    # Performs marginalization in log scale for stability
            marg_log = scipy.misc.logsumexp(dist) - math.log(len(dist))
            marg_logs[I - min_I + buff_zero] = marg_log

	# Grid search across the negative values
        I = 0

        for I in range(-1 * step, min_I - step, -1 * step):

            # Determines effect size
            dist = []

            # For each fraction of effective reagents, calculate and save likelihood
            for hit_rate in pos_hit_rate:

                # Finds likelihood of parameter combination
                like = likeGrid(rhos, I, hit_rate, back_dist)
                dist.append(like)

            # Saves raw likelihoods for additional analysis
            all_dist[I] = dist

	    # Performs marginalization in log scale for stability
            marg_log = scipy.misc.logsumexp(dist) - math.log(len(dist))
            marg_logs[I - min_I + buff_zero] = marg_log

        total_weight, min_ind, max_ind = findInterval(marg_logs, 0.95, step)

        # Finds maximum likelihood
        maxL = numpy.nanmax(marg_logs)
        maxI = numpy.nanargmax(marg_logs) + min_I - buff_zero

        # Saves most likely parameters, along with likelihood ratio
        geneDist[gene] = all_dist
        geneI[gene] = maxI

        # Calculates log-likelihood ratio
	geneL[gene] = 2 * (maxL - like0)
        geneInterval[gene] = (min_ind + min_I - buff_zero,
                                max_ind + min_I - buff_zero, total_weight)

        # If enrichments are not provided, returns 0s for downstream analysis
        if not rhos:
            geneI[gene] = 0
            geneL[gene] = 0
            geneInterval[gene] = (0, 0, 1)

    return geneI, geneL, geneInterval, geneDist


###############################################################################
# Function to calculate likelihoods

def retrieveLikelihoods(gene_rhos, back_rhos, nums, gene_span, scale, I_step):

    #
    gene_rhos_int, gene_span_int = intefy(gene_rhos, gene_span, scale)
    step_scaled = int(round(I_step * (10 ** scale)))
    back_rhos_calc = precalculateGrid(gene_span_int, back_rhos, scale)

    # Checks if single processor
    if nums == 1:
        single = True
    else:
        single = False

    # Checks and creates parallel environment
    try:
        if not single:
            import pp

    except ImportError:
        print('Parallel Python package pp not found. Defaulting to single core')
        single = True

    if not single:

        # Initiates parallel server
        job_server = pp.Server(nums)
        fn = pp.Template(job_server, trial, (likeGrid, findInterval,),
                ("numpy", "math", "scipy.misc"))

    # Single core version
    if single:
        geneI, geneL, geneInterval, geneDist = trial(gene_rhos_int, back_rhos_calc, gene_span_int, step_scaled)

    # Parallel version
    if not single:

        geneI = {}
        geneL = {}
        geneInterval = {}
        geneDist = {}

        GeneRhosSplit = []
        keys = gene_rhos_int.keys()

	# Determines the number of genes in each task
        n = int(len(keys) / nums) + 1
        keysSplit = [keys[i: i + n] for i in range(0, len(keys), n)]

        # Splits dictionary
        for keyList in keysSplit:
            dictChunk = {}

            for key in keyList:
                dictChunk[key] = gene_rhos_int[key]

            GeneRhosSplit.append(dictChunk)

	# Submits individual jobs for each gene list
        jobs = []
        for splits in GeneRhosSplit:
            jobs.append(fn.submit(splits, back_rhos_calc, gene_span_int, step_scaled))

	# Retrieves results from individual jobs
        for job in jobs:
            val = job()

            try:
                geneIi, geneLi, geneIntervali, geneDisti = val

            # Catches subprocess errors
            except TypeError:
                sys.exit('Subprocess failed')

            # Saves results
            geneI.update(geneIi)
            geneL.update(geneLi)
            geneInterval.update(geneIntervali)
            geneDist.update(geneDisti)

    return geneI, geneL, geneInterval, geneDist


###############################################################################
# Calculates permutations for casTLE p-values

def retrievePerm(draw_num, perm_num, back_rhos, tar_rhos, nums, scale, I_step):

    # Generates random 'genes' from all targeting elements
    perm_rhos = dict([(i, random.sample(tar_rhos, draw_num)) for i in range(perm_num)])

    # Calculates grid search
    gene_span = {}

    for gene, rhos in perm_rhos.items():

        max_I = int(1.2 * max(rhos + [-1, 0, 1]))
        min_I = int(1.2 * min(rhos + [-1, 0, 1]))

        gene_span[gene] = (min_I, max_I)

    # Finds the likelihood of these 'genes'
    permI, permL, permInterval, permDist = retrieveLikelihoods(perm_rhos, back_rhos,
                                        nums, gene_span, scale, I_step)

    return permI, permL, permInterval


###############################################################################
# Short script for calculating p-values based on empirical distribution

def rankLikelihoods(perm_rats, gene2rat):

    # Retrieve actual gene log-likelihood ratios
    genes_rats = gene2rat.items()
    genes = [x[0] for x in genes_rats]
    rats = [x[1] for x in genes_rats]

    # Ranks actual genes based on permutations
    rat_rank = numpy.searchsorted(sorted(perm_rats), rats, 'right')
    gene_rank = dict(zip(genes, rat_rank))

    geneP = {}

    # Converts ranks to p values
    for gene in gene_rank:
        geneP[gene] = 1 - (gene_rank[gene] - 1) / float(len(perm_rats))

    return geneP


###############################################################################
# Combines data from two screens

def retrieveCombo(data1, data2, gene_span, I_step):
    '''
    Function that calculates combo scores
    '''

    # Unpack data for both screens
    geneI1, geneL1, geneInterval1, geneDist1 = data1
    geneI2, geneL2, geneInterval2, geneDist2 = data2

    geneI = {}
    geneL = {}
    geneInterval = {}

    # Unify results
    for gene in gene_span:

        # If gene only in one screen, use those results
        if gene not in geneDist2:
            geneI[gene] = geneI1[gene]
            geneL[gene] = geneL1[gene]
            geneInterval[gene] = geneInterval1[gene]
            continue

        if gene not in geneDist1:
            geneI[gene] = geneI2[gene]
            geneL[gene] = geneL2[gene]
            geneInterval[gene] = geneInterval2[gene]
            continue

        # Retrieve likelihood grid for each
        dist1 = geneDist1[gene]
        dist2 = geneDist2[gene]

        # Calculate new null likelihood
        like01 = dist1[0][0]
        like02 = dist2[0][0]
        like0 = like01 + like02

        pos_I = sorted(dist1.keys())
        marg_logs = []

        for I in pos_I:

            if I == 0:
                marg_logs.append(like0)
                continue

            dist1I = [like01] + dist1[I]
            dist2I = [like02] + dist2[I]

            skip0 = True

            dist = []

            for like1 in dist1I:
                for like2 in dist2I:

                    if skip0:
                        skip0 = False
                        continue

                    dist.append(like1 + like2)

            marg_log = scipy.misc.logsumexp(dist) - math.log(len(dist))
            marg_logs.append(marg_log)

	# Retrieves the 95% confidence interval
        total_weight, min_ind, max_ind = findInterval(numpy.asarray(marg_logs), 0.95, 1)

        # Finds maximum likelihood
        maxL = max(marg_logs)

        # Saves most likely parameters, along with likelihood ratio
        geneI[gene] = pos_I[marg_logs.index(maxL)]
        geneL[gene] = 2 * (maxL - like0)
        geneInterval[gene] = (pos_I[min_ind], pos_I[max_ind], total_weight)

    return geneI, geneL, geneInterval


###############################################################################
# Defines span for gene combo while unifying gene IDs

def comboSpan(gene_rhos1, gene_rhos2, span):

    # Retrieves ID maps
    geneID2Name, geneID2Info, geneName2ID, geneEns2Name = retrieveInfo()

    gene_span = {}
    add_gene_rhos1 = {}
    add_gene_rhos2 = {}

    for gene, rhos1 in gene_rhos1.items():

        # Converts IDs
        geneID, name, entrez = retrieveIDs(gene, geneID2Name, geneName2ID, geneEns2Name)

        # Searches for ID in other screen
        if gene in gene_rhos2:
            unified = gene
            rhos2 = gene_rhos2[gene]

        elif geneID in gene_rhos2:
            unified = geneID
            rhos2 = gene_rhos2[geneID]

        elif name in gene_rhos2:
            unified = name
            rhos2 = gene_rhos2[name]

        elif entrez in gene_rhos2:
            unified = entrez
            rhos2 = gene_rhos2[entrez]

        else:
            unified = geneID
            rhos2 = []

        # Initiates grid
        max_I = int(span * max(rhos1 + rhos2 + [-1, 0, 1]))
        min_I = int(span * min(rhos1 + rhos2 + [-1, 0, 1]))

        # Saves grid by unified ID
        gene_span[unified] = (min_I, max_I)
        add_gene_rhos1[unified] = rhos1
        add_gene_rhos2[unified] = rhos2
        
    for gene, rhos2 in gene_rhos2.items():

        # Converts IDs
        geneID, name, entrez = retrieveIDs(gene, geneID2Name, geneName2ID, geneEns2Name)

        # Checks if gene already added
        if gene in add_gene_rhos1:
            continue

        elif geneID in add_gene_rhos1:
            continue

        elif name in add_gene_rhos1:
            continue

        elif entrez in add_gene_rhos1:
            continue

        else:
            rhos1 = []

        # Initiates grid
        max_I = int(span * max(rhos1 + rhos2 + [-1, 0, 1]))
        min_I = int(span * min(rhos1 + rhos2 + [-1, 0, 1]))

        # Saves grid by unified ID
        gene_span[geneID] = (min_I, max_I)
        add_gene_rhos1[geneID] = rhos1
        add_gene_rhos2[geneID] = rhos2
        
    return add_gene_rhos1, add_gene_rhos2, gene_span


###############################################################################
# Function to calculate permutations for casTLE pvalues on combinations

def comboPerm(draw_num1, draw_num2, perm_num, back_rhos1, back_rhos2,
                                tar_rhos1, tar_rhos2, nums, I_step, scale):

    # Creates 'genes'
    perm_rhos1 = dict([(i, random.sample(tar_rhos1, draw_num1)) for i in range(perm_num)])
    perm_rhos2 = dict([(i, random.sample(tar_rhos2, draw_num2)) for i in range(perm_num)])

    # Calculates grids in seperate function
    add_perm_rhos1, add_perm_rhos2, perm_span = comboSpan(perm_rhos1, perm_rhos2, 1.2)

    # Calculates likelihoods seperately
    print('Permuting first screen')
    data1 = retrieveLikelihoods(add_perm_rhos1, back_rhos1, nums,
                                        perm_span, scale, I_step)

    print('Permuting second screen')
    data2 = retrieveLikelihoods(add_perm_rhos2, back_rhos2, nums,
                                        perm_span, scale, I_step)

    # Combines likelihoods
    print('Combining permutations')
    permI, permL, permInterval = retrieveCombo(data1, data2, perm_span, I_step)

    return permI, permL, permInterval


###############################################################################
# Converts dictionary to integers

def intefy(gene_rhos, gene_span, scale):

    magnitude = int(10 ** scale)

    gene_rhos_int = {}
    gene_span_int = {}

    for gene, rhos in gene_rhos.items():

        rhos_int = []

        for rho in rhos:

            rho_int = int(round(magnitude * rho))
            rhos_int.append(rho_int)

        gene_rhos_int[gene] = rhos_int

    for gene, span in gene_span.items():
        min_I, max_I = span

        span_int = (int(round(magnitude * min_I)), int(round(magnitude * max_I)))

        gene_span_int[gene] = span_int

    return gene_rhos_int, gene_span_int


###############################################################################
# Precalculates background distribution on a grid

def precalculateGrid(gene_span, back_rhos, scale):

    magnitude = int(10 ** scale)
    mines, maxes = [], []

    for gene in gene_span:

        min_I, max_I = gene_span[gene]

        mines.append(min_I)
        maxes.append(max_I)

    max_all = 2 * max(maxes)
    min_all = 2 * min(mines)

    grid = range(min_all, max_all)

    back_dist_scaled = [rho * magnitude for rho in back_rhos]
    back_dist = scipy.stats.gaussian_kde(back_dist_scaled)

    back_rhos_calc = {}
    # Parameter to prevent underflow, may effect outlier tolerance
    underflow = 10 ** -10

    for point in grid:

        back = back_dist(point)
        back_rhos_calc[point] = back + underflow

    return back_rhos_calc


###############################################################################
# Script to process result records

def retrieveRecord(res_file, current_version):

    name = res_file[: -4]
    rec_file = name + '_record.txt'

    try:
        # Parses record file
        with open(rec_file, 'r') as rec_open:
            rec_csv = csv.reader(rec_open, delimiter='\t')
            script, version = rec_csv.next()

            if version != current_version:
                sys.exit('Error: Version number not current\n'
                            + 'Rerun analysis')

            if script != 'analyzeCounts.py':
                sys.exit('Error: Input is not a result file')

            last_time = rec_csv.next()[1]
            unt_file = rec_csv.next()[1]
            trt_file = rec_csv.next()[1]
            zero_files = rec_csv.next()[1]
            if zero_files:
                zero_files = eval(zero_files)
            file_out = rec_csv.next()[1]
            screen_type = rec_csv.next()[1]
            neg_name = rec_csv.next()[1]
            split_mark = rec_csv.next()[1]
            exclude = rec_csv.next()[1]
            if exclude:
                exclude = eval(exclude)
            thresh = int(rec_csv.next()[1])
            K = float(rec_csv.next()[1])
            back = rec_csv.next()[1]
            I_step = float(rec_csv.next()[1])
            scale = int(rec_csv.next()[1])
            draw_num = int(rec_csv.next()[1])

    except IOError:
        sys.exit('Error: Record of result file not found\n'
                    + 'Change file name or rerun analysis')

    # Saves parameters for processing
    stats = (script, version, last_time)
    files = (unt_file, trt_file, zero_files, file_out)
    info = (screen_type, neg_name, split_mark, exclude)
    param = (thresh, K, back, I_step, scale, draw_num)

    return stats, files, info, param


###############################################################################
# Script to save and use permutations

def calculatePval(res_file, permI, permL, erase, ratio_col):

    name = res_file[: -4]
    ref_file = name + '_ref.csv'

    # Unless prompted, remember previous permutations
    if erase:

        try:
            os.remove(ref_file)
            print('Warning: Previous reference removed')

        except OSError:

            print('Warning: No previous reference found')

    # Write new permutations to file
    with open(ref_file, 'a') as ref_open:
        ref_csv = csv.writer(ref_open, delimiter=',', lineterminator='\n')

        for number in permI:
            ref_csv.writerow([permI[number], permL[number]])

    perm_rats = []

    # Read back in both old and new permutations
    with open(ref_file, 'r') as ref_open:
        ref_csv = csv.reader(ref_open, delimiter=',', lineterminator='\n')

        for line in ref_csv:
            perm_rats.append(float(line[1]))

    gene2line = {}
    gene2rat = {}

    # Reads in previously calculated ratios
    with open(res_file, 'r') as res_open:
        res_csv = csv.reader(res_open, delimiter=',', lineterminator = '\n')
        header = res_csv.next()
        for line in res_csv:
            gene2line[line[0] + line[1]] = line
            gene2rat[line[0] + line[1]] = float(line[ratio_col])

    geneP = rankLikelihoods(perm_rats, gene2rat)

    all_perm_num = str(len(perm_rats))

    return gene2line, geneP, header, all_perm_num


###############################################################################
# Converts IDs

def retrieveIDs(gene, geneID2Name, geneName2ID, geneEns2Name):

    if gene in geneID2Name:
	geneID = gene
        entrez = gene
	name = geneID2Name[geneID]

    elif gene in geneName2ID:
	name = gene
	geneID = geneName2ID[gene]
	entrez = geneName2ID[gene]

    elif gene in geneEns2Name:
        name = geneEns2Name[gene]
        geneID = gene
        entrez = geneName2ID[name]

    else:
	geneID = gene
	name = gene
        entrez = gene

    return geneID, name, entrez


###############################################################################
