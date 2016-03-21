###############################################################################
# David Morgens
# 11/01/15
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
        digits = (num-1) - order_of_magnitude

    return round(x, digits)


###############################################################################
# Function which retrieves info about genes

def retrieveInfo(ref_base='GenRef'):

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
    with open(zero_unt_file, 'r') as zero_unt_open:

        zero_unt_csv = csv.reader(zero_unt_open, delimiter='\t')

        for line in zero_unt_csv:

            if int(line[1]) > thresh:
                zero_unt[line[0]] = int(line[1])

            else:
                zero_unt[line[0]] = thresh

    # Reads and filters in time zero treated file
    with open(zero_trt_file, 'r') as zero_trt_open:

        zero_trt_csv = csv.reader(zero_trt_open, delimiter='\t')

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

    with open(unt_file, 'r') as unt_open:
        unt_csv = csv.reader(unt_open, delimiter='\t')

        for line in unt_csv:

            # Skips blank lines
            if not line:
                continue

            # If no exclusion characters, save line
            if not exclude:
                untreated_raw[line[0]] = int(line[1])

            # If exclusion character if it does not contain substring
            else:
                if exclude not in line[0]:
                    untreated_raw[line[0]] = int(line[1])  

    # Stores treated counts as dictionary of name to count
    with open(trt_file, 'r') as trt_open:
        trt_csv = csv.reader(trt_open, delimiter='\t')

        for line in trt_csv:

            # Skips blank lines
            if not line:
                continue

            # If no exclusion characters, save line
            if not exclude:
                treated_raw[line[0]] = int(line[1])

            # If exclusion character if it does not contain substring
            else:
                if exclude not in line[0]:
                    treated_raw[line[0]] = int(line[1])              

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

def enrich_all(untreated, treated, neg_name, split_mark, K, time_zero):
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

    neg_shift = np.median(neg_raw)  # Calculates the shift as a median
    entry_rhos = {}

    # With the shift in hand, calculates the enrichment for each guide
    for entry in treated:

        # Note the calculated neg_shift is used now
        entry_rhos[entry] = enrich(treated[entry], total_trt,
                                    zero_trt[entry], total_zero_trt,
                                    untreated[entry], total_unt,
                                    zero_unt[entry], total_zero_unt,
                                    neg_shift, K)

    # Gathers rho values for each gene and seperates out negative controls
    gene_rhos = defaultdict(list)
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
            gene_ref[gene] += [(entry_rhos[entry], entry)]
            tar_rhos.append(entry_rhos[entry])

    return entry_rhos, gene_rhos, neg_rhos, tar_rhos, gene_ref


###############################################################################
# The likelihood function for casTLE

def likeEB(rhos, I, hit_rate, hit_like, back_likes, back_dist, off_likes, off_rate):

    '''
    Takes precomputed likelihoods and free parameter values and returns the
    log likelihood.
    '''

    like = 0  # Initiates log likelihood

    # Calculates negative rate by elimination
    back_rate = 1 - hit_rate - off_rate

    # Checks that rates make sense
    if back_rate < -0.01 or back_rate > 1.01:
        sys.exit('Error: Impossible off target rate')

    # Normalization constant
    hit_norm = hit_like * abs(I) + 1

    # For each entry, determines whether it falls within the hit region, then
    # calculate the appropriate likelihood
    if I < 0:
        for back_like, off_like, rho in zip(back_likes, off_likes, rhos):

            # Indicates enrichment is below the effect estimate, meaning most likely
            # true effect is I
            if rho < I:
                like += math.log(hit_rate * back_dist(rho - I) / hit_norm +
                                        back_rate * back_like + off_rate * off_like)

            # Indicates enrichment is the opposite sign, meaning most likely
            # true effect is 0
            elif rho > 0:
                like += math.log(hit_rate * back_like / hit_norm +
                                        back_rate * back_like + off_rate * off_like)

            # Indicates enrichment is within the bounded region, meaning most likely
            # true effect is rho
            else:
                like += math.log(hit_rate * hit_like / hit_norm +
                                        back_rate * back_like + off_rate * off_like)

    elif I > 0:
        for back_like, off_like, rho in zip(back_likes, off_likes, rhos):

            # Indicates enrichment is above the effect estimate, meaning most likely
            # true effect is I
            if rho > I:
                like += math.log(hit_rate * back_dist(rho - I) / hit_norm +
                                        back_rate * back_like + off_rate * off_like)

            # Indicates enrichment is the opposite sign, meaning most likely
            # true effect is 0
            elif rho < 0:
                like += math.log(hit_rate * back_like / hit_norm +
                                        back_rate * back_like  + off_rate * off_like)

            # Indicates enrichment is within the bounded region, meaning most likely
            # true effect is rho
            else:
                like += math.log(hit_rate * hit_like / hit_norm +
                                        back_rate * back_like  + off_rate * off_like)

    else:
        for back_like, off_like, rho in zip(back_likes, off_likes, rhos):

            # If I is zero, then true effect is most likely 0
            like += math.log(back_rate * back_like + off_rate * off_like)

    return like


###############################################################################
# Alternate likelihood function with no empirical baysian framework

def likeshRNA(rhos, I, hit_rate, hit_like, back_likes, back_dist, off_likes, off_rate):
    '''
    Takes precomputed likelihoods and free parameter values and returns the
    log likelihood.
    '''

    like = 0  # Initiates log likelihood    

    # Calculates negative rate by elimination
    back_rate = 1 - on_rate - back_like

    if back_rate < -0.01 or back_rate > 1.01:
        sys.exit('Impossible off target rate')

    # For each entry, determines whether it falls within the hit region, then
    # calculate the appropriate likelihood
    if I < 0:
        for back_like, off_like, rho in zip(back_likes, off_likes, rhos):

            if rho <= 0 and rho >= I:
                like += math.log(hit_rate * 1.0 / abs(I) +
                                        back_rate * back_like + off_rate * off_like)

            else:
                like += math.log(back_rate * back_like + off_rate * off_like)

    elif I > 0:
        for back_like, off_like, rho in zip(back_likes, off_likes, rhos):

            if rho >= 0 and rho <= I:
                like += math.log(hit_rate * 1.0 / abs(I) +
                                        back_rate * back_like + off_rate * off_like)

            else:
                like += math.log(back_rate * back_like + off_rate * off_like)

    else:
        for back_like, off_like, rho in zip(back_likes, off_likes, rhos):
            like += math.log(back_rate * back_like + off_rate * off_like)

    return like


###############################################################################
# Alternative likelihood function using a shifted noise distribution

def casLike(rhos, I, hit_rate, hit_like, back_likes, back_dist, off_likes, off_rate):
    '''
    Takes precomputed likelihoods and free parameter values and returns the
    log likelihood.
    '''

    like = 0  # Initiates log likelihood

    # Calculates negative rate by elimination
    back_rate = 1 - hit_rate - off_rate

    # Checks that rates make sense
    if back_rate < -0.01 or back_rate > 1.01:
        sys.exit('Impossible off target rate')

    # Calculates the likelihood of the hit distribution
    for back_like, off_like, rho in zip(back_likes, off_likes, rhos):
        
        like += math.log(hit_rate * back_dist(rho - I) +
                            back_rate * back_like + off_rate * off_like)

    return like 


###############################################################################
# Finds confidence/credible intervals

def findInterval(data, target, start):

    # Normalizes as fraction of weight
    total = float(sum(data))
    new_data = []

    for datum in data:

        if datum < 0:
            sys.exit('Negative likelihood error')

        new_data.append(datum / total)

    # Initiates search at the starting value
    total_weight = new_data[start] + 0
    min_ind = start
    max_ind = start

    # Continues adding to the interval until the target weight is reached
    while 1:

        # Ends if total weight exceeded
        if total_weight >= target:
            break

	# Checks it has not reached the edge
        if not min_ind == 0:
            left = data[min_ind - 1]
        else:
            left = 0

        if not max_ind == (len(data) - 1):
            right = data[max_ind + 1]
        else:
            right = 0

	# Picks the larger of the neighboring likelihoods
        if left >= right:
            total_weight += left
            min_ind -= 1

        else:
            total_weight += right
            max_ind += 1

	# If it gets stuck, error out
        if right == 0 and left == 0:
            sys.exit('Interval calculation failure')

    return total_weight, min_ind, max_ind


###############################################################################
# The subprocess code, which runs the parallel processes

def trial(gene_rhos, back_dist, off_dist, off_rate, like_fun, gene_span, I_step):
    '''
    Function that runs subprocesses
    '''

    # Precalculate likelihood of observing true effect
    hit_like = back_dist(0)

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
        back_likes = back_dist(rhos) + 0.0000001
        off_likes = off_dist(rhos) + 0.0000001

        # Sets up the grid search,
        min_I, max_I = gene_span[gene]
        pos_hit_rate = numpy.linspace(0.1, 0.9, 9)

        # Finds the likelihood of the null model, where no elements work slash
        # the gene has no effect
        I = 0
        hit_rate = 0

        like0 = like_fun(rhos, I, hit_rate, hit_like, back_likes,
                                back_dist, off_likes, off_rate)

        # Initiates grid of likelihoods
        pos_I = [0]
        marg_dist_log = [like0]
        all_dist = {}
        all_dist[0] = [like0]

	# Grid search across the positive values
        while I < max_I:

            # Determines effect size
            I += I_step
            pos_I.append(I)
            dist = []

            # For each fraction of effective reagents, calculate and save likelihood
            for hit_rate in pos_hit_rate:

                # Finds likelihood of parameter combination
                like = like_fun(rhos, I, hit_rate, hit_like, back_likes,
                                                back_dist, off_likes, off_rate)
                dist.append(like)

            # Saves raw likelihoods for additional analysis
            all_dist[round(I, 2)] = dist

	    # Performs marginalization in log scale for stability
            marg_log = scipy.misc.logsumexp(dist) - math.log(len(dist))
            marg_dist_log.append(marg_log)

	# Grid search across the negative values
        I = 0

        while I > min_I:

            # Determines effect size
            I -= I_step
            pos_I.insert(0, I)
            dist = []

            # For each fraction of effective reagents, calculate and save likelihood
            for hit_rate in pos_hit_rate:

                # Finds likelihood of parameter combination
                like = like_fun(rhos, I, hit_rate, hit_like, back_likes,
                                        back_dist, off_likes, off_rate)
                dist.append(like)

            # Saves raw likelihoods for additional analysis
            all_dist[round(I, 2)] = dist

	    # Performs marginalization in log scale for stability
            marg_log = scipy.misc.logsumexp(dist) - math.log(len(dist))
            marg_dist_log.insert(0, marg_log)

	# Calculates likelihood fractions for interval calculation
        norm_marg_dist_log = marg_dist_log - scipy.misc.logsumexp(marg_dist_log)
        norm_marg_dist = numpy.exp(norm_marg_dist_log)

	# Sanity check
        if sum(norm_marg_dist) <= 0.95 or sum(norm_marg_dist) >= 1.05:
            print sum(norm_marg_dist)
            sys.exit('Error: Unstable computation')

	# Retrieves the 95% credible interval
        total_weight, min_ind, max_ind = findInterval(norm_marg_dist, 0.95,
							norm_marg_dist.argmax())

        # Finds maximum likelihood
        maxL = max(marg_dist_log)

        # Saves most likely parameters, along with likelihood ratio
        geneDist[gene] = all_dist
        geneI[gene] = pos_I[marg_dist_log.index(maxL)]

        # Calculates log-likelihood ratio
	geneL[gene] = 2 * (maxL - like0)
        geneInterval[gene] = (pos_I[min_ind], pos_I[max_ind], total_weight)

        # If enrichments are not provided, returns 0s for downstream analysis
        if not rhos:
            geneI[gene] = 0
            geneL[gene] = 0
            geneInterval[gene] = (0, 0, 1)

    return geneI, geneL, geneInterval, geneDist

###############################################################################
# Function to calculate likelihoods

def retrieveLikelihoods(gene_rhos, back_rhos, off_rhos, off_rate, like_fun,
                                nums, gene_span, I_step):

    # Checks if single processor
    if nums == 1:
        single = True
    else:
        single = False

    # Calculates gaussian kernel estimations of the background distributions
    back_dist = st.gaussian_kde(back_rhos)
    off_dist = st.gaussian_kde(off_rhos)

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
        fn = pp.Template(job_server, trial, (like_fun, findInterval,),
                ("numpy", "math", "scipy.misc"))

    # Single core version
    if single:
        geneI, geneL, geneInterval, geneDist = trial(gene_rhos, back_dist,
                                                        off_dist, off_rate,
                                                        like_fun, gene_span, I_step)

    # Parallel version
    if not single:

        geneI = {}
        geneL = {}
        geneInterval = {}
        geneDist = {}

        GeneRhosSplit = []
        keys = gene_rhos.keys()

	# Determines the number of genes in each task
        n = int(len(keys) / nums) + 1
        keysSplit = [keys[i: i + n] for i in range(0, len(keys), n)]

        # Splits dictionary
        for keyList in keysSplit:
            dictChunk = {}

            for key in keyList:
                dictChunk[key] = gene_rhos[key]

            GeneRhosSplit.append(dictChunk)

	# Submits individual jobs for each gene list
        jobs = []
        for splits in GeneRhosSplit:
            jobs.append(fn.submit(splits, back_dist, off_dist, off_rate, like_fun,
                                                gene_span, I_step))

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
# Calculates p-values with MW and KS statistics

def calculatePvals(gene_rhos, back_rhos):

    rhoMW = {}
    rhoKS = {}

    # For each gene, calls the Kolgormov-Smirnoff and Mann-Whitney U test
    for gene in gene_rhos:
        rhoKS[gene] = scipy.stats.mstats.ks_twosamp(numpy.ma.array(data=gene_rhos[gene],
								mask=False),
                            			    numpy.ma.array(data=back_rhos,
								mask=False))[1]

        rhoMW[gene] = scipy.stats.mstats.mannwhitneyu(numpy.ma.array(data=gene_rhos[gene], 									mask=False),
                            			    numpy.ma.array(data=back_rhos,
								mask=False))[1]
    return rhoMW, rhoKS


###############################################################################
# Creates parallel environment to calculate pvalues

def retrievePvals(gene_rhos, back_rhos, nums):

    # Checks if parallel process requested
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

    # Initiates parallel server
    if not single:
        job_server = pp.Server(nums)
        fn = pp.Template(job_server, calculatePvals, (),
                ("numpy", "scipy.stats.mstats"))

    # Single core version
    if single:
        rhoMW, rhoKS = calculatePvals(gene_rhos, back_rhos)

    # Parallel version
    if not single:

        rhoMW = {}
        rhoKS = {}

        GeneRhosSplit = []
        keys = gene_rhos.keys()

	# Determines the number of genes in each task
        n = int(len(keys)/nums)+1
        keysSplit = [keys[i:i + n] for i in range(0, len(keys), n)]

        # Splits dictionary
        for keyList in keysSplit:

            dictChunk = {}
            for key in keyList:
                dictChunk[key] = gene_rhos[key]

            GeneRhosSplit.append(dictChunk)

	# Submits individual jobs for each gene list
        jobs = []
        for splits in GeneRhosSplit:
            jobs.append(fn.submit(splits, back_rhos))

	# Retrieves results from individual jobs
        for job in jobs:
            val = job()

            try:
                rhoMWi, rhoKSi = val

            # Catches subprocess errors
            except TypeError:
                sys.exit('Subprocess failed')

            # Saves results
            rhoMW.update(rhoMWi)
            rhoKS.update(rhoKSi)

    return rhoMW, rhoKS


###############################################################################
# Calculates permutations for casTLE p-values

def retrievePerm(draw_num, perm_num, back_rhos, tar_rhos, off_rate, like_fun,
			nums, ref_file, gene2rat, I_step, erase):

    # Generates random 'genes' from all targeting elements
    perm_rhos = dict([(i, random.sample(tar_rhos, draw_num)) for i in range(perm_num)])

    # Calculates grid search
    gene_span = {}

    for gene, rhos in perm_rhos.items():

        max_I = 2 * max(rhos + [0])
        min_I = 2 * min(rhos + [0])

        gene_span[gene] = (min_I, max_I)

    # Finds the likelihood of these 'genes'
    permI, permL, permInterval, geneDist = retrieveLikelihoods(perm_rhos, back_rhos,
                                        tar_rhos, off_rate, like_fun,
                                        nums, gene_span, I_step)

    # Unless prompted, remember previous permutations
    if erase:
        print('Previous reference removed')
        os.remove(ref_file)

    # Write new permutations to file
    with open(ref_file, 'a') as ref_open:
        ref_csv = csv.writer(ref_open, delimiter=',', lineterminator='\n')

        for number in permI:
            ref_csv.writerow([permI[number], permL[number], permInterval[number]])

    all_perm_rat = []

    # Read back in both old and new permutations
    with open(ref_file, 'r') as ref_open:
        ref_csv = csv.reader(ref_open, delimiter=',', lineterminator='\n')

        for line in ref_csv:
            all_perm_rat.append(float(line[1]))

    # Retrieve actual gene log-likelihood ratios
    genes_rats = gene2rat.items()
    genes = [x[0] for x in genes_rats]
    rats = [x[1] for x in genes_rats]

    # Ranks actual genes based on permutations
    rat_rank = numpy.searchsorted(sorted(all_perm_rat), rats, 'right')
    gene_rank = dict(zip(genes, rat_rank))

    geneP = {}

    # Converts ranks to p values
    for gene in gene_rank:
        geneP[gene] = 1 - (gene_rank[gene] - 1) / float(len(all_perm_rat))

    return geneP, str(len(all_perm_rat))


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
        like0 = dist1[0][0] + dist2[0][0]

        # Initiate new grid
        min_I, max_I = gene_span[gene]
        marg_dist_log = [like0]
        pos_I = [0]

        # Search over positive effect values
        I = 0
        while I < max_I:
            I += I_step
            pos_I.append(I)

            dist = []

            # Retrieves the likelihood for each combination
            skip = True

            for E_like in ([dist1[0][0]] + dist1[round(I, 2)]):
                for T_like in ([dist2[0][0]] + dist2[round(I, 2)]):

                    if skip:
                        skip = False
                        continue

                    dist.append(E_like + T_like)

	    # Performs marginalization in log scale for stability
            marg_log = scipy.misc.logsumexp(dist) - math.log(len(dist))
            marg_dist_log.append(marg_log)

        # Search over negative effect values
        I = 0
        while I > min_I:
            I -= I_step
            pos_I.insert(0, I)

            dist = []

            # Retrieves the likelihood for each combination
            skip = True

            for E_like in ([dist1[0][0]] + dist1[round(I,2)]):
                for T_like in ([dist2[0][0]] + dist2[round(I,2)]):

                    if skip:
                        skip = False
                        continue

                    dist.append(E_like + T_like)

	    # Performs marginalization in log scale for stability
            marg_log = scipy.misc.logsumexp(dist) - math.log(len(dist))
            marg_dist_log.insert(0, marg_log)

	# Calculates likelihood fraction for interval calculation
        norm_marg_dist_log = marg_dist_log - scipy.misc.logsumexp(marg_dist_log)
        norm_marg_dist = numpy.exp(norm_marg_dist_log)

	# Sanity check
        if sum(norm_marg_dist) <= 0.95 or sum(norm_marg_dist) >= 1.05:
            sys.exit('Error: Unstable computation')

	# Retrieves the 95% confidence interval
        total_weight, min_ind, max_ind = findInterval(norm_marg_dist, 0.95,
							norm_marg_dist.argmax())

        # Finds maximum likelihood
        maxL = max(marg_dist_log)

        # Saves most likely parameters, along with likelihood ratio
        geneI[gene] = pos_I[marg_dist_log.index(maxL)]
        geneL[gene] = 2 * (maxL - like0)
        geneInterval[gene] = (pos_I[min_ind], pos_I[max_ind], total_weight)

    return geneI, geneL, geneInterval


###############################################################################
# Defines span for gene combo while unifying gene IDs

def comboSpan(gene_rhos1, gene_rhos2):

    # Retrieves ID maps
    geneID2Name, geneID2Info, geneName2ID, geneEns2Name = retrieveInfo()

    gene_span = {}
    add_gene_rhos1 = {}
    add_gene_rhos2 = {}

    for gene, rhos1 in gene_rhos1.items():

        # Converts IDs
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

        # Searches for ID in other screen
        if gene in gene_rhos2:
            rhos2 = gene_rhos2[gene]

        elif geneID in gene_rhos2:
            rhos2 = gene_rhos2[geneID]

        elif name in gene_rhos2:
            rhos2 = gene_rhos2[name]

        else:
            rhos2 = []

        # Initiates grid
        max_I = 2 * max(rhos1 + rhos2 + [0])
        min_I = 2 * min(rhos1 + rhos2 + [0])

        # Saves grid by unified ID
        gene_span[geneID] = (min_I, max_I)
        add_gene_rhos1[geneID] = rhos1
        add_gene_rhos2[geneID] = rhos2
        
    for gene, rhos2 in gene_rhos2.items():

        # Converts IDs
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

        # Checks if gene already added    
        if geneID in add_gene_rhos1:
            continue

        else:
            rhos1 = []

        # Initiates grid
        max_I = 2 * max(rhos1 + rhos2 + [0])
        min_I = 2 * min(rhos1 + rhos2 + [0])

        # Saves grid by unified ID
        gene_span[geneID] = (min_I, max_I)
        add_gene_rhos1[geneID] = rhos1
        add_gene_rhos2[geneID] = rhos2
        
    return add_gene_rhos1, add_gene_rhos2, gene_span


###############################################################################
# Function to calculate permutations for casTLE pvalues on combinations

def comboPerm(draw_num1, draw_num2, perm_num, back_rhos1, back_rhos2,
                                tar_rhos1, tar_rhos2, off_rate1, off_rate2,
                                like_fun1, like_fun2,
                                nums, ref_file, gene2rat, I_step, erase):

    # Creates 'genes'
    perm_rhos1 = dict([(i, random.sample(tar_rhos1, draw_num1)) for i in range(perm_num)])
    perm_rhos2 = dict([(i, random.sample(tar_rhos2, draw_num2)) for i in range(perm_num)])

    # Calculates grids in seperate function
    add_perm_rhos1, add_perm_rhos2, perm_span = comboSpan(perm_rhos1, perm_rhos2)

    # Calculates likelihoods seperately
    data1 = retrieveLikelihoods(add_perm_rhos1, back_rhos1,
                                        tar_rhos1, off_rate1,
                                        like_fun1, nums,
                                        perm_span, I_step)

    data2 = retrieveLikelihoods(add_perm_rhos2, back_rhos2,
                                        tar_rhos1, off_rate1,
                                        like_fun2, nums,
                                        perm_span, I_step)

    # Combines likelihoods
    permI, permL, permInterval = retrieveCombo(data1, data2, perm_span, I_step)

    # Unless prompted, remember previous permutations
    if erase:
        print('Previous reference removed')
        os.remove(ref_file)

    # Write new permutations to file
    with open(ref_file, 'a') as ref_open:
        ref_csv = csv.writer(ref_open, delimiter=',', lineterminator='\n')

        for number in permI:
            ref_csv.writerow([permI[number], permL[number], permInterval[number]])

    all_perm_rat = []

    # Read back in both old and new permutations
    with open(ref_file, 'r') as ref_open:
        ref_csv = csv.reader(ref_open, delimiter=',', lineterminator='\n')

        for line in ref_csv:
            all_perm_rat.append(float(line[1]))

    # Actual gene statistics
    genes_rats = gene2rat.items()
    genes = [x[0] for x in genes_rats]
    rats = [x[1] for x in genes_rats]

    # Ranks actual genes based on permutations
    rat_rank = numpy.searchsorted(sorted(all_perm_rat), rats, 'right')
    gene_rank = dict(zip(genes, rat_rank))

    geneP = {}

    # Converts ranks to p values
    for gene in gene_rank:
        geneP[gene] = 1 - (gene_rank[gene] - 1) / float(len(all_perm_rat))

    return geneP, str(len(all_perm_rat))

###############################################################################
