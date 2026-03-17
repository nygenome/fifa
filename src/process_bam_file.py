#!/bin/env python

##########################################################################
#	USAGE: extract features from BAM file
#   DESCRIPTION: Functions for extracting features for variants in VCF file
#   from user provided BAM file
##########################################################################


################################################################################
### MODULES ####################################################################

import pysam 
import pysamstats
import numpy as np # type: ignore
import pandas as pd # type: ignore
import pickle
import csv
import os, sys
import re
from statistics import mean, median
from scipy.stats import entropy # type: ignore
from Bio.SeqUtils import gc_fraction # type: ignore
from Bio import SeqIO # type: ignore
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from threading import Lock
import traceback
import subprocess
import multiprocessing as mp
from multiprocessing import Pool
import time
import logging

################################################################### /MODULES ###
################################################################################


global contig_sequences
global logger
lock = Lock()

def safe_mean(values):
    """
    Return mean of array accounting for None values or empty arrays
    """
    if values:
        filtered_values = [v for v in values if v is not None]
        return mean(filtered_values) if filtered_values is not None else None
    return None 

def safe_median(values):
    """
    Return median of array accounting for None values or empty arrays
    """
    if values:
        filtered_values = [v for v in values if v is not None]
        return median(filtered_values) if filtered_values is not None else None
    return None

def is_read_filtered(read):
    """
    Check if read passes minimal SAM flag filters
    and our immitation of the Mutect2 Filters
    """
    
    mutect_filtered = False
    sam_flags = False

    # Check SAM flags
    if read.is_unmapped or read.is_secondary or read.is_duplicate:
        sam_flags = True
        mutect_filtered = True
    elif not read.is_proper_pair or read.mate_is_unmapped:
        mutect_filtered = True

    # Check read properties
    elif read.query_length <= 1 or read.mapq <= 10 or read.mapq == 255 or read.has_tag('SA') or not read.cigar:
        mutect_filtered = True

    # Check CIGAR string
    elif not re.fullmatch(r'(?=.*[MD=X])([0-9]+[MIDSHP=X])+', read.cigarstring) or \
         re.search(r'(D|I)[0-9]*(?=(D|I))', read.cigarstring) or \
         re.match(r'^[0-9]*[SH]*[0-9]*D', read.cigarstring) or re.match(r'D[0-9]*[SH]*$', read.cigarstring):
        mutect_filtered = True

    # Check alignment conditions
    elif read.query_alignment_start < 0 or read.query_alignment_end < read.query_alignment_start or \
         read.infer_read_length() <= 0:
        mutect_filtered = True

    return sam_flags, mutect_filtered

def get_query_position(read, chromosome_position):
    """
    Extension of the pysam.AllignedSegment.get_reference_positions()
    method to account for possible indel errors. 

    Returns:
        int: variant position in the read. Returns None for deletions
    """    
    ref_poses = read.get_reference_positions()
    if chromosome_position in ref_poses:
        pos_in_read = ref_poses.index(chromosome_position)
        return pos_in_read
    return None


def process_cigar_tupples(read, reference_pos):
    """
    Process cigar string to get distance to effective 5' and 3' end
    and position of the variant in the read (both with and without
    counting clipped bases). 
    Iterates through cigar string operators, while keeping track of the corresponding
    position in the reference sequence. 

    Args:
        read (pysam.AllignedSegment)
        reference_pos (int): position of the variant on the reference sequence

    Returns:
        int: position of the variant in the read, excluding hard clipped bases
        int: position of the variant in the read, including hard clipped bases
        int: distance from variant to read effective 5' end
        int: distance from variant to read effective 3' end
        int: length of clipped bases (soft and hard clipped)
        
        Returns None if read is missing cigarstring or if the variant is a deletion
    """    
    if not(read.cigarstring):
        return (None, None, None, None, None)
    
    ref_pos_in_read = read.reference_start
    position_in_query = 0 # position of the variant in the read, excluding hard clipped bases
    real_position_in_read = 0 # position of the variant in the read, including hard clipped bases
    distance_to_5prime = distance_to_3prime = 0 
    
    clipped_length = 0
    
    for operation, length in read.cigartuples:
        if operation == 0: # Match or mismatch
            if ref_pos_in_read <= reference_pos:  
                if reference_pos < ref_pos_in_read + length:
                    real_position_in_read += reference_pos - ref_pos_in_read
                    position_in_query += reference_pos - ref_pos_in_read
                    distance_to_5prime += reference_pos - ref_pos_in_read
                    distance_to_3prime += length - (reference_pos - ref_pos_in_read)
                else: 
                    real_position_in_read += length
                    position_in_query += length
                    distance_to_5prime += length
            else:
                distance_to_3prime += length 
            ref_pos_in_read += length

        elif operation == 1:  # Insertion
            if ref_pos_in_read <= reference_pos:
                distance_to_5prime += length
                real_position_in_read += length
                position_in_query += length
            else:
                distance_to_3prime += length 
        
        elif operation == 2:  # Deletion
            if ref_pos_in_read <= reference_pos < ref_pos_in_read + length:
                return (None, None, None, None, None) 
            
            ref_pos_in_read += length
        
        elif operation == 4:  # Soft clipping
            clipped_length += length
            if ref_pos_in_read <= reference_pos: 
                real_position_in_read += length
                position_in_query += length

        elif operation == 5:  # Hard clipping
            clipped_length += length
            if ref_pos_in_read <= reference_pos: 
                real_position_in_read += length
    
    return position_in_query, real_position_in_read, distance_to_5prime, distance_to_3prime, clipped_length

def get_mismatch_and_insertion_positions(read):
    """
    Counts total number of insertions plus mismatches, 
    and extracts base qualities at those positions with 
    pysam.AllignedSegment.query_qualities

    Returns:
        int: number of mismatches in the read
        list: base qualities of bases where a mismatch occurs. 
    """    
    mismatch_positions = []
    current_position = 0
    
    # Parse the MD tag for mismatches
    md_tag = read.get_tag('MD')
    i = 0
    while i < len(md_tag):
        if md_tag[i].isdigit():
            num = ''
            while i < len(md_tag) and md_tag[i].isdigit():
                num += md_tag[i]
                i += 1
            current_position += int(num)
        elif md_tag[i] == '^':
            while i < len(md_tag) and not md_tag[i].isdigit():
                i += 1
        else:
            mismatch_positions.append(current_position)
            current_position += 1
            i += 1
            
    mismatch_base_quals = [read.query_qualities[pos] for pos in mismatch_positions] if mismatch_positions else None
    
    # Parse the CIGAR string to count insertions in the total mismatches count
    cigar = read.cigarstring
    num_mismatches = len(mismatch_positions)
    i = 0
    while i < len(cigar):
        num = ''
        while cigar[i].isdigit() and i < len(cigar):
            num += cigar[i]
            i += 1
        if cigar[i] == 'I':
            for x in range(1, int(num)+1):
                num_mismatches += 1
            i += 1
        else:
            i += 1
    return num_mismatches, mismatch_base_quals

def get_ref_seq(fasta_file):
    contig_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_sequences[record.id] = record.seq
    return contig_sequences

## This feels like it could be done much more efficiently. 
## I essentially need the tri- and penta- nucleotide sequence around the variant 
## but the reference sequence. 

## Why am I reading everything into memory here... this seems inefficient.

def one_hot_encode(contig, position, contig_sequences):
    """
    Hot encode 5 base region around variant

    Args:
        contig, position (int)
        contig_sequences (dict): dict from parsing the reference fasta file

    Returns:
        array: hot encoded array of the 5 bases around variant
    """    

    seq = []
    ref_sequence = contig_sequences[contig]
    for i in range(-2,3):
        index = position + i
        if 0 <= index < len(ref_sequence):
            seq.append(ref_sequence[index])
        else: 
            seq.append('N') # For cases where 5-base window exceeds range of the chromosome
    
    mapping = dict(zip("ACGTN", range(5)))
    
    one_hot_encoded = np.zeros((5, 5)) 
    
    for i, base in enumerate(seq):
        if base in mapping:
            one_hot_encoded[i] = np.eye(5)[mapping[base]]
    
    return one_hot_encoded


################################################################################
###GENERATE BASE METRICS########################################################

def update_metrics(metrics, position, metrics_data, prefix, tumor_other_bases_count, num_total_reads):
    """
    Factored code for updating pooled dictionary of features for all variants (metrics). 
    Method takes in a dictionary of arrays of counts for the ref or var allele, and 
    handles mean / median calculations for features of both alleles. 
    Will append pooled dictionary with newly calculated features. 

    Args:
        metrics (dict): pooled dictionary of features for all variants. 
        position (str): variant position (key for dictionary metrics)
        metrics_data (dict): dictionary of lists of extracted features for the ref/var allele. 
        prefix (str): var or ref prefix
        tumor_other_bases_count (int): count of other bases (not var/ref) at that position
        num_total_reads (int): number of total reads mapping to that position
    """    

    metrics[position].update({
        f'{prefix}_avg_base_quality': safe_mean(metrics_data.get('base_qualities', [])),
        f'{prefix}_med_base_quality': safe_median(metrics_data.get('base_qualities', [])),
        f'{prefix}_avg_mapping_quality': safe_mean(metrics_data.get('mapping_quality', [])),
        f'{prefix}_med_mapping_quality': safe_median(metrics_data.get('mapping_quality', [])),
        f'{prefix}_avg_se_mapping_quality': safe_mean(metrics_data.get('single_end_mapq', [])),
        f'{prefix}_avg_num_mismatches_as_fraction': safe_mean(metrics_data.get('avg_num_mismatches', [])),
        f'{prefix}_avg_sum_mismatch_qualities': safe_mean(metrics_data.get('avg_sum_mismatch_base_quals', [])),
        f'{prefix}_med_frag_len': safe_median(metrics_data.get('read_frag_length', [])),
        f'{prefix}_avg_clipped_length': safe_mean(metrics_data.get('clipped_length', [])),
        f'{prefix}_avg_pos_as_fraction': safe_mean(metrics_data.get('avg_pos_as_fraction', [])),
        f'{prefix}_avg_distance_to_effective_3p_end': safe_mean(metrics_data.get('distances_to_3p_end', [])),
        f'{prefix}_avg_distance_to_effective_5p_end': safe_mean(metrics_data.get('distances_to_5p_end', [])),
        f'{prefix}_med_distance_to_effective_3p_end': safe_median(metrics_data.get('distances_to_3p_end', [])),
        f'{prefix}_med_distance_to_effective_5p_end': safe_median(metrics_data.get('distances_to_5p_end', [])),
        f'{prefix}_normmed_distance_to_effective_3p_end': safe_median(metrics_data.get('distances_to_3p_end', [])) / safe_median(metrics_data.get('read_frag_length', [])),
        f'{prefix}_normmed_distance_to_effective_5p_end': safe_median(metrics_data.get('distances_to_5p_end', [])) / safe_median(metrics_data.get('read_frag_length', [])),
        'tumor_other_bases_count': tumor_other_bases_count,
        'num_total_reads': num_total_reads,
        f'{prefix}_num_minus_strand': metrics_data.get('num_minus_strand', 0),
        f'{prefix}_num_plus_strand': metrics_data.get('num_plus_strand', 0)
    })

def get_base_metrics(bamfile, rec, ref_seq, FFPE, label):
    global contig_sequences

    if len(dict(rec.samples)[FFPE].get('AD')) <= 1 :
        logger.info(f'SKIPPING {position} because len(AD) < 1')
        return {}

    metrics={}

    position = '{0}:{1}_{2}>{3}'.format(rec.chrom, rec.pos, rec.ref, rec.alts[0])

    num_vars_filtered = 0 
    tumor_other_bases_count = 0
    num_total_reads = 0 

    if label[0] and label[1]:
        metrics[position]= {
            'Label' : 1 if rec.info.get(label[0]) == label[1] else 0
        }
    else:
        metrics[position]= {}

    sample_data = dict(rec.samples)[FFPE]    
    metrics[position].update({
        'tumor_ref_count': 0,
        'tumor_var_count': 0,
        'tumor_depth' : sample_data.get('DP') if sample_data.get('DP') else sum(sample_data.get('AD'))
        #'tumor_VAF' : dict(rec.samples)[FFPE].get('AF') if isinstance(dict(rec.samples)[FFPE].get('AF'), float) else dict(rec.samples)[FFPE].get('AF')[0] 
    })
        
    if 'AF' in sample_data.keys():
        metrics[position]['tumor_VAF'] = sample_data.get('AF') if isinstance(sample_data.get('AF'), float) else sample_data.get('AF')[0] 
    else:
        metrics[position]['tumor_VAF'] = None


    ## Should change because we're not using hot encodings anymore
    ## And also change to extract tri-/penta-nucleotide sequence
    one_hot_encoded = one_hot_encode(rec.chrom, rec.start, contig_sequences)
    keys = ['left_two_base', 'left_one_base', 'hot_encoded_ref_base', 'right_one_base', \
            'right_two_base']
    for key, value in zip(keys, one_hot_encoded):
        metrics[position][key] = value
    
    mapping = dict(zip("ACGTN", range(5)))
    metrics[position]['hot_encoded_var_base'] = np.eye(5)[mapping[rec.alts[0]]]
    
    ## Define the metrics
    ref_metrics = {
        'distances_to_3p_end': [], 'distances_to_5p_end': [], 'base_qualities': [],
        'mapping_quality': [], 'single_end_mapq': [], 'read_frag_length': [],
        'clipped_length': [], 'avg_pos_as_fraction': [], 'avg_num_mismatches': [],
        'avg_sum_mismatch_base_quals': [], 'num_minus_strand': 0, 'num_plus_strand': 0
    }
        
    var_metrics = {
        'distances_to_3p_end': [], 'distances_to_5p_end': [], 'base_qualities': [],
        'mapping_quality': [], 'single_end_mapq': [], 'read_frag_length': [],
        'clipped_length': [], 'avg_pos_as_fraction': [], 'avg_num_mismatches': [],
        'avg_sum_mismatch_base_quals': [], 'num_minus_strand': 0, 'num_plus_strand': 0
    }

    ## Suggestion from André : 
    ## class for these stats - make a class/ object  for that dictionary that I just initialize
    ## attribute not use string for the key
    ## move mean / median methods to class 
    
    for pileupcolumn in bamfile.pileup(contig=rec.chrom, start=rec.pos - 1, stop=rec.pos, min_base_quality=0, min_mapping_quality=0):
        if pileupcolumn.pos == rec.pos - 1:
            for pileupread in pileupcolumn.pileups:
                query_index = pileupread.query_position
                read = pileupread.alignment
                
                pos_in_read, read_index, distance_5prime, distance_3prime, clipped_length = \
                    process_cigar_tupples(read, rec.start)
                
                #Handle Edge cases
                if pos_in_read != query_index:
                    logger.info(f'{FFPE} {position}: Please verify why manually calculated position in read does not match '
                            f"pysam's query position\n {pos_in_read} != {query_index}")
                    continue

                if len(rec.alts[0]) <= 0:  
                    logger.info(f'{FFPE} sample has no alternative alleles at: {position}')
                    continue
                
                if query_index is None or query_index < 0 :
                    print(f'{FFPE} {position}: Query index is None or < 0')
                    # Still never understood why that is. Probably has to do with how Mutect2 processes reads
                    # vs. Pysam (but this can be something to ask Will about)
                    continue

                base = read.query_sequence[query_index]
                is_ref = base == rec.ref
                is_var = base == rec.alts[0]

                num_total_reads += 1

                if is_read_filtered(read)[1]:
                    num_vars_filtered += 1 

                if is_ref or is_var:
                    metrics[position]['tumor_ref_count' if is_ref else 'tumor_var_count'] += 1
                    
                    # Will update features for the ref/var allele accordingly 
                    metrics_list = ref_metrics if is_ref else var_metrics 

                    metrics_list['base_qualities'].append(read.query_qualities[query_index])

                    if read.is_reverse:
                        metrics_list['num_minus_strand'] += 1
                    else:
                        metrics_list['num_plus_strand'] += 1 

                    metrics_list['read_frag_length'].append(read.infer_read_length())

                    metrics_list['avg_pos_as_fraction'].append(read_index / (read.infer_read_length() / 2))
                    metrics_list['distances_to_5p_end'].append(distance_5prime) 
                    metrics_list['distances_to_3p_end'].append(distance_3prime)
                    metrics_list['clipped_length'].append(clipped_length)

                    
                    if read.is_paired:
                        metrics_list['mapping_quality'].append(read.mapping_quality)
                    else:
                        metrics_list['single_end_mapq'].append(read.mapping_quality)

                    if read.has_tag('MD'):
                        num_mismatches, mismatch_base_quals = get_mismatch_and_insertion_positions(read)
                        metrics_list['avg_num_mismatches'].append(num_mismatches / read.infer_read_length())
                        if mismatch_base_quals:
                            metrics_list['avg_sum_mismatch_base_quals'].append(sum(mismatch_base_quals))      
                else:
                    tumor_other_bases_count += 1
            
        
    if metrics[position]['tumor_ref_count'] > 0:
        update_metrics(metrics, position, ref_metrics, 'tumor_ref', \
                        tumor_other_bases_count, num_total_reads)
    
    if metrics[position]['tumor_var_count'] > 0:
        update_metrics(metrics, position, var_metrics, 'tumor_var', \
                        tumor_other_bases_count, num_total_reads)
        
        #this is a typo oups
        tumor_var_num_plus_strand = var_metrics.get('num_minus_strand', 0)
        tumor_var_num_minus_strand = var_metrics.get('num_plus_strand', 0)

        # Calculate FDeamC score for deamination variants only
        if tumor_var_num_plus_strand + tumor_var_num_minus_strand > 0: 
            if (rec.ref, rec.alts[0]) == ('C', 'T'):
                metrics[position]['FDeamC'] = float(tumor_var_num_plus_strand /\
                                                    (tumor_var_num_plus_strand + \
                                                        tumor_var_num_minus_strand))
            elif (rec.ref, rec.alts[0]) == ('G', 'A'):
                metrics[position]['FDeamC'] = float(tumor_var_num_minus_strand /\
                                                    (tumor_var_num_plus_strand + \
                                                        tumor_var_num_minus_strand))
            else:
                metrics[position]['FDeamC'] = 0
            
            metrics[position]['SOB'] = float((tumor_var_num_plus_strand - tumor_var_num_minus_strand) \
                                                / (tumor_var_num_plus_strand + tumor_var_num_minus_strand))
        else:
            metrics[position]['FDeamC'] = 0
            metrics[position]['SOB'] = 0
        
    ref_base_qualities = ref_metrics.get('base_qualities', None)
    var_base_qualities = var_metrics.get('base_qualities', None)
        
    if safe_median(ref_base_qualities) and safe_median(var_base_qualities):
        metrics[position]['tumor_ref_base_quality_frac'] = \
        safe_median(var_base_qualities) / safe_median(ref_base_qualities)
    else:
        metrics[position]['tumor_ref_base_quality_frac'] = 0
                
    metrics[position]['tumor_reads_filtered'] = num_vars_filtered
    
    if metrics[position]['tumor_VAF'] is None:
        metrics[position]['tumor_VAF'] = metrics[position]['tumor_var_count'] / metrics[position]['tumor_depth'] if metrics[position]['tumor_depth'] > 0 else 0 
    
    return metrics

########################################################GENERATE BASE METRICS###
################################################################################

################################################################################
###GENERATE WINDOW-BASED METRICS################################################

def get_coverage_in_window(bamfile, fastafile, chrom, left, right):
    """
    Extract coverage using BAM file for a specified window in the chromosome.
    Window is denoted by left and right endpoints, and cannot exceed the range of the 
    contig in the reference sequence.

    If left and right endpoints are outside of contig range, function select instead
    accordingly to the leftmost and rightmost ends of the contig, narrowing the window range. 

    Args:
        left (int): leftmost end of the window
        right (int): rightmost end of the window

    Returns:
        int: median coverage in the specified window
        float: standard deviation of coverege in thewindow
    """    
    window_coverage = pysamstats.load_pileup(
        type="coverage_ext", alignmentfile=bamfile, fafile=fastafile, 
        chrom=chrom, start=left, end=right, truncate=False
    )
    
    # Ensure that left and right endpoints meet contig ranges
    indeces = np.where((window_coverage.pos >= left) & (window_coverage.pos <= right))[0]
    index_left, index_right = (indeces[0], indeces[-1]) if indeces.size > 0 else (0, 0)
    
    window_data = window_coverage[index_left:index_right].reads_all
    
    if window_data.size <= 1:
        return None, None
    
    return np.median(window_data), np.var(window_data)
    
def get_read_fractions(bamfile, chrom, left, right):
    """
    Extracts counts of reads that are duplicated,
    poor mapq, or failing other features within a
    specified window around the variant.
    Leftmost and rightmost positions denote the ends of 
    the basepair window. 

    Args:
        left (int): leftmost end of the window
        right (int): rightmost end of the window

    Returns:
        int: median fragment lengths in the window
        int: fraction of duplicated reads
        int: franction of reads that are improperly paired
        int: median mapq of reads in the window
        int: fraction of reads that are filtered by mutect2
    """    
    frag_lenths, mapqs = [], []
    num_mapq0 = num_duplicates = num_improper_paired = num_total_reads = 0
    num_filtered_mutect = num_filtered_SAM = 0
    
    for read in bamfile.fetch(chrom, left, right):
        num_total_reads += 1
        
        filtered_SAM_flags, filtered_mutect2 = is_read_filtered(read)
        num_filtered_SAM += filtered_SAM_flags
        
        if read.query_alignment_end > bamfile.get_reference_length(chrom) or filtered_mutect2:
            num_filtered_mutect += 1
        else:
            frag_lenths.append(abs(read.template_length))
            mapqs.append(read.mapping_quality)
            
        num_mapq0 += read.mapping_quality == 0
        num_improper_paired += not read.is_proper_pair
        num_duplicates += read.is_duplicate
    
    
    return (
        safe_median(frag_lenths),
        num_duplicates / num_total_reads,
        num_mapq0 / num_total_reads,
        num_improper_paired / num_total_reads,
        safe_median(mapqs),
        num_filtered_mutect / num_total_reads
    )

def get_metrics_500_bp_window(rec, bamfile, fastafile): 
    """
    Calculate metrics in the +/-500 bp window around a given variant. 

    Args:
        rec (pysam.VariantRecord): a single variant recorded in the VCF file
        bamfile (pysam.AlignmentFile): user inputted BAM file
        vcffile (pysam.VariantFile): user inputted VCF file
        fastafile (pysam.FastaFile): user inputted reference sequence 

    Returns:
        dict: all windowed features extracted for a single variant record from the VCF file
    """    

    position_metrics = {} 
    chrom = rec.chrom
    left = max(rec.start - 500, 0)
    right = min(rec.stop + 500, bamfile.get_reference_length(chrom))

    left_window = [max(0, left - 500), left]
    right_window = [right, min(right + 500, bamfile.get_reference_length(chrom))]
    
    alignment_seq = fastafile.fetch(chrom, left, right)

    position_metrics['window_gc_cont'] = gc_fraction(alignment_seq)
    position_metrics['window_seq_entropy'] = entropy(np.unique(list(alignment_seq), return_counts=True)[1] / len(alignment_seq), base=2)
    
    position_metrics['window_median_cov'], position_metrics['window_cov_variance'] = get_coverage_in_window(bamfile, fastafile, chrom, left, right)
    
    if position_metrics['window_median_cov'] is None or position_metrics['window_median_cov'] == 0:
        logger.info("COVERAGE ERROR. No coverage in flanking window.")
        position_metrics['window_min_cov_ratio'] = position_metrics['window_max_cov_ratio'] = None
    else: 
        left_window_median_cov = get_coverage_in_window(bamfile, fastafile, chrom, *left_window)[0]
        right_window_median_cov = get_coverage_in_window(bamfile, fastafile, chrom, *right_window)[0]
        
        left_window_median_cov = left_window_median_cov if left_window_median_cov is not None else 0
        right_window_median_cov = right_window_median_cov if right_window_median_cov is not None else 0

        window_coverage_ratios = (left_window_median_cov / position_metrics['window_median_cov'],
                              right_window_median_cov / position_metrics['window_median_cov'])
    
        position_metrics['window_min_cov_ratio'] = min(window_coverage_ratios)
        position_metrics['window_max_cov_ratio'] = max(window_coverage_ratios)
    
    (position_metrics['window_median_frag_len'], 
     position_metrics['window_dup_frac'], 
     position_metrics['window_multi_frac'], 
     position_metrics['window_improper_frac'], 
     position_metrics['window_median_mapq'], 
     position_metrics['window_read_filter_frac']) = get_read_fractions(bamfile, chrom, left, right)
    
    return position_metrics

################################################################################
###############################################GENERATE WINDOW-BASED METRICS####

def get_chromosome_metrics(chrom, vcf_path, bam_path, ref_seq, sample, label):
    window_metrics = {}
    base_metrics = {}
    bamfile = pysam.AlignmentFile(bam_path, "rb", reference_filename=ref_seq)
    fastafile = pysam.FastaFile(ref_seq)
    vcffile = pysam.VariantFile(vcf_path) 

    for rec in vcffile.fetch(contig=chrom):

        position = '{0}:{1}_{2}>{3}'.format(rec.chrom, rec.pos, rec.ref, rec.alts[0])
        if rec.ref not in ['A', 'C', 'T', 'G'] or rec.alts[0] not in ['A', 'C', 'T', 'G']:
            continue

        window_metrics[position] = get_metrics_500_bp_window(rec, bamfile, fastafile)
        base_metrics.update(get_base_metrics(bamfile, rec, ref_seq, sample, label))

    vcffile.close()
    bamfile.close()
    fastafile.close()
    return base_metrics, window_metrics

#################################################################################
### SPLIT VARIANTS BY CHROMOSOME ################################################

def split_vars_by_chromosome(bam_path, vcf_path, ref_seq, FFPE, label, num_processes):
    logger.info("Grouping Variants by chromosome...")
    async_results = []
    
    is_compressed = vcf_path.endswith('.gz')
    has_index = os.path.exists(vcf_path + '.csi') or os.path.exists(vcf_path + '.tbi')
    if not is_compressed:
        logger.error("VCF File is not compressed. Please bgzip the VCF file before running feature extraction.")
        exit(1)
    if not has_index:
        logger.error("VCF File is not indexed. Please index the VCF file before running feature extraction.")
        exit(1) 
    compressed_vcf = vcf_path
    
     ## Change here reminder! 
    vcffile = pysam.VariantFile(compressed_vcf) 
    
    pool = Pool(processes=num_processes) 
    ## This here is the point I'm confused about

    for contig in vcffile.header.contigs:
        async_result = pool.apply_async(get_chromosome_metrics, (contig, compressed_vcf, bam_path, ref_seq, FFPE, label))
        async_results.append(async_result)
    vcffile.close()
    pool.close()
    pool.join()
    
    base_metrics = {}
    window_metrics = {}
    
    # Merge together all results from the different chromosomes
    for async_result in async_results:
        base, window = async_result.get()
        base_metrics.update(base)
        window_metrics.update(window)
    
    num_vars = len(window_metrics.keys())
    logger.info(f"Finished metrics for {num_vars} in {FFPE}")

    ## remove here
    # os.remove(compressed_vcf)
    # if remove_index:
    #    os.remove(index_file)

    return base_metrics, window_metrics


    ## Old code to gzip and index VCF file

    '''
    logger.info("compressing VCF")
    os.makedirs('/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/temp', exist_ok=True) # Need to change this

    if vcf_path.endswith('.gz'):
        compressed_vcf = os.path.join('/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/temp/', os.path.basename(vcf_path))
        print(compressed_vcf)
        os.system(f'cp {vcf_path} {compressed_vcf}')  
    else:
        compressed_vcf = os.path.join('/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/temp/', os.path.basename(vcf_path) +'.gz')
        print(compressed_vcf)
        subprocess.run(f"bcftools view -O z -o {compressed_vcf} {vcf_path}", shell=True)
    logger.info("VCF File Compressed") 

    index_file = f"{compressed_vcf}.csi" 
    if not os.path.exists(index_file):
        remove_index = True
        subprocess.run(f"bcftools index {compressed_vcf}", shell=True)
        logger.info("VCF File Indexed")
    else:
        remove_index = False
        logger.info("Index file already exists")
    '''

    #pysam.tabix_compress(vcf_path, os.path.join(vcf_path, '.gz'), force=TRUE)
    #index = pysam.tabix_index(vcf_path, preset="vcf", force=True, keep_original=True)


#################################################################################
################################################ SPLIT VARIANTS BY CHROMOSOME ###


def get_mobster_tail_scores(sample, vcf_path):
    prefix = os.path.basename(vcf_path).replace(".vcf.gz", "")
    outfile = os.path.join(os.getcwd(), f'{prefix}_mobster.csv')
    mobster_scores = {}

    subprocess.run(
        ['Rscript', 'run_mobster.R', sample, vcf_path, outfile],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
        )
    
    if not os.path.isfile(outfile):
        logger.info(f"MOBSTER did not run succesfully on {sample}")
        vcffile = pysam.VariantFile(vcf_path)
        for rec in vcffile.fetch():
            position = '{0}:{1}_{2}>{3}'.format(rec.chrom, rec.pos, rec.ref, rec.alts[0])
            if rec.ref not in ['A', 'C', 'T', 'G'] or rec.alts[0] not in ['A', 'C', 'T', 'G']:
                continue
            mobster_scores[position] = {'Tail': 1}
        
    else:    
        with open(outfile, newline='') as mfile:
            reader = csv.reader(mfile, delimiter=',')
            header = next(reader)
            for row in reader: 
                sample, chrom, pos, REF, ALT, Tail = row
                position = '{0}:{1}_{2}>{3}'.format(chrom, pos, REF, ALT)
                if REF not in ['A', 'C', 'T', 'G'] or ALT not in ['A', 'C', 'T', 'G']:
                    continue
                mobster_scores[position] = {'Tail': Tail}
        logger.info(f"Finished MOBSTER calculations for {sample}")
    
    if os.path.isfile(outfile):
        os.remove(outfile)

    return mobster_scores

def extract_all_features(bam_path, vcf_path, ref_seq, sample, label, num_threads):    
    global contig_sequences 
    contig_sequences = get_ref_seq(ref_seq)

    start = time.time()
    
    with ThreadPoolExecutor(max_workers=3) as executor:
        future_window_base_metrics = executor.submit(split_vars_by_chromosome, bam_path, vcf_path, ref_seq, sample, label, num_threads)
        future_mobster_metrics = executor.submit(get_mobster_tail_scores, sample, vcf_path)
        
        base_metrics, window_metrics = future_window_base_metrics.result()
        mobster_metrics = future_mobster_metrics.result()

    end = time.time() - start

    logger.info(f"Finished all metrics for {sample} in {end}")
    return  mobster_metrics, window_metrics, base_metrics


#################################################################################
### PROCESS INPUT FILES #########################################################

def get_column_idx(header):
    # Find the index of the columns matching the patterns and store the results in a dictionary
    patterns = {
        'SAMPLE': re.compile(r'\s*(FFPE|Sample)\s*', re.IGNORECASE),
        'COHORT': re.compile(r'COHORT', re.IGNORECASE),
        'VCF_FILE': re.compile(r'VCF[_\s]*FILE', re.IGNORECASE),
        'BAM_FILE': re.compile(r'BAM[_\s]*FILE', re.IGNORECASE),
        'REF_SEQ': re.compile(r'REF[_\s]*SEQ', re.IGNORECASE)
    }
    columns = {}
    for keyword, pattern in patterns.items():
        columns[keyword] = next((i for i, col in enumerate(header) if pattern.match(col)), -1)
    
    matched_columns = [
        columns['SAMPLE'],
        columns['COHORT'],
        columns['VCF_FILE'],
        columns['BAM_FILE'],
        columns['REF_SEQ']
    ]
    if -1 in matched_columns:
        logger.error("Error: One or more required columns not found in pairs file.\nRequired columns are: \nFFPE, Cohort, VCF_FILE, BAM_FILE, REF_SEQ")
        exit(1)
    return matched_columns
    
def process_row(sample, cohort, vcf_path, bam_path, ref_seq, output_file, label, num_threads):   
    try:
        if os.path.isfile(vcf_path) and os.path.isfile(bam_path) and os.path.isfile(ref_seq):
            logger.info(f'Processing BAM file for sample: {sample}')
            flattened_data = []
            mobster_metrics, window_metrics, base_metrics = extract_all_features(bam_path, vcf_path, ref_seq, sample, label, num_threads)
        
            metrics = {
                variant: {**base_metrics.get(variant, {}), **window_metrics.get(variant, {}), **mobster_metrics.get(variant,{})}
                for variant in set(base_metrics) & set(window_metrics) & set(mobster_metrics)
            }
            ## Minor note to myself, but sets are unordered in python, so set(base_metrics), set(window_metrics), etc 
            ## will not preserve the order the variants were in (which is in the order they appear on the chromosome). 
            ## To modify this, I can look into : from collections import OrderedDict

            for variant, metric in metrics.items(): 
                row = {'Cohort': cohort, 'Sample': sample, 'Variant': variant}
                row.update(metric)
                flattened_data.append(row)

            df = pd.DataFrame(flattened_data)
            with lock:
                ## Check how fast this is vs. just using a print statement
                ## Minimal compute / be as fast as possible with lock

                df.to_csv(output_file, index=False)
                logger.info(f"Features extracted and stored: \n{output_file}")
        else:
            logger.error("There is an issue with one of your input files.")
            if not (os.path.isfile(vcf_path)):
                logger.error(f" The error is with your VCF File path: {vcf_path}")
            if not (os.path.isfile(bam_path)):
                logger.error(f" The error is with your BAM File path: {bam_path}")
            if not (os.path.isfile(ref_seq)):
                logger.error(f" The error is with your reference seq: {ref_seq}")
            raise FileNotFoundError("One or more input files are missing. Please verify that the paths are correct.")
    except Exception as e:
        logger.error(f"Error processing {sample}: {e}")
        traceback.print_exc()

def process_bam_file(outpath, label, num_threads, sample, vcf_file, 
bam_file, cohort=None, ref_seq=None):
    global logger
    logger = logging.getLogger(__name__)
    try:
        if os.path.isfile(outpath):
            output_file = outpath 
            
            if not (outpath.endswith('_extracted_features.csv')):
                logger.info("WARNING: If output path is a file, it must end with '_extracted_features.csv' due to dependencies "
                "in other parts of the code.")
                base_path, _ = os.path.splitext(outpath)
                output_file = base_path + "_extracted_features.csv"
                logger.info("New Output file : " + output_file)
        else:
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            output_file=os.path.join(outpath, f"{sample}_extracted_features.csv")
            
        process_row(sample, cohort, vcf_file, bam_file, ref_seq, output_file, label, num_threads)
    except Exception as e:
        logger.error(f"An error occurred: {e}")
    pass

#################################################################################
######################################################### PROCESS INPUT FILES ###