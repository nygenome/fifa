#!/usr/bin/env python

##########################################################################
#    USAGE: Extract features from BAM file for a set of variants in a VCF file,
#    with most-up-to-date paralellization scheme
##########################################################################

################################################################################
### MODULES ####################################################################

import pysam 
import pysamstats
import numpy as np # type: ignore
import pandas as pd # type: ignore
import csv
import os, sys
import re
from scipy.stats import entropy # type: ignore
from Bio.SeqUtils import gc_fraction # type: ignore
from concurrent.futures import ProcessPoolExecutor
import traceback
import subprocess
from multiprocessing import Queue, Lock, Process, Pool, Manager
import time
import logging
from metrics_dictionary import MetricsDictionary
from metrics_dictionary import safe_median

################################################################### /MODULES ###
################################################################################

global logger

def is_read_filtered(read):
    """
    Check if read passes
    our immitation of the Mutect2 Filters
    """
    
    mutect_filtered = False

    # Check SAM flags
    if read.is_unmapped or read.is_secondary or read.is_duplicate:
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

    return mutect_filtered

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
        return (None, None, None, None)
    
    ref_pos_in_read = read.reference_start
    real_position_in_read = 0 # position of the variant in the read, including hard clipped bases
    distance_to_5prime = distance_to_3prime = 0 
    
    clipped_length = 0
    
    for operation, length in read.cigartuples:
        if operation == 0: # Match or mismatch
            if ref_pos_in_read <= reference_pos:  
                if reference_pos < ref_pos_in_read + length:
                    real_position_in_read += reference_pos - ref_pos_in_read
                    distance_to_5prime += reference_pos - ref_pos_in_read
                    distance_to_3prime += length - (reference_pos - ref_pos_in_read)
                else: 
                    real_position_in_read += length
                    distance_to_5prime += length
            else:
                distance_to_3prime += length 
            ref_pos_in_read += length

        elif operation == 1:  # Insertion
            if ref_pos_in_read <= reference_pos:
                distance_to_5prime += length
                real_position_in_read += length
            else:
                distance_to_3prime += length 
        
        elif operation == 2:  # Deletion
            if ref_pos_in_read <= reference_pos < ref_pos_in_read + length:
                return (None, None, None, None) 
            
            ref_pos_in_read += length
        
        elif operation == 4:  # Soft clipping
            clipped_length += length
            if ref_pos_in_read <= reference_pos: 
                real_position_in_read += length

        elif operation == 5:  # Hard clipping
            clipped_length += length
            if ref_pos_in_read <= reference_pos: 
                real_position_in_read += length
    
    return real_position_in_read, distance_to_5prime, distance_to_3prime, clipped_length

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

    coverage = pysamstats.load_coverage(bamfile, chrom=chrom, start=left, end=right, truncate=True)
    num_total_reads = coverage['reads_all'].sum()
    ## In theory according to the documentation you should be able to get the num_mapq0, etc from 
    ## pysamstats, but for some reason I can't find those attributes in the recarray

    if num_total_reads == 0:
        return (0,0,0,0,0)
    
    num_improper_paired = num_total_reads - coverage['reads_pp'].sum()
    frag_lenths, mapqs = [], []
    num_mapq0 = num_duplicates = num_filtered_mutect = 0
    
    for read in bamfile.fetch(chrom, left, right):        
        if read.query_alignment_end > bamfile.get_reference_length(chrom) or is_read_filtered(read):
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

################################################################################
###############################################GENERATE WINDOW-BASED METRICS####

#################################################################################
### SPLIT VARIANTS FOR PROCESSING ###############################################

def process_variant(queue, sample, cohort, bam_path, ref_seq, iolock, final_dictionary):
    bamfile = pysam.AlignmentFile(bam_path, "rb", reference_filename=ref_seq)
    fastafile = pysam.FastaFile(ref_seq)

    while True:
        rec, sample_data, vcf_index = queue.get()
        
        if rec is None:
            bamfile.close()
            fastafile.close()
            break
        chrom = rec['CHROM']
        pos = rec['POS']
        ref = rec['REF']
        alt = rec['ALT']

        variant_id = '{0}:{1}_{2}>{3}'.format(chrom, pos, ref, alt)

        metrics = MetricsDictionary(cohort=cohort, sample=sample, index=vcf_index) 
        try :
            if 'Label' in rec.keys():
                metrics.set_metric('Label', rec['Label'])

            if len(sample_data.get('AD')) <= 1 :
                continue

            dp = sample_data.get('DP')
            metrics.set_metric('tumor_depth', sum(dp) if isinstance(dp, (list, tuple)) else dp)

            if 'AF' in sample_data.keys():
                af = sample_data.get('AF')
                metrics.set_metric('tumor_VAF', af if isinstance(af, float) else af[0])

            for pileupcolumn in bamfile.pileup(contig=chrom, start=pos - 1, stop=pos, min_base_quality=0, min_mapping_quality=0):
                if pileupcolumn.pos == pos - 1:
                    for pileupread in pileupcolumn.pileups:
                        metrics.increment_metric('num_total_reads')
                        query_index = pileupread.query_position
                        read = pileupread.alignment
                        
                        ## This is probably inefficient... there must be ways to do this with pysam without
                        ## the need for extra methods 
                        read_index, distance_5prime, distance_3prime, clipped_length = \
                            process_cigar_tupples(read, pos)

                        if query_index is None or query_index < 0 :
                            # Pretty sure that this is a problem when a read spans an indel
                            continue

                        base = read.query_sequence[query_index]
                        is_ref = base == ref
                        is_var = base == alt

                        
                        if is_read_filtered(read):
                            metrics.increment_metric('tumor_reads_filtered')

                        if is_ref or is_var: 
                            prefix='tumor_ref' if is_ref else 'tumor_var'
                            metrics.increment_metric('tumor_ref_count' if is_ref else 'tumor_var_count')
                            metrics.increment_metric(f'{prefix}_num_minus_strand' if read.is_reverse else f'{prefix}_num_plus_strand')
                            
                            metrics.add_metric(f'{prefix}_base_qualities', read.query_qualities[query_index])
                            metrics.add_metric(f'{prefix}_read_frag_length', read.infer_read_length())
                            metrics.add_metric(f'{prefix}_avg_pos_as_fraction', (read_index / (read.infer_read_length() / 2)))
                            metrics.add_metric(f'{prefix}_distances_to_5p_end', distance_5prime)
                            metrics.add_metric(f'{prefix}_distances_to_3p_end', distance_3prime)
                            metrics.add_metric(f'{prefix}_clipped_length', clipped_length)

                            if read.is_paired:
                                ## this is because we used to track the mapq of unpaired reads seperatly 
                                metrics.add_metric(f'{prefix}_mapping_quality', read.mapping_quality)

                            if read.has_tag('MD'):
                                num_mismatches, mismatch_base_quals = get_mismatch_and_insertion_positions(read)
                                metrics.add_metric('avg_num_mismatches', num_mismatches / read.infer_read_length())
                                if mismatch_base_quals:
                                    metrics.add_metric('avg_sum_mismatch_base_quals',(sum(mismatch_base_quals)))   
                        else:
                            metrics.increment_metric('tumor_other_bases_count')                   
            
            metrics.aggregate_base_metrics(ref, alt) 

            ## Get Window-Based Metrics
            left = max(pos - 500, 0)
            right = min(pos + 500, bamfile.get_reference_length(chrom))

            left_window = [max(0, left - 500), left]
            right_window = [right, min(right + 500, bamfile.get_reference_length(chrom))]
            
            alignment_seq = fastafile.fetch(chrom, left, right)

            metrics.set_metric('window_gc_cont', gc_fraction(alignment_seq))
            metrics.set_metric('window_seq_entropy', entropy(np.unique(list(alignment_seq), return_counts=True)[1] / len(alignment_seq), base=2))
            
            median_cov, cov_variance = get_coverage_in_window(bamfile, fastafile, chrom, left, right)
            if median_cov is None or median_cov == 0:
                metrics.set_metric('window_min_cov_ratio', None)
                metrics.set_metric('window_max_cov_ratio', None)
                continue
            metrics.set_metric('window_median_cov', median_cov)
            metrics.set_metric('window_cov_variance', cov_variance)
            
            left_window_median_cov = get_coverage_in_window(bamfile, fastafile, chrom, *left_window)[0]
            right_window_median_cov = get_coverage_in_window(bamfile, fastafile, chrom, *right_window)[0]

            left_window_median_cov = left_window_median_cov if left_window_median_cov is not None else 0
            right_window_median_cov = right_window_median_cov if right_window_median_cov is not None else 0

            metrics.update_coverage_ratios(left=left_window_median_cov, right=right_window_median_cov)

            fractions = get_read_fractions(bamfile, chrom, left, right)
            metrics.set_metric('window_median_frag_len', fractions[0])
            metrics.set_metric('window_dup_frac', fractions[1])
            metrics.set_metric('window_multi_frac', fractions[2])
            metrics.set_metric('window_improper_frac', fractions[3])
            metrics.set_metric('window_median_mapq', fractions[4])
            metrics.set_metric('window_read_filter_frac', fractions[5])

            ## Just changed to extract tri-/penta-nucleotide sequence (although this assumes that the variant is not
            ## at the last position in the chromosome, is that ok ? )

            sequence = fastafile.fetch(chrom, pos - 3, pos + 1)
            metrics.set_metric('trinucleotide_context', sequence[1:3])
            metrics.set_metric('pentanucleotide_context', sequence)
        
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            iolock.acquire()
            final_dictionary[variant_id] = {}
            iolock.release()
            continue

        iolock.acquire()
        final_dictionary[variant_id] = metrics.get_all_metrics()
        iolock.release()

def read_vcf(sample, label, vcf_path, queue, num_threads):
    vcffile = pysam.VariantFile(vcf_path) 
    for index, rec in enumerate(vcffile.fetch()):
        if rec.ref not in ['A', 'C', 'T', 'G'] or rec.alts[0] not in ['A', 'C', 'T', 'G']:
            continue

        sample_data = dict(rec.samples[sample])

        rec_dict = {'CHROM' : str(rec.chrom),
               'POS' : int(rec.pos),
               'REF' : str(rec.ref),
               'ALT' : str(rec.alts[0])}
        
        if label[0] and label[1]:
            rec_dict['Label'] = 1 if rec.info.get(label[0]) == label[1] else 0
        
        queue.put((rec_dict, sample_data, index)) 

    for i in range(int(num_threads)):
        queue.put((None, None, None))

    vcffile.close()

def get_mobster_tail_scores(sample, vcf_path, out_path, mobster_scores):
    outfile = os.path.join(os.path.dirname(out_path), f'{sample}_mobster.csv')

    script_path = os.path.abspath(__file__)
    directory_name = os.path.dirname(script_path)

    process = subprocess.run(
        # ['Rscript', directory_name + '/run_mobster.R', sample, vcf_path, outfile],
        [directory_name + '/run_mobster.R', sample, vcf_path, outfile],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
        )
    
    logger.info(process.stdout)
    logger.error(process.stderr)
    
    if not os.path.isfile(outfile):
        logger.info(f"MOBSTER did not run succesfully on {sample}")
    else:    
        with open(outfile, newline='') as mfile:
            logger.info(f"Finished MOBSTER calculations for {sample}")
            reader = csv.reader(mfile, delimiter=',')
            header = next(reader)
            for row in reader: 
                sample, chrom, pos, REF, ALT, Tail = row
                if REF not in ['A', 'C', 'T', 'G'] or ALT not in ['A', 'C', 'T', 'G']:
                    continue
                variant_id = '{0}:{1}_{2}>{3}'.format(chrom, pos, REF, ALT)
                mobster_scores[variant_id] = {'Tail': Tail}
        mfile.close()
        os.remove(outfile)

def extract_all_features(bam_path, vcf_path, ref_seq, sample, cohort, label, num_threads, output_file):
    is_compressed = vcf_path.endswith('.gz')
    has_index = os.path.exists(vcf_path + '.csi') or os.path.exists(vcf_path + '.tbi')
    if not is_compressed:
        logger.error("VCF File is not compressed. Please bgzip the VCF file before running feature extraction.")
        exit(1)
    if not has_index:
        logger.error("VCF File is not indexed. Please index the VCF file before running feature extraction.")
        exit(1) 

    start = time.time()
    
    queue = Queue(maxsize=500) 

    iolock = Lock()
    to_df = []
    num_vars = 0

    with Manager() as manager:
        mobster_scores = manager.dict()
        all_features = manager.dict()

        pool = [Process(target=process_variant, args=(queue, sample, cohort, bam_path, ref_seq, iolock, all_features)) for i in range(int(num_threads))]
        pool.insert(0, Process(target=get_mobster_tail_scores, args=(sample, vcf_path, output_file, mobster_scores)))
        for P in pool:
            P.start()

        read_vcf(sample, label, vcf_path, queue, num_threads)

        for P in pool:
            P.join()
        
        ## Takes care of cases when MOBSTER doesn't run succesfully
        result = {variant: {**all_features.get(variant, {}),**(mobster_scores.get(variant, {'Tail': 1}))}
            for variant in all_features.keys()}

        num_vars = len(all_features)
        to_df = [{'Variant': variant, **metric} for variant, metric in result.items()]
    
    pd.DataFrame(to_df).sort_values("vcf_index").drop("vcf_index", axis=1).to_csv(output_file, index=False)

    end = time.time() - start
    logger.info(f"Finished all metrics for {num_vars} vars in {sample} in {end}")
    logger.info(f"Features stored: \n{output_file}")

    ## Make a bunch of runs with different numbers of CPUs and see how it scales
    ## ask Jen about nice seff command
    ## Add num threads to the options 

#################################################################################
################################################ SPLIT VARIANTS for PROCESSING ##

#################################################################################
### PROCESS INPUT FILES #########################################################
    
def process_sample(sample, cohort, vcf_path, bam_path, ref_seq, output_file, label, num_threads): 
    try:
        if os.path.isfile(vcf_path) and os.path.isfile(bam_path) and os.path.isfile(ref_seq):
            logger.info(f'Processing BAM file for sample: {sample}')
            extract_all_features(bam_path, vcf_path, ref_seq, sample, cohort, label, num_threads, output_file)
                
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
bam_file, ref_seq, cohort=None):
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
        
        process_sample(sample, cohort, vcf_file, bam_file, ref_seq, output_file, label, num_threads)
    except Exception as e:
        logger.error(f"An error occurred: {e}")

#################################################################################
######################################################### PROCESS INPUT FILES ###