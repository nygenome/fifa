#!/bin/env python

##########################################################################
#	USAGE: extract features from BAM file
#   DESCRIPTION: Functions for extracting features for variants in VCF file
#   from user provided BAM file
#   Created by Valentina Grether
##########################################################################

import pysam 
import pysamstats
import numpy as np # type: ignore
import pandas as pd # type: ignore
import pickle
import csv
import os, sys
import re
from statistics import mean, median
from sklearn.preprocessing import OneHotEncoder # type: ignore
from scipy.stats import entropy # type: ignore
from Bio.SeqUtils import gc_fraction # type: ignore
from Bio import SeqIO # type: ignore
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
from threading import Lock
import traceback
import subprocess
import multiprocessing as mp
from multiprocessing import Pool
import time

global contig_sequences
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

def get_ref_seq(fasta_file):
    contig_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_sequences[record.id] = record.seq
    return contig_sequences


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


def update_metrics(metrics, position, metrics_data, prefix):
    metrics[position].update({
        f'{prefix}_med_frag_len': safe_median(metrics_data.get('read_frag_length', []))
    })

def get_base_metrics(bamfile, rec, ref_seq, FFPE, label):
    if len(dict(rec.samples)[FFPE].get('AD')) <= 1 :
        print(f'SKIPPING {position} because len(AD) < 1')
        return {}

    metrics={}

    position = '{0}:{1}_{2}>{3}'.format(rec.chrom, rec.pos, rec.ref, rec.alts[0])
    metrics[position]= {}

    metrics[position].update({
        'tumor_ref_count': 0,
        'tumor_var_count': 0 })
    
    ref_metrics = {
        'read_frag_length': []
    }
        
    var_metrics = {
        'read_frag_length': []
    }
    
    for pileupcolumn in bamfile.pileup(contig=rec.chrom, start=rec.pos - 1, stop=rec.pos,min_base_quality=0, min_mapping_quality=0):
        if pileupcolumn.pos == rec.pos - 1:
            for pileupread in pileupcolumn.pileups:
                query_index = pileupread.query_position
                read = pileupread.alignment
            
                if len(rec.alts[0]) <= 0:  
                    print(f'{FFPE} sample has no alternative alleles at: {position}')
                    continue
                
                if query_index is None or query_index < 0 :
                    print(f'{FFPE} {position}: Query index is None or < 0')
                    continue

                base = read.query_sequence[query_index]
                is_ref = base == rec.ref
                is_var = base == rec.alts[0]

                if is_ref or is_var:
                    metrics[position]['tumor_ref_count' if is_ref else 'tumor_var_count'] += 1
                    metrics_list = ref_metrics if is_ref else var_metrics 
                    metrics_list['read_frag_length'].append(abs(read.template_length))
                    if rec.chrom == 'chr1': 
                        print(abs(read.template_length))
    
    if metrics[position]['tumor_ref_count'] > 0:
        update_metrics(metrics, position, ref_metrics, 'tumor_ref')
    else:
        metrics[position]['tumor_ref_med_frag_len'] = 0
    
    if metrics[position]['tumor_var_count'] > 0:
        update_metrics(metrics, position, var_metrics, 'tumor_var')
    else:
        metrics[position]['tumor_var_med_frag_len'] = 0
    
    metrics[position].pop("tumor_ref_count")
    metrics[position].pop("tumor_var_count")

    return metrics
    

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
    frag_lenths = []    
    for read in bamfile.fetch(chrom, left, right):
        filtered_SAM_flags, filtered_mutect2 = is_read_filtered(read)
        if read.query_alignment_end > bamfile.get_reference_length(chrom) or filtered_mutect2:
            pass
        else:
            frag_lenths.append(abs(read.template_length))
    
    return (
        safe_median(frag_lenths)
    )

def metrics_500_bp_window(rec, bamfile): 
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

    position_metrics['window_median_frag_len'] = get_read_fractions(bamfile, chrom, left, right)
    
    return position_metrics

def get_chromosome_metrics(chrom, vcf_path, bam_path, ref_seq, sample, label):
    window_metrics = {}
    base_metrics = {}
    bamfile = pysam.AlignmentFile(bam_path, "rb", reference_filename=ref_seq)
    fastafile = pysam.FastaFile(ref_seq)
    vcffile = pysam.VariantFile(vcf_path) 

    for rec in vcffile.fetch(contig=chrom):
        if rec.chrom != chrom:
            continue
        position = '{0}:{1}_{2}>{3}'.format(rec.chrom, rec.pos, rec.ref, rec.alts[0])
        if rec.ref not in ['A', 'C', 'T', 'G'] or rec.alts[0] not in ['A', 'C', 'T', 'G']:
            continue
        window_metrics[position] = metrics_500_bp_window(rec, bamfile)
        base_metrics.update(get_base_metrics(bamfile, rec, ref_seq, sample, label))

    vcffile.close()
    bamfile.close()
    fastafile.close()
    return base_metrics, window_metrics

def get_features_in_parallel(bam_path, vcf_path, ref_seq, FFPE, label, num_processes=4):
    print("Starting window metrics")
    start = time.time()
    async_results = []
    print("compressing VCF")
    os.makedirs('/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/temp', exist_ok=True)
    if vcf_path.endswith('.gz'):
        compressed_vcf = os.path.join('/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/temp/', os.path.basename(vcf_path))
        print(compressed_vcf)
        os.system(f'cp {vcf_path} {compressed_vcf}')  
    else:
        compressed_vcf = os.path.join('/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/temp/', os.path.basename(vcf_path) +'.gz')
        print(compressed_vcf)
        subprocess.run(f"bcftools view -O z -o {compressed_vcf} {vcf_path}", shell=True)
    print("compressed") 

    index_file = f"{compressed_vcf}.csi" 
    subprocess.run(f"bcftools index -f {compressed_vcf}", shell=True)
    print("indexed")

    vcffile = pysam.VariantFile(compressed_vcf) 

    pool = Pool(processes=num_processes)

    for contig in vcffile.header.contigs:
        async_result = pool.apply_async(get_chromosome_metrics, (contig, compressed_vcf, bam_path, ref_seq, FFPE, label))
        async_results.append(async_result)
    vcffile.close()
    pool.close()
    pool.join()

    window_metrics = {}
    base_metrics = {}
    for async_result in async_results:
        base, window = async_result.get()
        base_metrics.update(base)
        window_metrics.update(window)
    
    print("got all results")
    print("now removing compressed vcf")
    os.remove(compressed_vcf)
    print("and removing the index file we made lol")
    os.remove(index_file)
    
    end = time.time() - start
    num_vars = len(window_metrics.keys())
    print(f"Finished metrics (in parallel) for {num_vars} in {FFPE} in {end}")
    return base_metrics, window_metrics

def get_metrics_in_parallel(bam_path, vcf_path, ref_seq, sample, label):    
    global contig_sequences 
    contig_sequences = get_ref_seq(ref_seq)

    with ThreadPoolExecutor(max_workers=4) as executor:
        future_window_base_metrics = executor.submit(get_features_in_parallel, bam_path, vcf_path, ref_seq, sample, label)
    
    base_metrics, window_metrics = future_window_base_metrics.result()

    return  base_metrics, window_metrics

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
        print("Error: One or more required columns not found in pairs file.\nRequired columns are: \nFFPE, Cohort, VCF_FILE, BAM_FILE, REF_SEQ")
        exit(1)
    return matched_columns
    
def process_row(sample, cohort, vcf_path, bam_path, ref_seq, output_file, label):   
    try:
        if os.path.isfile(vcf_path) and os.path.isfile(bam_path) and os.path.isfile(ref_seq):
            print(f'Processing BAM file for sample: {sample}')
            flattened_data = []
            base_metrics, window_metrics = get_metrics_in_parallel(bam_path, vcf_path, ref_seq, sample, label)
            print(len(base_metrics))
            print(len(window_metrics))

            metrics = {
                variant: {**base_metrics.get(variant, {}), **window_metrics.get(variant, {})}
                for variant in set(base_metrics) & set(window_metrics)
            }

            for variant, metric in metrics.items(): 
                row = {'Cohort': cohort, 'Sample': sample, 'Variant': variant}
                row.update(metric)
                flattened_data.append(row)

            df = pd.DataFrame(flattened_data)
            with lock:
                df.to_csv(output_file, index=False)
                print(f"Features extracted and stored: \n{output_file}")

        else:
            print("VCF: ")
            print(os.path.isfile(vcf_path))
            print("BAM : ")
            print(os.path.isfile(bam_path))
            print("REF : ")
            print(os.path.isfile(ref_seq))

    except Exception as e:
        print(f"Error processing {sample}: {e}")
        traceback.print_exc()

def process_pairs_path(pairs_path, output_file, label):
    """
    Main function for processing BAM file for variants in VCF file. 
    Pairsfile must include required columns 

    Args:
        pairs_path (str): path to input pairs file
        outpath (str): path for output csv file with features for each variant across samples
        max_workers (int): maximum number of threads. Default = 4
    
    """   
    with open(pairs_path, newline='') as pairsfile:
        reader = csv.reader(pairsfile, delimiter=',')
        
        matched_columns = get_column_idx(next(reader))

        with ProcessPoolExecutor() as executor:  
            futures = [
                    executor.submit(
                        process_row, 
                        row[matched_columns[0]], 
                        row[matched_columns[1]], 
                        row[matched_columns[2]], 
                        row[matched_columns[3]], 
                        row[matched_columns[4]], 
                        output_file, 
                        label) for row in reader]

            for future in futures:
                future.result()
    pairsfile.close()

def process_bam_file(outpath, label, sample=None, cohort=None, vcf_file=None, 
bam_file=None, ref_seq=None, pairs_path=None):
    try:
        if os.path.isfile(outpath):
            output_file = outpath 
        else:
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            output_file=os.path.join(outpath, f"{sample}_extracted_features.csv")
            
        if pairs_path is not None:
            process_pairs_path(pairs_path, output_file, label)
        else:
            process_row(sample, cohort, vcf_file, bam_file, ref_seq, output_file, label)
    except Exception as e:
        print(f"An error occurred: {e}")
    pass

