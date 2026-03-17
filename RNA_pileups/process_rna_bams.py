#!/bin/env python

##########################################################################
#	USAGE: extract features from BAM file
#   DESCRIPTION: Functions for extracting features for variants in VCF file
#   from user provided BAM file
#   Created by Valentina Grether
##########################################################################

import pysam  # type: ignore
import pysamstats # type: ignore
import numpy as np # type: ignore
import pandas as pd # type: ignore
import pickle
import csv
import os, sys
import re
from statistics import mean, median
from sklearn.preprocessing import OneHotEncoder # type: ignore
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
from threading import Lock
import traceback
import subprocess
import multiprocessing as mp
from multiprocessing import Pool
import time
from collections import OrderedDict

lock = Lock()

def pileup(bamfile, rec, ref_seq, FFPE, label):
    if len(dict(rec.samples)[FFPE].get('AD')) <= 1 :
        print(f'SKIPPING {position} because len(AD) < 1')
        return {}
    
    metrics={}
    position = '{0}:{1}_{2}>{3}'.format(rec.chrom, rec.pos, rec.ref, rec.alts[0])

    if label[0] and label[1]:
        metrics[position]= {
            'EBM_LABEL' : 1 if int(rec.info.get(label[0])) == int(label[1]) else 0
        }
    else:
        metrics[position]= {}
                
    metrics[position].update({
        'tumor_ref_count': 0,
        'tumor_var_count': 0,
        'tumor_other_bases_count': 0
    })
        
    for pileupcolumn in bamfile.pileup(contig=rec.chrom, start=rec.pos - 1, stop=rec.pos, min_base_quality=0, min_mapping_quality=0):
        if pileupcolumn.reference_pos == rec.pos - 1:
            metrics[position].update({
                'total_coverage': pileupcolumn.nsegments #num_total_reads
            })

            if pileupcolumn.nsegments is None:
                print(f'{position} has no coverage')
                continue

            for pileupread in pileupcolumn.pileups:
                query_index = pileupread.query_position
                 
                if query_index is None or query_index < 0 :
                    print(f'{FFPE} {position}: Query index is None or < 0')
                    continue

                base = pileupread.alignment.query_sequence[query_index]
                if base == rec.ref:
                    metrics[position]['tumor_ref_count'] += 1
                elif base == rec.alts[0]:
                    metrics[position]['tumor_var_count'] += 1
                else:
                    metrics[position]['tumor_other_bases_count'] += 1

    return metrics

def get_chromosome_metrics(chrom, vcf_path, bam_path, ref_seq, sample, label):
    pileup_counts = OrderedDict()
    bamfile = pysam.AlignmentFile(bam_path, "rb", reference_filename=ref_seq)
    fastafile = pysam.FastaFile(ref_seq)
    vcffile = pysam.VariantFile(vcf_path) 

    for rec in vcffile.fetch(contig=chrom):
        if rec.chrom != chrom:
            continue
        position = '{0}:{1}_{2}>{3}'.format(rec.chrom, rec.pos, rec.ref, rec.alts[0])
        if rec.ref not in ['A', 'C', 'T', 'G'] or rec.alts[0] not in ['A', 'C', 'T', 'G']:
            continue
        pileup_counts.update(pileup(bamfile, rec, ref_seq, sample, label))

    vcffile.close()
    bamfile.close()
    fastafile.close()
    return pileup_counts

def get_counts(bam_path, vcf_path, ref_seq, FFPE, label, num_processes=4):
    print("Starting Pileup Counts")
    start = time.time()
    async_results = []
    print("Compressing VCF")
    os.makedirs('temp', exist_ok=True)
    if vcf_path.endswith('.gz'):
        compressed_vcf = os.path.join(os.getcwd(), 'temp/', os.path.basename(vcf_path))
        print(compressed_vcf)
        os.system(f'cp {vcf_path} {compressed_vcf}')  
    else:
        compressed_vcf = os.path.join(os.getcwd(), 'temp/', os.path.basename(vcf_path) +'.gz')
        print(compressed_vcf)
        subprocess.run(f"bcftools view -O z -o {compressed_vcf} {vcf_path}", shell=True)
    print("Compressed!") 

    index_file = f"{compressed_vcf}.csi" 
    if not os.path.exists(index_file):
        remove_index = True
        subprocess.run(f"bcftools index {compressed_vcf}", shell=True)
        print("Indexed!")
    else:
        remove_index = False
        print("Index file already exists")

    vcffile = pysam.VariantFile(compressed_vcf) 

    pool = Pool(processes=num_processes)

    for contig in vcffile.header.contigs:
        async_result = pool.apply_async(get_chromosome_metrics, (contig, compressed_vcf, bam_path, ref_seq, FFPE, label))
        async_results.append(async_result)
    vcffile.close()
    pool.close()
    pool.join()

    pileup_counts = {}
    for async_result in async_results:
        pileup_counts.update(async_result.get())
    
    print("Got all results")
    print("Now removing compressed VCF")
    os.remove(compressed_vcf)
    if remove_index:
        print("and removing the index file we made lol")
        os.remove(index_file)
    
    end = time.time() - start
    num_vars = len(pileup_counts.keys())
    print(f"Finished metrics (in parallel) for {num_vars} in {FFPE} in {end}")
    return pileup_counts

def get_column_idx(header):
    # Find the index of the columns matching the patterns and store the results in a dictionary
    patterns = {
        'SAMPLE': re.compile(r'\s*(FFPE|Sample)\s*', re.IGNORECASE),
        'VCF_FILE': re.compile(r'VCF[_\s]*FILE', re.IGNORECASE),
        'BAM_FILE': re.compile(r'BAM[_\s]*FILE', re.IGNORECASE),
        'REF_SEQ': re.compile(r'REF[_\s]*SEQ', re.IGNORECASE)
    }
    columns = {}
    for keyword, pattern in patterns.items():
        columns[keyword] = next((i for i, col in enumerate(header) if pattern.match(col)), -1)
    
    matched_columns = [
        columns['SAMPLE'],
        columns['VCF_FILE'],
        columns['BAM_FILE'],
        columns['REF_SEQ']
    ]
    if -1 in matched_columns:
        print("Error: One or more required columns not found in pairs file.\nRequired columns are: \nFFPE, VCF_FILE, BAM_FILE, REF_SEQ")
        exit(1)
    return matched_columns
    
def process_row(sample, vcf_path, bam_path, ref_seq, outpath, label):   
    try:
        if os.path.isfile(vcf_path) and os.path.isfile(bam_path) and os.path.isfile(ref_seq):
            print(f'Processing BAM file for sample: {sample}')
            flattened_data = []
            metrics = get_counts(bam_path, vcf_path, ref_seq, sample, label)

            for variant, metric in metrics.items(): 
                row = {'Sample': sample, 'Variant': variant}
                row.update(metric)
                flattened_data.append(row)

            df = pd.DataFrame(flattened_data)
            with lock:
                output_file = os.path.join(outpath, f"{sample}_variant_counts.csv")
                df.to_csv(output_file, index=False)
                print(f"Features extracted and stored: \n{output_file}")

    except Exception as e:
        print(f"Error processing {sample}: {e}")
        traceback.print_exc()

def process_pairs_path(pairs_path, outpath, label):
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
                        outpath, 
                        label) for row in reader]

            for future in futures:
                future.result()
    pairsfile.close()

def process_bam_file(outpath, label, sample=None, vcf_file=None, 
bam_file=None, ref_seq=None, pairs_path=None):
    try:
        if os.path.isfile(outpath):
            outpath = os.path.basename(outpath)

        if not os.path.exists(outpath):
            os.makedirs(outpath)
            
        if pairs_path is not None:
            process_pairs_path(pairs_path, outpath, label)
        else:
            process_row(sample, vcf_file, bam_file, ref_seq, outpath, label)
    except Exception as e:
        print(f"An error occurred: {e}")
    pass

