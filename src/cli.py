#!/usr/bin/env python

import os 
import logging
import argparse as ap
from process_bam_file import process_bam_file as regular_process_bam_file
from parallelizing_bam_metrics import process_bam_file as parallel_process_bam_file
import train_new_ebm
import train_with_hyperparameter
from merge_models import merge_ebms
from classify_with_scaling import predict
from recover_annotations import predict as predict_with_rna

VERSION = 0.1

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    parser = ap.ArgumentParser(description='FIFA for Filtering FFPE Artifacts\n'
                                           'Version {0}\n'
                                           'Developed by Valentina Grether'.format(VERSION),
                               formatter_class=ap.RawTextHelpFormatter)

    subparsers = parser.add_subparsers(dest='subcommand', required=True)

    ## Extract features from BAM/VCF files for making predictions or retraining
    feature_parser = subparsers.add_parser('extract', help='Extract features from VCF/BAM files for making predictions or retraining')
    feature_parser.add_argument('-s', '--sample', required=True, 
                            help='Name of sample being submitted to be processed.')
    feature_parser.add_argument('-v', '--vcffile', required=True, 
                            help='Path to sample\'s VCF file, with variants to extract features for.')
    feature_parser.add_argument('-b', '--bamfile', required=True, 
                            help='Path to sample\'s BAM file, with reads to extract features for.')
    feature_parser.add_argument('-r', '--refseq', required=True, 
                            help='Path to reference file used for alignment of the sample\'s BAM file.')
    feature_parser.add_argument('-c', '--cohort', default='NA', 
                            help='Sample\'s Cohort')
    feature_parser.add_argument('-n', '--num_threads', default=2, type=int,
                                help='Number of threads allocated')            
    feature_parser.add_argument('-o', '--output_path', default=os.getcwd(),
                                help='Path to save output files with variants and all extracted features. (default: current directory)')
    feature_parser.add_argument('-l', '--label', nargs=2, default=[None,None],
                                help='VCF INFO Field with Variant Label and Real/Truth/1 Label. (default: None, None)')
    
    ## Argument to choose between my original paralelization scheme or new one after talking with Andre
    ## helpul for testing resource allocation 
    feature_parser.add_argument('-p ', '--original_parallel', action='store_true',
                                help='Boolean: Whether we\'re using the original (true) or new script (false)')
    
    ## Train new EBM model on user-inputted cohort
    retrain_parser = subparsers.add_parser('retrain', help='Train new FIFA (EBM) model with user inputted data')
    retrain_parser.add_argument('-o', '--output_path', default=os.path.join(os.getcwd(), 'fifa_model.pkl'), 
                                help='Path to save new model. (default: "fifa_model.pkl" in current directory)')
    retrain_parser.add_argument('-d', '--directory', help='Path to directory containing all samples for training a new cohort\'s model.' \
    ' Each sample file should be a CSV with extracted features for all variants in the sample.')
    retrain_parser.add_argument('-l', '--labels_path', default=None, 
                                help='Path to CSV file with true labels for all variants (True Labels should be: Real or 1)')
    retrain_parser.add_argument('-hp ', '--hyperparameter', action='store_true',
                                help='Boolean: Whether to conduct hyperparameter grid-search when training new EBM model')
    
    ## Merge multiple EBM models together
    merge_parser = subparsers.add_parser('merge', help='Merge together multiple EBM models and test performance on HCC1395 chr1 variants')
    merge_parser.add_argument('-o', '--output_path', default=None, 
                                help='Path to save merged FIFA model. If not provided, the merged model is not saved.', required=False)
    merge_parser.add_argument('-m', '--input_models', nargs='+', 
                                help='Path file(s) of EBM models to be merged. '
                                'If intending to use original FIFA models, please specify: NYGC1, NYGC2, HTMCP, or BLGSP', required=True)
    
    ## Make predictions on new samples
    predict_parser = subparsers.add_parser('predict', help='Create predictions on labeled or unlabeled data')
    predict_parser.add_argument('-s', '--sample', required=True,
                            help='Sample Name \nFor submitting individual sample for processing.')
    predict_parser.add_argument('-v', '--vcffile', required=True,
                            help='Path to sample\'s VCF file, with variants to classify.')
    predict_parser.add_argument('-f', '--features_path', nargs='+', 
                                    help='Paths to CSV files with features for each variant', required=True)
    predict_parser.add_argument('-o', '--output_dir', default=os.getcwd(), 
                                help='Directory to save VCFs with model predictions. (default: current directory)')
    predict_parser.add_argument('-m', '--model_files', nargs='+', 
                                help='Path file(s) of EBM models to be used. '
                                '(If submitting multiple models, they will be merged.)', required=True)
    predict_parser.add_argument('-r', '--rna_annotations', default=None, 
                                help='Path to file with RNA annotations', required=False)
    
    args = parser.parse_args()

    if args.subcommand == 'extract':
        try:
            logger.info('Extracting features from BAM/VCF files')
            if args.original_parallel:
                logger.info('Using original parallelization scheme. Less efficient than latest version (omit -p flag)')
                regular_process_bam_file(outpath=args.output_path, label=args.label, num_threads=args.num_threads, sample=args.sample, cohort=args.cohort, 
                vcf_file=args.vcffile, bam_file=args.bamfile, ref_seq=args.refseq)
            else:
                parallel_process_bam_file(outpath=args.output_path, label=args.label, num_threads=args.num_threads, sample=args.sample, cohort=args.cohort, 
                vcf_file=args.vcffile, bam_file=args.bamfile, ref_seq=args.refseq)
        except Exception as e:
            logger.error(f"An error occurred: {e}")

    elif args.subcommand == 'retrain':
        logger.info('Training new EBM model with user-inputted cohort')
        if args.hyperparameter:
            logger.info('Conducting 5-fold hyperparameter grid-search')
            train_with_hyperparameter.retrain(args.directory, args.labels_path, args.output_path)
        else:
            train_new_ebm.retrain(args.directory, args.labels_path, args.output_path)
        
    elif args.subcommand == 'predict':
        logger.info('Classify variants in a VCF')
        try:
            if isinstance(args.model_files, list):
                model_file = list(args.model_files)
            else:
                model_file = args.model_files
            if(args.rna_annotations):
                logger.info('Will include RNA annotations (for rescuing FN predictions)')
                predict_with_rna(extracted_features_paths=list(args.features_path), outpath=args.output_dir, model_file=model_file, 
                sample=args.sample, vcf_path=args.vcffile, rna_path=args.rna_annotations)   
            else:
                predict(extracted_features_paths=list(args.features_path), outpath=args.output_dir, model_file=model_file, 
                sample=args.sample, vcf_path=args.vcffile)   
        except Exception as e:
            logger.error(f"An error occurred: {e}")
    
    elif args.subcommand == 'merge':
        logger.info('Merging multiple EBM models')
        merge_ebms(args.input_models, args.output_path)
        pass
    else:
        logger.error('Please choose proper function.')
