#!/bin/env python

##########################################################################
#	USAGE: classify variants in a test set
#
#   DESCRIPTION: Given a file with extracted features and 1+ EBM models, 
#   make predictions on the variants and write out VCF file with EBM annotations
#
#   MERGING EBM MODELS: Import EBM model(s) from pickle file(s). If multiple model paths are provided,
#   the models are merged together. To make models compatible for merging, features 
#   with zero importance are removed, and features are aligned so that all models have the 
#   same feature set. 
#   
#   GENERATING PREDICTIONS: Predictions are added as an annotation in the INFO/FIFA_LABEL field
#   of the output VCF file, and probability of being a "real" variant added to INFO/FIFA_PROB field. 
#
#   Created by Valentina Grether
##########################################################################

import os
import pickle
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
from scipy.stats.mstats import winsorize
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, balanced_accuracy_score, f1_score, recall_score
from interpret.glassbox import ExplainableBoostingClassifier
from interpret.glassbox import merge_ebms
import multiprocessing as mp
from multiprocessing import Pool
from json import dump
from scipy.stats import trim_mean
from sklearn.preprocessing import FunctionTransformer, StandardScaler
import sys
import glob
import re
import pysam
import csv
from process_bam_file import get_column_idx
from helper_funcs import convert_hot_encodings, scale_features
import logging

global logger
global predictions_path

np.random.seed(42)

def make_predictions(ebm, variants):
    """
    Scale features, and make predictions on the variants. 

    Args:
        ebm (ExplainableBoostingClassifier): ebm model 
        variants (pd.DataFrame): data frame with variants and all features

    """   
    label_encoder = LabelEncoder() 

    if 'Label' in variants.columns:
        y_test = label_encoder.fit_transform(variants['Label'].apply(lambda x: 1 if x == 'Real' else 0))
        variants.drop(columns='Label', inplace=True)

    ## Small check to make sure that I converted the hot-encondings we were originally using 
    ## to the sequence context features 
    if 'hot_encoded_ref_base' in variants.columns:
        variants = convert_hot_encodings(variants)
    
    # Transform features 
    variants = scale_features(variants)
    
    variant_ids_test = variants[['Sample', 'Variant', 'Cohort']].reset_index(drop=True)

    X_test = variants.drop(['Sample', 'Variant', 'Cohort'], axis=1)
    del variants 

    if not (set(ebm.feature_names_in_).issubset(set(X_test.columns))):
        print("Warning: Feature mismatch between model and test data. Please check input features.")
        sys.exit(1)

    y_pred = pd.DataFrame(ebm.predict(X_test), columns=['Predicted'])
    y_probs = pd.DataFrame(ebm.predict_proba(X_test)[:, 1], columns=['Probability'])
    predictions = pd.concat([variant_ids_test, y_pred, y_probs], axis=1)
    return predictions

def generate_output_vcf_file(predictions, vcf_path): 
    """
    Generates output VCF file with FIFA predictions added as an annotation in the INFO field,
    And probability of being a "real" variant added as an independent annotation. 

    Args:
        predictions (pd.DataFrame): data frame with variants and model predictions
        vcf_path (str): path to the original VCF file
    """   

    global predictions_path 
    outpath = os.path.join(predictions_path, re.sub(r'\.vcf(\.gz)?$', '.fifa.vcf', os.path.basename(vcf_path)))
    
    if os.path.isfile(outpath):
        os.remove(outpath)

    try:
        origvcf = pysam.VariantFile(vcf_path, "r")
        origvcf.header.info.add("FIFA_LABEL", "1", "Integer", "Binary FIFA Classification (0-Artifact or 1-Real)")
        origvcf.header.info.add("FIFA_PROB", "1", "Float", "FIFA - Probability of Variant being Real, based on FIFA's Classification")
        vcf_out = pysam.VariantFile(outpath, 'w', header=origvcf.header)
        
        prediction_dict = dict(zip(predictions['Variant'], predictions['Predicted']))
        probability_dict = dict(zip(predictions['Variant'], predictions['Probability']))
        positions = set()
    
        for variant in origvcf:
            position = f"{variant.chrom}:{variant.pos}_{variant.ref}>{variant.alts[0]}"
            
            if position in positions:
                print(f"Duplicate entry found for: {position}")
                continue

            pred = prediction_dict.get(position)
            prob = probability_dict.get(position)
            
            if (pred is not None) and (prob is not None) and variant.info.get('HighConfidence'):
                positions.add(position)
                variant.info['FIFA_LABEL'] = int(pred)
                variant.info['FIFA_PROB'] = round(prob, 4)
            
            vcf_out.write(variant)

        logger.info(f"Successfully generated: {outpath}")

    except Exception as e:
        print(f"An error occurred: {e}")

def parse_pairs_path(pairs_path, combined_df): 
    with open(pairs_path, newline='') as pairsfile:
        reader = csv.reader(pairsfile, delimiter=',')
        matched_columns = get_column_idx(next(reader))
        
        for row in reader: 
            sample = row[matched_columns[0]]
            vcf_path = row[matched_columns[2]]
            predictions = combined_df[combined_df['Sample'] == sample]
            generate_output_vcf_file(predictions, vcf_path)

def get_zero_score_names(models):
    # Accesory function to identify non-contributory features across multiple EBM models
    # These features need to be removed to make models compatible for merging 
    zero_score_names = set()
    for model in models:
        zero_score_names.update(
            name for name, scores in zip(model.term_names_, model.term_scores_)
            if np.sum(np.abs(scores)) == 0
        )
    return zero_score_names

def remove_zero_score_terms_by_name(model, zero_score_names):
    return model.remove_terms(zero_score_names)

def predict(extracted_features_paths, outpath, model_file, pairs_path=None, sample=None, vcf_path=None):
    global predictions_path 
    global logger 

    logger = logging.getLogger(__name__)
    try:
        if os.path.isdir(outpath):
            predictions_path = outpath
        elif os.path.isfile(outpath):
            predictions_path = os.path.basename(outpath)
        else: # not os.path.exists(predictions_path):
            os.makedirs(outpath, exist_ok=True)
            predictions_path = outpath
    except Exception as e:
        print(f"An error occurred: {e}")
    
    try: 
        variants = pd.concat(
            (pd.read_csv(path).fillna(0) for path in extracted_features_paths), ignore_index=True)
    except Exception as e:
        print(f"An error occurred: {e}")

    try:
        # Merge cohort-level EBM models together 
        # Remove non-contributory features and align features in the same order
        # across models to make them compatible for merging
        if isinstance(model_file, list):
            models = []
            for file in model_file:
                ebm = ExplainableBoostingClassifier()
                with open(file, 'rb') as f:
                    ebm = pickle.load(f)
                models.append(ebm)

            zero_score_names = get_zero_score_names(models)
            for i, model in enumerate(models):
                    models[i] = remove_zero_score_terms_by_name(model, zero_score_names)
            ebm1 = models[0]
            for model in models[1:]:
                model_feature_index = {name: i for i, name in enumerate(model.feature_names_in_)}
                reorder_indices = [model_feature_index[name] for name in ebm1.feature_names_in_]
                
                model.feature_names_in_ = [model.feature_names_in_[i] for i in reorder_indices]
                model.term_names_ = [model.term_names_[i] for i in reorder_indices]
                model.feature_types_in_ = [model.feature_types_in_[i] for i in reorder_indices]
            
            ebm =  merge_ebms(models)

        # Or load single EBM model
        elif os.path.isfile(model_file):
            ebm = ExplainableBoostingClassifier()
            with open(model_file, 'rb') as f:
                ebm = pickle.load(f)
            
    except Exception as e:
        print(f"An error occurred: {e}")
        print(f"The path for your ebm model is incorrect. \n{model_file}")
    
    combined_df = make_predictions(ebm, variants)

    if(pairs_path):
        if not os.path.isfile(pairs_path):
            print(f"The path for your sample files is incorrect")
            exit(1)
        parse_pairs_path(pairs_path, combined_df) 
    else: 
        generate_output_vcf_file(combined_df, vcf_path)
