import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from interpret.glassbox import ExplainableBoostingClassifier
from json import dump
from sklearn.preprocessing import FunctionTransformer, StandardScaler
import sys
import glob

## A random assortment of helper functions
## So that I'm not repeating code

#Convert hot encodings to seq context features
def convert_hot_encodings(variants):
    one_hot_to_base = {
        '[1. 0. 0. 0. 0.]': 'A', #A
        '[0. 1. 0. 0. 0.]': 'C', #C
        '[0. 0. 1. 0. 0.]': 'G', #G 
        '[0. 0. 0. 1. 0.]': 'T', #T
        '[0. 0. 0. 0. 1.]': 'N', #N
        '[0. 0. 0. 0. 0.]': 'N'
    }

    categorical_columns = ['left_two_base', 'left_one_base', 'hot_encoded_ref_base', 'right_one_base', 'right_two_base', 'hot_encoded_var_base']
    
    if set(categorical_columns).issubset(variants.columns): 
        for col in categorical_columns:
            variants[col] = variants[col].apply(lambda x: one_hot_to_base.get(str(x), x))

        variants = variants[~variants[categorical_columns].apply(lambda x: x.str.contains('N')).any(axis=1)]
        variants['trinucleotide_context'] = variants['left_one_base'] + '[' + variants['hot_encoded_ref_base'] + '>' + variants['hot_encoded_var_base'] + ']' + variants['right_one_base']
        variants['pentanucleotide_context'] = variants['left_two_base'] + variants['left_one_base'] + '[' + variants['hot_encoded_ref_base'] + '>' + variants['hot_encoded_var_base'] + ']' + variants['right_one_base'] + variants['right_two_base']

        variants = variants.drop(categorical_columns, axis=1)

        if 'tumor_ref_avg_se_mapping_quality' in variants.columns and 'tumor_var_avg_se_mapping_quality' in variants.columns:
            variants = variants.drop(['tumor_ref_avg_se_mapping_quality',
                                'tumor_var_avg_se_mapping_quality'], axis=1)

    categorical_columns = ['trinucleotide_context', 'pentanucleotide_context']

    for col in categorical_columns:
        variants[col] = variants[col].astype('category')
        variants = variants[variants[col].notnull()]
    return variants

def scale_features(variants):
    features_to_scale = [
    "tumor_reads_filtered",
    "tumor_ref_avg_base_quality",
    "tumor_ref_avg_clipped_length",
    "tumor_ref_avg_num_mismatches_as_fraction",
    "tumor_ref_avg_sum_mismatch_qualities",
    "tumor_ref_med_base_quality",
    "tumor_ref_med_frag_len",
    "tumor_VAF",
    "tumor_var_avg_base_quality",
    "tumor_var_avg_clipped_length",
    "tumor_var_avg_num_mismatches_as_fraction",
    "tumor_var_avg_sum_mismatch_qualities",
    "tumor_var_med_base_quality",
    "tumor_var_med_frag_len",
    "window_dup_frac",
    "window_improper_frac",
    "window_median_frag_len",
    "window_read_filter_frac"
    ]
    features_to_keep_both = ['tumor_depth']

    for feature in set([*features_to_scale, *features_to_keep_both]):
        variants[f'{feature}_scaled'] = variants.groupby('Sample')[feature].transform(
            lambda x: np.arcsinh(x))

    variants.drop(features_to_scale, axis=1, inplace=True)
    return variants

def load_variants(directory, labels_path):
    ## Need to check that features are being saved with _extracted_features.csv suffix
    pattern='*_extracted_features.csv'

    if os.path.isfile(directory): # in case I messed up and passed a single file instead of a directory
        pattern=''

    files = glob.glob(directory + pattern, recursive=True)

    try: 
        features = pd.concat(
            (pd.read_csv(path).fillna(0) for path in files), ignore_index=True)
    except Exception as e:
        print(f"An error occurred: {e}")
        print(f"There is an issue with one of your files.")

    labels = pd.read_csv(labels_path, sep='\t').fillna(0)
    labels['Label'] = labels['Label'].apply(lambda x: 1 if str(x).lower() in ['real', '1'] else 0)
    labels = labels[['Sample', 'Variant', 'Label']] # In case the Cohort column got saved accidentally

    variants = pd.merge(features, labels, on=['Sample', 'Variant'], how='left').fillna(0)

    ## Convert hot encodings to tri- and penta-nuc context features (if applicable); 
    ## Scale relevant features 
    if 'hot_encoded_ref_base' in variants.columns:
        variants = convert_hot_encodings(variants)
    variants = scale_features(variants)
    return variants

# Accesory function to identify non-contributory features across multiple EBM models
# These features need to be removed to make models compatible for merging
def get_zero_score_names(models):
    zero_score_names = set()
    for model in models:
        zero_score_names.update(
            name for name, scores in zip(model.term_names_, model.term_scores_)
            if np.sum(np.abs(scores)) == 0
        )
    return zero_score_names
