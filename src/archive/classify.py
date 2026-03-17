#!/bin/env python

##########################################################################
#	USAGE: classify variants in a test set
#
#   DESCRIPTION: Original script for classifying variants with an XGBoost model.
#   Also just a graveyard of functions that I didn't want to lose.
#
#   Created by Valentina Grether
##########################################################################

import os
import pickle
import numpy as np
import xgboost as xgb
import pandas as pd
import vtreat
import csv
import re
import pysam
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.preprocessing import LabelEncoder
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import MinMaxScaler,RobustScaler,QuantileTransformer
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats.mstats import winsorize
from process_bam_file import get_column_idx

global predictions_path

np.random.seed(42)

def eval_preds(preds, labels, variant_ids):
    threshold = 0.5
    binary_preds = (preds >= threshold).astype(int)

    # Convert predictions to integers just to be safe 
    binary_preds = [int(pred) for pred in binary_preds]
    labels = [int(label) for label in labels]

    # Check for None values
    assert all(pred is not None for pred in binary_preds), "Predictions contain None values"
    assert all(label is not None for label in labels), "Labels contain None values"

    accuracy = sum(1 for i in range(len(binary_preds)) if binary_preds[i] == labels[i]) / float(len(binary_preds))
    precision = sum(1 for i in range(len(binary_preds)) if binary_preds[i] == 1 == labels[i]) / sum(1 for i in range(len(binary_preds)) if binary_preds[i] == 1)
    
    return accuracy, precision

def make_predictions(bst, variants, cols_to_scale):
    """
    Scale features, and make predictions on the variants. 
    Write out model predictions on each variant to a csv file. 

    Args:
        bst (xgboost.XGBClassifier): xgboost model 
        variants (pd.DataFrame): data frame with variants and all features
        labeled (bool): True if dataset includes Real/Artifact labels, 
            False if making predictions on unlabeled data. 
        cols_to_scale (lst): List of columns to scale with MinMax Scaler. 
            Default is all numerical columns. 
    """   
    labeled = 'Label' in variants.columns

    one_hot_to_base = {
    '[1. 0. 0. 0. 0.]': 'A', #A
    '[0. 1. 0. 0. 0.]': 'C', #C
    '[0. 0. 1. 0. 0.]': 'G', #G 
    '[0. 0. 0. 1. 0.]': 'T', #T
    '[0. 0. 0. 0. 1.]': 'N', #N
    '[0. 0. 0. 0. 0.]': 'N'
    }

    categorical_columns = ['left_two_base', 'left_one_base', 'hot_encoded_ref_base', 'right_one_base', 'right_two_base', 'hot_encoded_var_base']

    for col in categorical_columns:
        variants[col] = variants[col].apply(lambda x: one_hot_to_base.get(str(x), x))
    variants = variants[~variants[categorical_columns].apply(lambda x: x.str.contains('N')).any(axis=1)]

    variants.loc[:, 'trinucleotide_context'] = (
        variants['left_one_base'] + '[' + 
        variants['hot_encoded_ref_base'] + '>' + 
        variants['hot_encoded_var_base'] + ']' + 
        variants['right_one_base']
    )

    variants.loc[:, 'pentanucleotide_context'] = (
        variants['left_two_base'] + 
        variants['left_one_base'] + '[' + 
        variants['hot_encoded_ref_base'] + '>' + 
        variants['hot_encoded_var_base'] + ']' + 
        variants['right_one_base'] + 
        variants['right_two_base']
    )

    variants = variants.drop(categorical_columns, axis=1)

    categorical_columns = ['trinucleotide_context', 'pentanucleotide_context']
    variants[categorical_columns] = variants[categorical_columns].astype('category')

    variant_ids_test = variants[['Sample', 'Variant', 'Cohort']].reset_index(drop=True)

    model_feature_names = bst.feature_names
    X_test = variants.drop(['Sample', 'Variant', 'Cohort'], axis=1)

    for col in set(model_feature_names) - set(X_test.columns):
        X_test[col] = np.nan

    X_test = X_test[model_feature_names]

    dtest = xgb.DMatrix(X_test, enable_categorical=True, missing=np.nan)
    binary_preds = (bst.predict(dtest) >= 0.5).astype(int)
    xgb_preds = pd.DataFrame(binary_preds, columns=['xgb_preds']).reset_index(drop=True)
    return pd.concat([variant_ids_test, xgb_preds], axis=1)


## Puting this here for now so that I have it stored, but obviously from classify with scaling script
def make_predictions(ebm, variants):
    """
    Scale features, and make predictions on the variants. 
    Write out model predictions on each variant to a csv file. 

    Args:
        ebm (ExplainableBoostingClassifier): ebm model 
        variants (pd.DataFrame): data frame with variants and all features

    """   
    label_encoder = LabelEncoder() 
    

    if 'Label' in variants.columns:
        y_test = label_encoder.fit_transform(variants['Label'].apply(lambda x: 1 if x == 'Real' else 0))
        variants.drop(columns='Label', inplace=True)

    one_hot_to_base = {
    '[1. 0. 0. 0. 0.]': 'A', #A
    '[0. 1. 0. 0. 0.]': 'C', #C
    '[0. 0. 1. 0. 0.]': 'G', #G 
    '[0. 0. 0. 1. 0.]': 'T', #T
    '[0. 0. 0. 0. 1.]': 'N', #N
    '[0. 0. 0. 0. 0.]': 'N'
    }

    categorical_columns = ['left_two_base', 'left_one_base', 'hot_encoded_ref_base', 'right_one_base', 'right_two_base', 'hot_encoded_var_base']

    for col in categorical_columns:
        variants[col] = variants[col].apply(lambda x: one_hot_to_base.get(str(x), x))

    variants = variants[~variants[categorical_columns].apply(lambda x: x.str.contains('N')).any(axis=1)]
    variants['trinucleotide_context'] = variants['left_one_base'] + '[' + variants['hot_encoded_ref_base'] + '>' + variants['hot_encoded_var_base'] + ']' + variants['right_one_base']
    variants['pentanucleotide_context'] = variants['left_two_base'] + variants['left_one_base'] + '[' + variants['hot_encoded_ref_base'] + '>' + variants['hot_encoded_var_base'] + ']' + variants['right_one_base'] + variants['right_two_base']

    variants = variants.drop(categorical_columns, axis=1)

    categorical_columns = ['trinucleotide_context', 'pentanucleotide_context']
    for col in categorical_columns:
        variants[col] = variants[col].astype('category')
        variants = variants[variants[col].notnull()]

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
    variants_over_30 = variants[variants['tumor_VAF'] > 0.30]
    variants_under_30 = variants[variants['tumor_VAF'] <= 0.30]   

    del variants

    for feature in set([*features_to_scale, *features_to_keep_both]):
        variants_under_30[f'{feature}_scaled'] = variants_under_30.groupby('Sample')[feature].transform(
            lambda x: np.arcsinh(x))

    variants_under_30.drop(features_to_scale, axis=1, inplace=True)

    variant_ids_test = variants_under_30[['Sample', 'Variant', 'Cohort']].reset_index(drop=True)

    X_test = variants_under_30.drop(['Sample', 'Variant', 'Cohort'], axis=1)
    
    if not (set(ebm.feature_names_in_).issubset(set(X_test.columns))):
        print("Warning: Feature mismatch between model and test data. Please check input features.")
        sys.exit(1)

    y_pred = pd.DataFrame(ebm.predict(X_test), columns=['Predicted'])
    under_30 = pd.concat([variant_ids_test, y_pred], axis=1)
    over_30 = variants_over_30[['Sample', 'Variant', 'Cohort']].reset_index(drop=True)
    over_30['Predicted'] = 1
    del variants_over_30
    del variants_under_30
    return pd.concat([under_30, over_30], ignore_index=True)

def generate_output_vcf_file(predictions, sample, vcf_path): 
    global predictions_path 
    binary_preds = predictions['xgb_preds']
    outpath = os.path.join(predictions_path, re.sub(r'\.vcf(\.gz)?$', '.xgboost.vcf', os.path.basename(vcf_path)))
    
    if os.path.isfile(outpath):
        os.remove(outpath)

    try:
        origvcf = pysam.VariantFile(vcf_path, "r")
        origvcf.header.info.add("XGB", "1", "Integer", "XGBoost Classification (0-Artifact or 1-Real)")
        vcf_out = pysam.VariantFile(outpath, 'w', header=origvcf.header)
        
        prediction_dict = dict(zip(predictions['Variant'], predictions['xgb_preds']))
        positions = set()
    
        for variant in origvcf:
            position = f"{variant.chrom}:{variant.pos}_{variant.ref}>{variant.alts[0]}"
            
            if position in prediction_dict:
                if position in positions:
                    print(f"Duplicate entry found for: {position}")
                else:
                    positions.add(position)
                    variant.info['XGB'] = int(prediction_dict[position])
                    vcf_out.write(variant)
        
        print(f"Successfully generated: {outpath}")

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
            generate_output_vcf_file(predictions, sample, vcf_path)
            
def predict(extracted_features_paths, outpath, model_file, pairs_path=None, sample=None, vcf_path=None, cols_to_scale=None):
    global predictions_path 
    try:
        if os.path.isdir(outpath):
            predictions_path = outpath
        elif os.path.isfile(outpath):
            predictions_path = os.path.basename(outpath)
        if not os.path.exists(predictions_path):
            os.makedirs(predictions_path)
    except Exception as e:
        print(f"An error occurred: {e}")
    
    try: 
        variants = pd.concat(
            (pd.read_csv(path).fillna(0) for path in extracted_features_paths), ignore_index=True)
    except Exception as e:
        print(f"An error occurred: {e}")
        print(f"There is an issue with the following file of extracted features.\n{path}")
        
    try:
        bst = xgb.Booster()
        bst.load_model(model_file)
    except Exception as e:
        print(f"An error occurred: {e}")
        print(f"The path for your xgboost model is incorrect. \n{path}")
    
    
    combined_df = make_predictions(bst, variants, cols_to_scale)
     # Bueno valu te recuerdo.... 
     # we determine whether accuracy predictions can be calculated 
     # based on whether variants do/do not have labels
    if(pairs_path):
        if not os.path.isfile(pairs_path):
            print(f"The path for your sample files is incorrect")
            exit(1)
        parse_pairs_path(pairs_path, combined_df) 
    else: 
        generate_output_vcf_file(combined_df, sample, vcf_path)


def treat_and_scale_numerical(data, scale_columns):
    """
    Scales select numerical features in the data. 

    Args:
        data (pd.DataFrame): Input DataFrame of variants and features to be scaled.
        scale_columns (lst): List of numerical features to scale (exclude categorical columns)

    Returns:
        pd.DataFrame: A DataFrame combining scaled and unscaled numerical features for all variants
    """    
    scale_columns = [col for col in scale_columns if col in data.columns]
    categorical_columns = ['left_two_base', 'left_one_base', 'hot_encoded_ref_base', 'right_one_base', 'right_two_base', 'hot_encoded_var_base']
    
    categorical_columns = ['trinucleotide_context', 'pentanucleotide_context']

    for col in categorical_columns:
        if col in data.columns: 
            data[col] = data[col].astype('category')
        
    def scale_features(group):
        if scale_columns:
            scaler = MinMaxScaler() 
            group[scale_columns] = group[scale_columns].apply(lambda x: winsorize(x.to_numpy(), limits=[0.05, 0.05]))
            group[scale_columns] = scaler.fit_transform(group[scale_columns])
        return group
    
    scaled_data = data.groupby('Cohort').apply(scale_features).reset_index(drop=True)
    return scaled_data

#if cols_to_scale:
    #    scale_columns = pd.read_csv(cols_to_scale, header=None, squeeze=True)
    #    if (set(scale_columns) - set(model_feature_names)) or (set(scale_columns) - set(variants.columns)):
    #        print("Warning. Trying to scale columns not found in the model and/or inputted feature list")
    #else: 
    #    scale_columns = variants.select_dtypes(include=np.number).columns.tolist()
    
    #X_test = variants.drop(['Label'], axis=1) if labeled else variants
    #y_test = LabelEncoder().fit_transform(variants['Label']) if labeled else None
    #sample_ids = X_test['Sample']
    #variant_ids_test = pd.DataFrame(X_test['Variant'], columns=['Variant'])
    
    #for col in (set(model_feature_names) - set(X_test.columns)):
    #    X_test[col] = np.nan
    
    #X_test = X_test[model_feature_names]
    
    #preds = bst.predict(X_test)
    #binary_preds = (preds >= 0.5).astype(int)

    #if labeled:
     #   accuracy, precision = eval_preds(preds, y_test.astype(int), variant_ids_test)
    #    print("Accuracy of XGBoost model: %f" % accuracy)
    #    print("Precision of XGBoost model: %f" % precision)

    #sample_ids.reset_index(drop=True, inplace=True)
    #variant_ids_test.reset_index(drop=True, inplace=True)
    #binary_preds = pd.DataFrame(binary_preds, columns=['Model Predictions'])

    #if labeled:
    #    labels_df = pd.DataFrame(y_test, columns=['True Labels'])
    #    combined_df = pd.concat([sample_ids, variant_ids_test, labels_df, binary_preds], axis=1)
    #else:
    #    combined_df = pd.concat([sample_ids, variant_ids_test, binary_preds], axis=1)