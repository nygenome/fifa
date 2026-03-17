## Need to change data + model paths to relative paths / 
## Add them to the github repo 

import os
import pickle
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score, precision_score, f1_score, recall_score
from interpret.glassbox import ExplainableBoostingClassifier
import multiprocessing as mp
from multiprocessing import Pool
from json import dump
from scipy.stats import trim_mean
from sklearn.preprocessing import FunctionTransformer, StandardScaler
import sys
import glob
import logging
from interpret.glassbox import merge_ebms
from classify_with_scaling import get_zero_score_names, remove_zero_score_terms_by_name # in future can move these functions to a utils file
from helper_funcs import convert_hot_encodings, scale_features

np.random.seed(42)
label_encoder = LabelEncoder() 
logger = logging.getLogger(__name__)

def merge_ebms(input_models, output_path=None, output_pdfs=None):
    models = []

    for i, model_path in enumerate(input_models.split(',')):
        ## Will change these to point to the models on the github
        if model_path == "NYGC1":
            with open('/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/ebm_hyperparams_DLBCL.pkl', 'rb') as file:
                model = pickle.load(file)
            models.append(model)
        elif model_path == "NYGC2":
            with open('/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/ebm_hyperparams_ROT.pkl', 'rb') as file:
                model = pickle.load(file)
            models.append(model)
        elif model_path == "CGCI-BLGSP":
            with open('/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/ebm_hyperparams_CGCI-BLGSP.pkl', 'rb') as file:
                model = pickle.load(file)
            models.append(model)
        elif model_path == "CGCI-HTMCP":
            with open('/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/ebm_hyperparams_CGCI-HTMCP.pkl', 'rb') as file:
                model = pickle.load(file)
            models.append(model)
        elif not os.path.isfile(model_path):
            raise FileNotFoundError(f"Model file {model_path} does not exist.")
        else:
            model = pickle.load(open(model_path, 'rb'))
            models.append(model)

    if models.is_empty():
        raise ValueError("No models were loaded. Please check the input model paths.")

    zero_score_names = get_zero_score_names(models)

    if len(models) > 1:
        for model in models:
            model = remove_zero_score_terms_by_name(model, zero_score_names)
            model_feature_index = {name: i for i, name in enumerate(model.feature_names_in_)}
            reorder_indices = [model_feature_index[name] for name in models[0].feature_names_in_]
            
            model.feature_names_in_ = [model.feature_names_in_[i] for i in reorder_indices]
            model.term_names_ = [model.term_names_[i] for i in reorder_indices]
            model.feature_types_in_ = [model.feature_types_in_[i] for i in reorder_indices]

    ebm = merge_ebms(models)

    if output_path: 
        logger.info(f'Output path provided; Merged model saved to: {output_path}')
        with open(output_path, 'wb') as file:
            pickle.dump(ebm, file)


    ## Test in HCC1395
    # Obviously will need to change this path 
    directory = "/gpfs/commons/groups/compbio/projects/FFPE_filtering/seqc2/extracted_features/"
    pattern = "*_extracted_features.csv"

    matching_files = glob.glob(directory + pattern, recursive=True)
    try: 
        features = pd.concat((pd.read_csv(f).fillna(0) for f in matching_files), ignore_index=True)
        features = features.drop(['Cohort_x', 'Cohort_y', '0'], axis=1)
        features['Cohort'] = 'HCC1395'
        features['tumor_depth'] = 0
    except Exception as e:
        print(f"An error occurred: {e}")
    labels_path="/gpfs/commons/groups/compbio/projects/FFPE_filtering/seqc2/replicating_deepsom/real_labels.csv"
    labels = pd.read_csv(labels_path, sep=',').fillna(0)
    labels['Label'] = labels['Label'].apply(lambda x: 1 if x == 'Real' else 0)
    variants = pd.merge(features, labels, on=['Sample', 'Variant'], how='inner')

    variants = convert_hot_encodings(variants)
    variants = scale_features(variants)

    logger.info('Uploaded cell line test samples')

    ## Generate PDF plots of the performance of the merged model on HCC1395
    ## Do these need to be plots ? Or can they be just summary metrics ? 
    ## Could probably also generate a confusion matrix here too
    logger.info('Generating Plots showing performance of merged FIFA model')
    if os.path.isfile(output_pdfs):
        logger.warning(f"Output PDF path {output_pdfs} already exists. Overwriting.")

    ## Testing on chr1 variants from HCC1395 with Mutect2
    print(f'Metrics in HCC1395 Cell Line Replicates (chr1 variants)')
    for sample in variants[variants['Cohort'] == 'HCC1395']['Sample'].unique():
        testing_set = variants[(variants['Sample'] == sample) & (variants['Variant'].str.startswith('chr1:'))].drop(['Sample', 'Variant', 'Cohort'], axis=1)

        X_test = testing_set.drop(['Label'], axis=1)
        y_test = label_encoder.fit_transform(testing_set['Label'])
        
        y_pred = ebm.predict(X_test)
        
        print(f'Mutect2 calls from {sample}')
        f1 = f1_score(y_test, y_pred)
        print("F1: {:.3f}".format(f1))
        
        accuracy = accuracy_score(y_test, y_pred)
        print("Accuracy: {:.3f}".format(accuracy))
        
        precision = precision_score(y_test, y_pred)
        print("Precision: {:.3f}".format(precision))
        
        recall = recall_score(y_test, y_pred)
        print("Recall: {:.3f}".format(recall))

    ## Testing on DeepSomatic's chr 1 variants
    ## Again this is code that I can clean up

    features=pd.read_csv("/gpfs/commons/groups/compbio/projects/FFPE_filtering/seqc2/replicating_deepsom/deepsomatic-ffpe-wgs-case-study/1395_tumor_ffpe_wgs_extracted_features.csv").fillna(0)
    labels_path="/gpfs/commons/groups/compbio/projects/FFPE_filtering/seqc2/replicating_deepsom/deepsomatic-ffpe-wgs-case-study/github_download/labels.txt"
    column_names=['CHROM', 'POS', 'REF', 'ALT', 'Label']
    labels = pd.read_csv(labels_path, sep='\t', names=column_names).fillna(0)
    labels['Variant'] = labels.apply(lambda row: f"{row['CHROM']}:{row['POS']}_{row['REF']}>{row['ALT']}", axis=1)
    labels = labels.drop(['CHROM', 'POS', 'REF', 'ALT'], axis=1)
    labels['Sample'] = '1395_tumor_ffpe_wgs'
    labels['Label'] = labels['Label'].apply(lambda x: 1 if x == 'Real' else 0)
    variants = pd.merge(features, labels, on=['Sample', 'Variant'], how='inner')

    ## Convert hot encodings to tri- and penta-nuc context features (if applicable); 
    ## Scale relevant features 
    if 'hot_encoded_ref_base' in variants.columns:
        variants = convert_hot_encodings(variants)
    variants = scale_features(variants)

    for sample in variants[variants['Cohort'] == 'DeepSomatic']['Sample'].unique():
        testing_set = variants[(variants['Sample'] == sample) & (variants['Variant'].str.startswith('chr1:'))].drop(['Sample', 'Variant', 'Cohort'], axis=1)

        X_test = testing_set.drop(['Label'], axis=1)
        y_test = label_encoder.fit_transform(testing_set['Label'])
        y_pred = ebm.predict(X_test)

        print(f'Metrics in DeepSomatic\'s Held-Out Test Set (chr1 variants)')
        f1 = f1_score(y_test, y_pred)
        print("F1: {:.3f}".format(f1))
        
        accuracy = accuracy_score(y_test, y_pred)
        print("Accuracy: {:.3f}".format(accuracy))
        
        precision = precision_score(y_test, y_pred)
        print("Precision: {:.3f}".format(precision))
        
        recall = recall_score(y_test, y_pred)
        print("Recall: {:.3f}".format(recall))

