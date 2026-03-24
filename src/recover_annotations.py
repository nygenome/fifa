import os
import pickle
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, balanced_accuracy_score, f1_score, recall_score
from interpret.glassbox import ExplainableBoostingClassifier
from interpret.glassbox import merge_ebms
import multiprocessing as mp
from multiprocessing import Pool
from json import dump
from scipy.stats import trim_mean
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
    Write out model predictions on each variant to a csv file. 

    Args:
        ebm (ExplainableBoostingClassifier): ebm model 
        variants (pd.DataFrame): data frame with variants and all features

    """   
    label_encoder = LabelEncoder() 
    

    if 'Label' in variants.columns:
        y_test = label_encoder.fit_transform(variants['Label'].apply(lambda x: 1 if x == 'Real' else 0))
        variants.drop(columns='Label', inplace=True)

    ## Convert hot encodings to tri- and penta-nuc context features (if applicable); 
    ## Scale relevant features 
    if 'hot_encoded_ref_base' in variants.columns:
        variants = convert_hot_encodings(variants)
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

def generate_output_vcf_file(predictions, rna_annotations, vcf_path): 
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
        origvcf.header.info.add("RNA_seq", "3", "Integer", "Number of Reads in the RNA supporting the REF, ALT, or OTHER ALT alleles")
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
            
            if (pred is not None) and (prob is not None): #and variant.info.get('HighConfidence')
                positions.add(position)
                variant.info['FIFA_LABEL'] = int(pred)
                variant.info['FIFA_PROB'] = round(prob, 4)
                if position in rna_annotations['Variant'].values:
                    row = rna_annotations[rna_annotations['Variant'] == position][['REF', 'ALT', 'OTHER_ALT']]
                    
                    ref = int(row['REF'].fillna(0).values[0])
                    alt = int(row['ALT'].fillna(0).values[0])
                    other_alt = int(row['OTHER_ALT'].fillna(0).values[0])
                    variant.info['RNA_seq'] = [ref, alt, other_alt]

                    total = ref + alt + other_alt

                    if (total > 5) and (alt > 0.1 * total):
                        variant.info['FIFA_LABEL'] = 1 
            
            vcf_out.write(variant)

        logger.info(f"Successfully generated: {outpath}")

    except Exception as e:
        print(f"An error occurred: {e}")

## In the future, this code should be factored out from recover_annotations.py and 
## classify_with_scaling.py

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

def predict(extracted_features_paths, outpath, model_file, rna_path=None, sample=None, vcf_path=None):
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
        if isinstance(model_file, list):
            models = []
            for i, model_path in enumerate(model_file.split(',')):
                fifa_install_dir=os.path.dirname(os.getcwd())
                if model_path == "NYGC1":
                    path= os.path.join(fifa_install_dir, 'models/ebm_hyperparams_NYGC1.pkl')
                    with open(path, 'rb') as file:
                        model = pickle.load(file)
                    models.append(model)
                elif model_path == "NYGC2":
                    path= os.path.join(fifa_install_dir, 'models/ebm_hyperparams_NYGC2.pkl')
                    with open(path, 'rb') as file:
                        model = pickle.load(file)
                    models.append(model)
                elif re.match(r'^(?:CGCI-)?BLGSP$', model_path):
                    path= os.path.join(fifa_install_dir, 'models/ebm_hyperparams_CGCI-BLGSP.pkl')
                    with open(path, 'rb') as file:
                        model = pickle.load(file)
                    models.append(model)
                elif re.match(r'^(?:CGCI-)?HTMCP$', model_path):
                    path= os.path.join(fifa_install_dir, 'models/ebm_hyperparams_CGCI-HTMCP.pkl')
                    with open(path, 'rb') as file:
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

                ebm =  merge_ebms(models)
        elif os.path.isfile(model_file):
            ebm = ExplainableBoostingClassifier()
            with open(model_file, 'rb') as f:
                ebm = pickle.load(f)
        
            
    except Exception as e:
        print(f"An error occurred: {e}")
        print(f"The path for your ebm model is incorrect. \n{model_file}")
    
    combined_df = make_predictions(ebm, variants)
    rna_annotations = pd.read_csv(rna_path).fillna(0).reset_index()
    rna_annotations = rna_annotations[rna_annotations['SAMPLE'] == sample]

    generate_output_vcf_file(combined_df, rna_annotations, vcf_path)