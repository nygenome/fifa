import os
import pickle
import numpy as np
import sys
import pandas as pd
from sklearn.preprocessing import LabelEncoder
import glob
from sklearn.metrics import accuracy_score, precision_score, roc_auc_score, recall_score, f1_score
from interpret.glassbox import ExplainableBoostingClassifier
import logging
from helper_funcs import convert_hot_encodings, scale_features, load_variants

#global logger
np.random.seed(42)
label_encoder = LabelEncoder()

def retrain(directory, labels_path, output_path):
    """
    Train new EBM model with user-inputted cohort. 

    Args:
    directory (str): Directory with extracted features files for all samples.
    labels_path (str): Path to TSV file with variant labels
    output_path (str): Path to save new EBM model

    Returns:
    None: Saves new EBM model to output_path
    """
    variants = load_variants(directory, labels_path)
    
    ## Not sure when this error had come up 
    #if 'Cohort_x' in variants.columns:
    #    variants.drop(['Cohort_x', 'Cohort_y', 'sample'],axis=1,inplace=True)

    X = variants.drop(['Sample', 'Variant', 'Cohort', 'Label'], axis=1)
    y = label_encoder.fit_transform(variants['Label'])

    del variants

    weight=sum(y == 1) / len(y)
    w = np.array([(1-weight) if label == 1 else weight for label in y])
    ebm = ExplainableBoostingClassifier()
                                    
    ebm.fit(X, y, sample_weight=w)

    with open(output_path, 'wb') as file:
        pickle.dump(ebm, file)

