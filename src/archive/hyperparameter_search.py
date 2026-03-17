import optuna
from sklearn.model_selection import train_test_split
import os
import pickle
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.preprocessing import LabelEncoder
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import MinMaxScaler,RobustScaler,QuantileTransformer
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, balanced_accuracy_score, f1_score
from interpret.glassbox import ExplainableBoostingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
import time
import multiprocessing as mp
from multiprocessing import Pool
from json import dump
from sklearn.preprocessing import FunctionTransformer, StandardScaler
from concurrent.futures import ProcessPoolExecutor
import sys

np.random.seed(42)

def load_and_merge_features(features_path, tail_scores_path, cohort_name):
    features = pd.read_csv(features_path).fillna(0)
    features['Cohort'] = cohort_name
    tail_scores = pd.read_csv(tail_scores_path).reset_index(drop=True)
    tail_scores['POS'] = tail_scores['POS'].astype(int)
    tail_scores['Variant'] = tail_scores.apply(lambda row: f"{row['CHROM']}:{row['POS']}_{row['REF']}>{row['ALT']}", axis=1)
    tail_scores = tail_scores[['sample', 'Variant', 'Tail']]
    merged = pd.merge(features, tail_scores, left_on=['Sample', 'Variant'], right_on=['sample', 'Variant'], how='inner')
    merged.drop(columns='sample', inplace=True)
    merged['Label'] = merged['Label'].apply(lambda x: 1 if x == 'Real' else 0)
    return merged

def load_and_merge_labels(directory, labels_path):
    pattern='*_extracted_features.csv'
    files = glob.glob(directory + pattern, recursive=True)
    try:
        features = pd.concat((pd.read_csv(f).fillna(0) for f in files), ignore_index=True)
    except Exception as e:
        print(f"Error reading features: {e}\nProblematic dir: {directory}")
        return pd.DataFrame()
    labels = pd.read_csv(labels_path, sep='\t').fillna(0)
    labels['Label'] = labels['Label'].apply(lambda x: 1 if x == 'Real' else 0)
    return pd.merge(features, labels, on=['Sample', 'Variant'], how='inner')

def scaling_iteration():
    variants_BLGSP = load_and_merge_labels(
    "/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/blgsp/extracted_features/",
    "/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/blgsp/labels.txt"
    )

    variants_HTMCP = load_and_merge_labels(
        "/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/all_htmcp/extracted_features/",
        "/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/all_htmcp/labels.txt"
        )
    
    features_path = '/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/DLBCL/extracted_features/DLBCL_all_features.csv'
    variants_DLBCL = pd.read_csv(features_path).fillna(0)

    features_path = '/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/ROT/FFPE/extracted_features/ROT_all_features.csv'
    variants_ROT = pd.read_csv(features_path).fillna(0) 
    
    variants_HTMCP['Cohort'] = 'CGCI-HTMCP'
    columns_order = variants_BLGSP.columns
    variants_DLBCL = variants_DLBCL[columns_order]
    variants_ROT = variants_ROT[columns_order]
    variants_HTMCP = variants_HTMCP[columns_order]
    
    variants=pd.concat([variants_BLGSP, variants_DLBCL, variants_HTMCP, variants_ROT], ignore_index=True) #variants_HTMCP  ,variants_DLBCL
    del variants_BLGSP, variants_DLBCL, variants_HTMCP, variants_ROT
    
    ### REMINDER TO REMOVE THIS LINE NEXT TIME!!!!!!!
    #variants = variants[variants['tumor_VAF'] <= 0.3]

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
    '''
    variants = variants.drop(['tumor_ref_avg_se_mapping_quality',
                        'tumor_var_avg_se_mapping_quality'], axis=1)
    '''
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

    for feature in set([*features_to_scale, *features_to_keep_both]):
        variants[f'{feature}_scaled'] = variants.groupby('Sample')[feature].transform(
            lambda x: np.arcsinh(x))
    
    variants.drop(features_to_scale, axis=1, inplace=True)
    
    '''
    data_training = variants[(variants['Cohort'] == testing_cohort)]
    data_testing = variants[(variants['Cohort'] != testing_cohort)]
    variants = scale_data(data_training, data_testing)
    '''
    
    print("habemus datos")

    return variants

def objective(trial, X_train, X_test, y_train, y_test):    
    max_bins= trial.suggest_categorical("max_bins", [32, 256, 512, 1024])
    max_leaves= trial.suggest_categorical("max_leaves", [2, 3, 5])
    early_stopping_rounds= trial.suggest_categorical("early_stopping_rounds", [50, 100])
    early_stopping_tolerance= trial.suggest_categorical("early_stopping_tolerance", [0, 1e-5, 5e-5])
    learning_rate= trial.suggest_categorical("learning_rate", [0.005, 0.015])

    ebm = ExplainableBoostingClassifier(max_bins=max_bins, 
                                        max_leaves=max_leaves, 
                                        early_stopping_rounds=early_stopping_rounds, 
                                        early_stopping_tolerance=early_stopping_tolerance, 
                                        learning_rate=learning_rate)
    ebm.fit(X_train, y_train) 

    y_preds = ebm.predict(X_test)

    # Suggestion from ChatGPT : Use probabilities for AUC
    # ebm.predict_proba(X_test)[:, 1]  # Use probabilities for AUC

    bal_auc = f1_score(y_test, y_preds)
    return bal_auc

def optimize(X_train, X_test, y_train, y_test, fold):
    study_time = time.time()
    study = optuna.create_study(direction="maximize", study_name=f"Fold_{fold}")

    study.optimize(lambda trial: objective(trial, X_train, X_test, y_train, y_test), 
                   n_trials=10, #change here
                   n_jobs = 1, 
                   timeout=6000, 
                   gc_after_trial=True)
    
    print("Number of finished trials: ", len(study.trials))
    trial = study.best_trial

    study_end = time.time() - study_time
    print(f"Fold {fold} finished in {study_end} seconds.")

    return [trial.value, trial.params]

def run(training_cohort, model_file):
    label_encoder = LabelEncoder() 
    variants = scaling_iteration()

    print("final feature set merged!!!")

    training_set = variants[variants['Cohort'] == training_cohort].drop(['Sample', 'Variant', 'Cohort'], axis=1)
    testing_set = variants[variants['Cohort'] != training_cohort].drop(['Sample', 'Variant', 'Cohort'], axis=1)

    # Let's try 5-fold CV
    X = training_set.drop(['Label'], axis=1)
    y = label_encoder.fit_transform(training_set['Label'])
    
    params = {}
    k=5
    skf = StratifiedKFold(n_splits=k)

    begin= time.time()
    
    with ProcessPoolExecutor(max_workers=k) as executor:
        futures = [
            executor.submit(optimize, X.iloc[train], X.iloc[test], y[train], y[test], fold)
            for fold, (train, test) in enumerate(skf.split(X, y), 1)]
        results = [f.result() for f in futures]

    end = time.time() - begin

    best_performance = 0 
    print(f"{end} seconds elapsed for all folds")

    for result in results:
        performance = result[0]
        parameters = result[1]
        if performance > best_performance:
            best_performance = performance
            params = parameters

    print(f"Best Precision Score: {best_performance}. In total, testing took {end} seconds")

    #cont = pd.Series('continuous').repeat(59)
    #nom = pd.Series(['nominal', 'nominal'])
    #types= pd.concat((cont, nom), ignore_index=True)

    weight=sum(y == 1) / len(y)
    w = np.array([(1-weight) if label == 1 else weight for label in y])

    ebm = ExplainableBoostingClassifier(max_bins=params['max_bins'], 
                                        max_leaves=params['max_leaves'], 
                                        early_stopping_rounds=params['early_stopping_rounds'], 
                                        early_stopping_tolerance=params['early_stopping_tolerance'], 
                                        learning_rate=params['learning_rate'],
                                        random_state=42)


    ebm.fit(X, y, sample_weight=w)

    with open(model_file, 'wb') as file:
        pickle.dump(ebm, file)

    print("Model Saved")

    X_test = testing_set.drop(['Label'], axis=1)
    y_test = label_encoder.fit_transform(testing_set['Label'])
    y_pred = ebm.predict(X_test)

    print(f"Model trained with {training_cohort}:")

    auc = roc_auc_score(y_test, y_pred)
    print("AUC: {:.3f}".format(auc))

    accuracy = accuracy_score(y_test, y_pred)
    print("Accuracy: {:.3f}".format(accuracy))

    precision = precision_score(y_test, y_pred)
    print("Precision: {:.3f}".format(precision))

if __name__ == "__main__":
    model_file = sys.argv[1]
    training_cohort = sys.argv[2]
    run(training_cohort=training_cohort, model_file=model_file)

    