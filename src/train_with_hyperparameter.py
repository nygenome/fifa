import optuna
import os
import pickle
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, balanced_accuracy_score, f1_score
from interpret.glassbox import ExplainableBoostingClassifier
from sklearn.model_selection import StratifiedKFold
import time
import multiprocessing as mp
from multiprocessing import Pool
from json import dump
from concurrent.futures import ProcessPoolExecutor
import sys
from helper_funcs import load_variants

np.random.seed(42)

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

    f1 = f1_score(y_test, y_preds)
    return f1

def optimize(X_train, X_test, y_train, y_test, fold):
    study_time = time.time()
    study = optuna.create_study(direction="maximize", study_name=f"Fold_{fold}")

    study.optimize(lambda trial: objective(trial, X_train, X_test, y_train, y_test), 
                   n_trials=10, 
                   n_jobs = 1, 
                   timeout=6000, 
                   gc_after_trial=True)
    
    print("Number of finished trials: ", len(study.trials))
    trial = study.best_trial

    study_end = time.time() - study_time
    print(f"Fold {fold} finished in {study_end} seconds.")

    return [trial.value, trial.params]

def run(directory, labels_path, model_file):
    label_encoder = LabelEncoder() 
    variants = load_variants(directory=directory, labels_path=labels_path)
    
    X = variants.drop(['Sample', 'Variant', 'Cohort', 'Label'], axis=1)
    y = label_encoder.fit_transform(variants['Label'])
    
    del variants
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

    print("Final FIFA Model Saved")

def retrain(directory, labels_path, output_path): 
    run(directory=directory, labels_path=labels_path, model_file=output_path)