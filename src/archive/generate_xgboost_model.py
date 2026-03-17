import os
import sys
import pickle
import numpy as np
import xgboost as xgb
import pandas as pd
import vtreat
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.preprocessing import LabelEncoder
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import MinMaxScaler,RobustScaler,QuantileTransformer
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats.mstats import winsorize

def treat_and_scale_numerical(data, scale_columns):
    """
    Scales select numerical features in the data. 

    Parameters:
    data (pd.DataFrame): Input DataFrame of variants and features to be scaled.
    scale_columns (list of str): List of numerical features to be scaled
    categorical_columns (list of str): List of categorical features (exclude from scaling)

    Returns:
    pd.DataFrame: A DataFrame combining scaled and unscaled numerical features for all variants
    """
    scale_columns = [col for col in scale_columns if col in data.columns]
    categorical_columns = ['left_two_base', 'left_one_base', 'hot_encoded_ref_base', 'right_one_base', 'right_two_base', 'hot_encoded_var_base']
    
    for col in categorical_columns:
        data[col] = data[col].astype('category')
        
    def scale_features(group):
        if scale_columns:
            scaler = MinMaxScaler()  #RobustScaler(quantile_range=(5, 95)) #
            group[scale_columns] = group[scale_columns].apply(lambda x: winsorize(x.to_numpy(), limits=[0.05, 0.05]))
            group[scale_columns] = scaler.fit_transform(group[scale_columns])
        return group
    
    scaled_data = data.groupby('Cohort').apply(scale_features).reset_index(drop=True)
    
    return scaled_data

def perform_nfold_CV(X_train, y_train, n):
  param_distributions = {'n_estimators': [10, 50, 70], 
                        'max_depth': [None, 3, 5, 8, 10],
                        'eta': [0.01, 0.08, 0.1, 0.3], #learning_rate
                        'min_child_weight': [1, 3, 5, 7],
                        'subsample': [0.7, 0.8, 1.0],
                        'colsample_bytree': [0.8, 0.95, 1.0]
                      }

  model = xgb.XGBClassifier(objective='binary:logistic', eval_metric=['logloss', 'auc'], enable_categorical=True,  tree_method="hist")

  #Using N-fold CV
  random_search = RandomizedSearchCV(estimator=model, param_distributions=param_distributions, cv=n, verbose=3, n_jobs=4)
  
  # Perform the search
  random_search.fit(X_train, y_train)
  
  best_n_estimators = random_search.best_params_['n_estimators']
  best_max_depth = random_search.best_params_['max_depth']
  best_eta = random_search.best_params_['eta']
  best_min_child_weight = random_search.best_params_['min_child_weight']
  best_subsample = random_search.best_params_['subsample']
  best_colsample_bytree = random_search.best_params_['colsample_bytree']
  
  return best_n_estimators, best_max_depth, best_eta, best_min_child_weight, best_subsample, best_colsample_bytree

np.random.seed(42) # maybe can set this through passing in an argument? 
label_encoder = LabelEncoder()

cohort_file = sys.argv[1]
model_file = sys.argv[2]

variants = pd.read_csv(cohort_file)
variants['Cohort'] = 'ROT' #manually adding this in right now because welp 

#For adding in predictions from Mobster
#tail_file = sys.argv[3]
#tail_scores = pd.read_csv(tail_file).reset_index(drop=True)
#tail_scores['POS'] = tail_scores['POS'].astype(int)
#tail_scores['Variant'] = tail_scores.apply(lambda row: f"{row['CHROM']}:{row['POS']}_{row['REF']}>{row['ALT']}", axis=1)
#tail_scores = tail_scores[['sample', 'Variant', 'Tail']]
#variants = pd.merge(variants, tail_scores, left_on=['Sample', 'Variant'], right_on=['sample', 'Variant'], how='inner').drop('sample', axis=1)

variants = variants.fillna(0)
scale_columns = variants.select_dtypes(include=np.number).columns.tolist()

X = variants.drop(['Label'], axis=1)
y = label_encoder.fit_transform(variants['Label'])

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

X_train = treat_and_scale_numerical(X_train, scale_columns)
X_test = treat_and_scale_numerical(X_test, scale_columns)

variant_ids_train = X_train[['Sample', 'Variant', 'Cohort']]
variant_ids_test = X_test[['Sample', 'Variant', 'Cohort']]
X_train = X_train.drop(['Sample', 'Variant', 'Cohort'], axis=1)
X_test = X_test.drop(['Sample', 'Variant', 'Cohort'], axis=1)

dtrain = xgb.DMatrix(X_train, y_train, enable_categorical=True)
dtest = xgb.DMatrix(X_test, y_test, enable_categorical=True)

best_n_estimators, best_max_depth, best_eta, best_min_child_weight, \
best_subsample, best_colsample_bytree = perform_nfold_CV(X_train, y_train, 10)

params = {
    "max_depth": best_max_depth, 
    "eta": best_eta, 
    'min_child_weight': best_min_child_weight, 
    'subsample': best_subsample,
    'colsample_bytree': best_colsample_bytree,
    "objective": "binary:logistic",
    "enable_categorical": True,
    "tree_method": "hist"
}

evals = [(dtest, "validation"), (dtrain, "train")]

bst = xgb.XGBClassifier(**params, num_boost_round = best_n_estimators)
bst.fit(X_train, y_train, eval_set=[(X_train, y_train)])

bst.save_model(model_file)

#Test on Training Data
#Upper Bound for Accuracy 
preds = bst.predict(X_train)
threshold = 0.5
binary_preds = (preds >= threshold).astype(int)

labels = y_train.astype(int)
accuracy = sum(1 for i in range(len(binary_preds)) if int(binary_preds[i]) == labels[i]) / float(len(binary_preds))
print("Accuracy on training data=%f" % accuracy)

#Test on Testing Data
preds = bst.predict(X_test)
threshold = 0.5
binary_preds = (preds >= threshold).astype(int)

labels = y_test.astype(int)

accuracy = sum(1 for i in range(len(binary_preds)) if int(binary_preds[i]) == labels[i]) / float(len(binary_preds))
precision = sum(1 for i in range(len(binary_preds)) if int(binary_preds[i]) == 1 == labels[i]) / sum(1 for i in range(len(binary_preds)) if int(binary_preds[i]) == 1)
print("Accuracy on test data=%f" % accuracy)
print("Precision on test data=%f" % precision)
  
