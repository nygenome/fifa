import os
import pickle
import numpy as np
import xgboost as xgb
import pandas as pd
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import MinMaxScaler,RobustScaler,QuantileTransformer
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats.mstats import winsorize

from src.old_scripts.classify import treat_and_scale_numerical, eval_preds

def update_xgboost_model(bst, variants, scale_columns):
  """
  Option to retrain xgboost model and incorporate additional training data. 
  Features in the new data are first scaled (assume that all numerical columns are being scaled).

  Args:
  bst (xgboost.Booster): xgboost model 
  variants (pd.DataFrame): DataFrame with variants and all features
  scale_columns (lst): List of columns to scale with MinMax Scaler. 
  label (str): VCF INFO Field with Variant Label

  Returns:
  xgboost.Booster: New xgboost model, retrained with user input data
  """   
  X = variants.drop(['Label'], axis=1)
  
  y = variants['Label']

  X_train_new, X_test_new, y_train_new, y_test_new = train_test_split(X, y, test_size=0.2, random_state=42)

  X_train_new = treat_and_scale_numerical(X_train_new, scale_columns)
  X_test_new = treat_and_scale_numerical(X_test_new, scale_columns)

  variant_ids_train = X_train_new[['Sample', 'Variant', 'Cohort']]
  variant_ids_test = X_test_new[['Sample', 'Variant', 'Cohort']]
  X_train_new = X_train_new.drop(['Sample', 'Variant', 'Cohort'], axis=1)
  X_test_new = X_test_new.drop(['Sample', 'Variant', 'Cohort'], axis=1)

  dtrain = xgb.DMatrix(X_train_new, label=y_train_new, enable_categorical=True)
  dtest = xgb.DMatrix(X_test_new, label=y_test_new, enable_categorical=True)

  new_bst = xgb.train(
    #bueno esto lo resolvere mas adelante que se yo lol
    params={
        "process_type": "update",
        "tree_method": "hist",
        "updater": "refresh",
        "refresh_leaf": True
    },
    dtrain=dtrain,
    xgb_model=bst,
    num_boost_round=10, 
    evals=[(dtrain, "train")]
  )
  #new_bst = bst.fit(X_train_new, y_train_new, eval_set=[(X_train_new, y_train_new)])
  preds = new_bst.predict(dtrain)

  try:
    accuracy_train, precision_train = eval_preds(preds, y_train_new.astype(int), variant_ids_train)
  except Exception as e:
    print(f"Error during training evaluation: {e}\nLikely caused by no \'Real\' predictions")

  # Note: This will overwrite the predictions file so that only the predictions on the test set are saved.
  preds = new_bst.predict(dtest)
  try:
    accuracy_test, precision_test = eval_preds(preds, y_test_new.astype(int), variant_ids_test)
    print("Accuracy of new XGBoost model: %f" % accuracy_test)
    print("Precision of new XGBoost model: %f" % precision_test)
  except Exception as e:
    print(f"Error during testing evaluation: {e}\nLikely caused by no \'Real\' predictions")

  return new_bst

def retrain(model_path, extracted_features_paths, cols_to_scale, outpath):    
  path_checks = {
    model_path: "The path for your xgboost model is incorrect",
    extracted_features_paths: "The path for the file containing variants and their features is incorrect",
    cols_to_scale: "The path pointing to the columns to be scaled by the xgboost model cannot be found"
  }

  for path, error_message in path_checks.items():
    if not os.path.isfile(path):
      print(error_message)
      exit(1)

  bst = xgb.Booster()
  bst.load_model(model_path)

  variants = pd.read_csv(extracted_features_paths).reset_index(drop=True).fillna(0)
  variants = variants.drop(['Tail'], axis=1)
  #REMOVE THIS, JUST BECAUSE ORIGINAL FILE DIDN'T HAVE IT!!!
  #variants['Cohort'] = 'DLBCL'
  
  new_bst = update_xgboost_model(bst, variants, cols_to_scale)
  new_bst.save_model(outpath + '/retrained_xgboost_model.json')
