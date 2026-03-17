import os
import pickle
import numpy as np
import xgboost as xgb
import pandas as pd
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.preprocessing import LabelEncoder
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import MinMaxScaler,RobustScaler,QuantileTransformer
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
from scipy.stats.mstats import winsorize
import seaborn as sns
import glob
from sklearn.metrics import accuracy_score, precision_score, roc_auc_score
from interpret.glassbox import ExplainableBoostingClassifier
from interpret import show
from interpret.glassbox import merge_ebms
from interpret import set_visualize_provider
from interpret.provider import InlineProvider
set_visualize_provider(InlineProvider())


np.random.seed(42)
label_encoder = LabelEncoder()

features_path = '/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/train_xgboost_model/DLBCL_features_all_chroms.csv'
variants_DLBCL = pd.read_csv(features_path).fillna(0)
variants_DLBCL['Cohort'] = 'DLBCL' 

features_path = '/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/train_xgboost_model/ROT_features_all_chroms.csv'
variants_ROT = pd.read_csv(features_path).fillna(0) 
variants_ROT['Cohort'] = 'ROT' 

tail_scores = pd.read_csv("/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/train_xgboost_model/mobster_scores_DLBCL.csv").reset_index(drop=True)
tail_scores['POS'] = tail_scores['POS'].astype(int)
tail_scores['Variant'] = tail_scores.apply(lambda row: f"{row['CHROM']}:{row['POS']}_{row['REF']}>{row['ALT']}", axis=1)
tail_scores = tail_scores[['sample', 'Variant', 'Tail']]
variants_DLBCL = pd.merge(variants_DLBCL, tail_scores, left_on=['Sample', 'Variant'], right_on=['sample', 'Variant'], how='inner').drop('sample', axis=1)
variants_DLBCL['Label'] = variants_DLBCL['Label'].apply(lambda x: 1 if x == 'Real' else 0)

tail_scores = pd.read_csv("/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/train_xgboost_model/mobster_scores_ROT.csv").reset_index(drop=True)
tail_scores['POS'] = tail_scores['POS'].astype(int)
tail_scores['Variant'] = tail_scores.apply(lambda row: f"{row['CHROM']}:{row['POS']}_{row['REF']}>{row['ALT']}", axis=1)
tail_scores = tail_scores[['sample', 'Variant', 'Tail']]
variants_ROT = pd.merge(variants_ROT, tail_scores, left_on=['Sample', 'Variant'], right_on=['sample', 'Variant'], how='inner').drop('sample', axis=1)
variants_ROT['Label'] = variants_ROT['Label'].apply(lambda x: 1 if x == 'Real' else 0)

directory = "/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/blgsp/"
pattern = "*_extracted_features.csv"

matching_files = glob.glob(directory + pattern, recursive=True)

try: 
    features_BLGSP = pd.concat(
        (pd.read_csv(path).fillna(0) for path in matching_files), ignore_index=True)
except Exception as e:
    print(f"An error occurred: {e}")
    print(f"There is an issue with the following file of extracted features.\n{path}")

labels = pd.read_csv("/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/blgsp/labels.txt", sep='\t').fillna(0)
labels['Label'] = labels['Label'].apply(lambda x: 1 if x == 'Real' else 0)

variants_BLGSP = pd.merge(features_BLGSP, labels, on=['Sample', 'Variant'], how='inner')

directory = "/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/all_htmcp/"
matching_files = glob.glob(directory + pattern, recursive=True)
try: 
    features_HTMCP = pd.concat(
        (pd.read_csv(path).fillna(0) for path in matching_files), ignore_index=True)
except Exception as e:
    print(f"An error occurred: {e}")
    print(f"There is an issue with the following file of extracted features.\n{path}")

labels = pd.read_csv("/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/all_htmcp/labels.txt",
                    sep="\t").fillna(0)
labels['Label'] = labels['Label'].apply(lambda x: 1 if x == 'Real' else 0)

variants_HTMCP = pd.merge(features_HTMCP, labels, on=['Sample', 'Variant'], how='inner')

columns_order = variants_BLGSP.columns
variants_HTMCP = variants_HTMCP[columns_order]
variants_DLBCL = variants_DLBCL[columns_order]

variants=pd.concat([variants_BLGSP,variants_HTMCP,variants_ROT], ignore_index=True)

#Convert old hot encodings to our newer sequence context features

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

X = variants.drop(['Label'], axis=1)
y = label_encoder.fit_transform(variants['Label'])

