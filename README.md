# FIFA for FFPE Filtering

README last updated March 25, 2026

## Installation

Create environment
```bash
conda env create -f conda/FIFA_environment.yml
```

Install MOBSTER for R
```R
> install.packages("devtools")
> devtools::install_github("caravagnalab/mobster")
```

###
To access the help manual, run:

```bash
/path/to/fifa/src/cli.py --help
```

*Note: For convenience, the user can create a symlink between `/path/to/fifa/src/cli.py` and the command `fifa`*

---

### The tool can be run on three different modes: 
#### 1) Extract
For each SNV specified in a VCF file, extract features using the sample's VCF and BAM file. Creates a table of features that can be used to make predictions, or to train a new EBM model. 

```bash
/path/to/fifa/src/cli.py extract -s [SAMPLE NAME] -c [COHORT] -v [VCF PATH] -b [BAM PATH] -r [REF SEQ] -o [OUTPUT DIR] -n [NUMBER OF THREADS ALLOCATED]
```

_Optional Parameters:_ 

- Variants can be labeled with a user-specified flag. Use `-l` for this purpose:
  ```
  -l [INFO FIELD WITH VARIANT LABEL] [TRUTH VALUE]  (default: None, None)
  -l 'INFO/Label' 'Real'
  ```
  
- If intending to run the original feature extraction script, please include the flag "-p". The original script does not properly paralelize processes across resources, so in practice there's no real reason to use it other than for testing resource allocation / resource efficiency. 
  ```
  -p  # Include this flag to enable the original scheme.
  ```
---

#### 2) Predict

```bash
/path/to/fifa/src/cli.py predict -s [SAMPLE NAME] -v [VCF PATH] -f [FEATURES PATH] -o [OUTPUT DIR FOR VCF WITH ANNOTATIONS] -m [PATHS TO EBM MODELS]
```

#### Notes on Model Paths:
If multiple models are submitted, they will automatically be merged together. Current FIFA models are:

- **NYGC1**: `/path/to/fifa/models/ebm_hyperparams_NYGC1.pkl`
- **NYGC2**: `/path/to/fifa/models/ebm_hyperparams_NYGC2.pkl`
- **HTMCP**: `/path/to/fifa/models/ebm_hyperparams_CGCI-HTMCP.pkl`
- **BLGSP**: `/path/to/fifa/models/ebm_hyperparams_CGCI-BLGSP.pkl`

To merge one of the core FIFA models, the user can provide either the path or simply specify the cohort by name [NYGC1, NYGC2, HTMCP, BLGSP]

#### Optional Parameters:
- Add annotations from the RNA pileups to the output VCF. Will automatically rescue variants erroneously discarded by FIFA:
  ```bash
  -r [PATH TO RNA PILEUPS]
  ```
---

#### 3) Retrain
Train a new EBM model using a new cohort
*Note : Data must be labeled.*

```bash
/path/to/fifa/src/cli.py retrain -o [OUTPUT PATH] -d [DIRECTORY CONTAINING NEW SAMPLES] -l [CSV FILE WITH TRUE LABELS]
```

_Note:_ 
The program assumes that the directory contains a series of files in the format "[SAMPLE]_extracted_features.csv". For example, see :
```
/path/to/fifa/data/1395_tumor_ffpe_wgs_extracted_features.csv
```

_Optional Parameters:_ 
  ```bash
  -hp  # FLAG; CONDUCT HYPERPARAMETER GRID-SEARCH WHEN TRAINING THE NEW EBM
  ```

##### 4) Merge Models

Additional script for merging models together, and saving them as a new model. Produces preliminary metrics on the HCC1395 chr1 variants. 
```bash
/path/to/fifa/src/cli.py merge -o [OUTPUT PATH FOR FIFA MODEL] -m [PATHS TO INPUT MODELS] 
```
