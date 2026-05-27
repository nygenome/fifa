# FIFA for FFPE Filtering

Workflow from:

**An explainable boosting machine model for identifying artifacts caused by formalin-fixed paraffin embedding**

Valentina Grether, Zoe R. Goldstein, Jennifer M. Shelton, Timothy R. Chu, William F. Hooper, Heather Geiger, André Corvelo, Rachel Martini, Melissa B. Davis, Nicolas Robine, Will Liao

bioRxiv 2026.03.10.710815; doi: https://doi.org/10.64898/2026.03.10.710815

## Installation

### Install python3 modules
```bash
$ python3 -m \
  pip install \
  --no-cache-dir --upgrade \
  "pip<25" \
  "setuptools<82" \
  wheel \
  pybind11 \
  "Cython<3" \
  "numpy==1.24.4" \
  "pysam==0.15.4"

$ python3 \
-c "import pkg_resources, pybind11, numpy, pysam"

$ python3 -m \
  pip install \
  --no-build-isolation \
  --no-cache-dir \
  -r requirements.txt
```

### Install R packages 

Installation tested with R v 4.0.3.

```R
> Rscript install_packages.R
```

### Set fifa alias

```bash
echo "alias fifa='python3 PATH_TO_REPO/fifa/src/cli.py'"
```

### Docker

Fifa is available with its required modules and packages as a prebuilt docker image. The DockerFile is distributed in this repo.

```
us.gcr.io/nygc-comp-s-fd4e/fifa:0.1
```



## USAGE
To access the help manual, run:

```bash
fifa --help
```

*Note: For convenience, the user can create a symlink between `/path/to/fifa/src/cli.py` and the command `fifa`*

---

### Subcommands

Fifa can be run on three different modes. WDL workflows are available to run fifa.

```
prediction_wkf.wdl
merging_wkf.wdl
# retraining_wkf.wdl
```

#### 1) Extract
For each SNV specified in a VCF file, extract features using the sample's VCF and BAM file. Creates a table of features that can be used to make predictions, or to train a new EBM model. 

```bash
fifa extract --help
usage: fifa extract [-h] -s SAMPLE -v VCFFILE -b BAMFILE -r REFSEQ [-c COHORT] [-n NUM_THREADS] [-o OUTPUT_PATH] [-l LABEL LABEL] [-p]

optional arguments:
  -h, --help            show this help message and exit
  -s SAMPLE, --sample SAMPLE
                        Name of sample from VCF being submitted to be processed.
  -v VCFFILE, --vcffile VCFFILE
                        Path to sample's VCF file, with variants to extract features for.
  -b BAMFILE, --bamfile BAMFILE
                        Path to sample's BAM file, with reads to extract features for.
  -r REFSEQ, --refseq REFSEQ
                        Path to reference file used for alignment of the sample's BAM file.
  -c COHORT, --cohort COHORT
                        Sample's Cohort
  -n NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads allocated
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Path to save output files with variants and all extracted features. (default: current directory)
  -l LABEL LABEL, --label LABEL LABEL
                        VCF INFO Field with Variant Label and Real/Truth/1 Label. (default: None, None)
  -p , --original_parallel
                        Boolean: Whether we're using the original (true) or new script (false)
```

_Optional Parameters:_ 

- Variants can be labeled with a user-specified flag. Use `-l` for this purpose:
  ```
  -l [INFO FIELD WITH VARIANT LABEL] [TRUTH VALUE]  (default: None, None)
  -l 'INFO/Label' 'Real'
  ```
  
- If intending to run the original feature extraction script, please include the flag `-p`. The original script does not properly paralelize processes across resources, so in practice there's no real reason to use it other than for testing resource allocation / resource efficiency. 
  ```
  -p  # Include this flag to enable the original scheme.
  ```
---

#### 2) Predict

```bash
fifa predict --help
usage: fifa predict [-h] -s SAMPLE -v VCFFILE -f FEATURES_PATH [FEATURES_PATH ...] [-o OUTPUT_DIR] -m MODEL_FILES [MODEL_FILES ...]
                      [-r RNA_ANNOTATIONS]

optional arguments:
  -h, --help            show this help message and exit
  -s SAMPLE, --sample SAMPLE
                        Sample Name For submitting individual sample for processing.
  -v VCFFILE, --vcffile VCFFILE
                        Path to sample's VCF file, with variants to classify.
  -f FEATURES_PATH [FEATURES_PATH ...], --features_path FEATURES_PATH [FEATURES_PATH ...]
                        Paths to CSV files with features for each variant
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to save VCFs with model predictions. (default: current directory)
  -m MODEL_FILES [MODEL_FILES ...], --model_files MODEL_FILES [MODEL_FILES ...]
                        Path file(s) of EBM models to be used. (If submitting multiple models, they will be merged.)
  -r RNA_ANNOTATIONS, --rna_annotations RNA_ANNOTATIONS
                        Path to file with RNA annotations
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
fifa retrain --help
usage: fifa retrain [-h] [-o OUTPUT_PATH] [-d DIRECTORY] [-l LABELS_PATH] [-hp]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Path to save new model. (default: "fifa_model.pkl" in current directory)
  -d DIRECTORY, --directory DIRECTORY
                        Path to directory containing all samples for training a new cohort's model. Each sample file should be a CSV with
                        extracted features for all variants in the sample.
  -l LABELS_PATH, --labels_path LABELS_PATH
                        Path to CSV file with true labels for all variants (True Labels should be: Real or 1)
  -hp , --hyperparameter
                        Boolean: Whether to conduct hyperparameter grid-search when training new EBM model
```

_Note:_ 
The program assumes that the directory contains a series of files in the format `[SAMPLE]_extracted_features.csv`. For example, see :
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
fifa merge --help
usage: fifa merge [-h] [-o OUTPUT_PATH] -m INPUT_MODELS [INPUT_MODELS ...]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Path to save merged FIFA model. If not provided, the merged model is not saved.
  -m INPUT_MODELS [INPUT_MODELS ...], --input_models INPUT_MODELS [INPUT_MODELS ...]
                        Path file(s) of EBM models to be merged. If intending to use original FIFA models, please specify: NYGC1, NYGC2,
                        HTMCP, or BLGSP
```
