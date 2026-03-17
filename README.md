# FIFA for FFPE Filtering

README last updated Dec. 16, 2025

To access the help manual, run:

```bash
/path/to/fifa/src/cli.py help
```

*Note: Ideally, the repo would automatically create symlink can be created between `/path/to/fifa/src/cli.py` and the command `fifa` for more convenience. But baby steps*

---

### The tool can be run on three different modes: 
#### 1) Extract
For each SNV specified in a VCF file, extract features using the sample's VCF and BAM file. Creates a table of features that can be used to make predictions, or to train a new EBM model. 

```bash
/path/to/fifa/src/cli.py extract -s [SAMPLE NAME] [COHORT] [VCF PATH] [BAM PATH] [REF SEQ] -o [OUTPUT DIR] -n [NUMBER OF THREADS ALLOCATED]
```

_Optional Parameters:_ 

- Variants can be labeled with a user-specified flag. Use `-l` for this purpose:
  ```
  -l [INFO FIELD WITH VARIANT LABEL] [TRUTH VALUE]  (default: None, None)
  -l 'INFO/Label' 'Real'
  ```
  
- If intending to run my original paralelization scheme (from before my talk with Andre) include the flag "-p". In practice there's no real reason to do this other than for testing resource allocation / resource efficiency. 
  ```
  -p  # Include this flag to enable the original scheme.
  ```
---

#### 2) Predict

```bash
/path/to/fifa/src/cli.py predict -s [SAMPLE NAME] [VCF PATH] -f [FEATURES PATH] -o [OUTPUT DIR FOR VCF WITH ANNOTATIONS] -m [PATHS TO EBM MODELS]
```

#### Notes on Model Paths:
If multiple models are submitted, they will automatically be merged together. Current FIFA models are:

- **DLBCL**: `/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/models/ebm_hyperparams_DLBCL.pkl`
- **ROT**: `/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/models/ebm_hyperparams_ROT.pkl`
- **HTMCP**: `/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/models/ebm_hyperparams_CGCI-HTMCP.pkl`
- **BLGSP**: `/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/models/ebm_hyperparams_CGCI-BLGSP.pkl`

(In the future I'll upload these to the repo)

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
/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/blgsp/extracted_features
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
