# Scripts and Notebooks

This folder contains Jupyter notebooks used for development, testing, exploratory 
analysis, and model training for the NMD-Scanner project. These notebooks are 
**not part of the Python package** but document the code used during the research and 
implementation process.

## Contents

### 1. `create_test_VCF.ipynb`
Creates synthetic VCF files containing possible edge-case variants used to test the robustness 
of the function incorporating variants into the sequence.

The notebook generates variants such as:
- start codon substitutions, deletions, insertions, duplications  
- stop codon substitutions, length-changing variants at the end of the CDS  
- nonsense mutations  
- frameshift and in-frame indels  
- variants before and after the CDS  
- complex or multi-base substitutions

---

### 2. `train.ipynb`
Notebook used to train predictive models to evaluate NMD efficiency. The workflow includes:

#### **(A) Benchmarking against NMDEff**
- Runs the original NMDEff implementation (<https://github.com/hjkng/nmdeff>)  
- Extracts NMD efficiency scores using their official scripts  
- Prepares the same variants using the NMD-Scanner to compare feature sets and predictions  

#### **(B) Model exploration**
Multiple machine-learning models were tested, including:
- Random Forest
- XGBoost
- LightGBM
- Gradient Boosting
- Ridge
- Lasso
- SVR 

#### **(C) Final model selection**
- After evaluation, the Random Forest Regressor gave the best performance based on cross-validation metrics.
- The trained model is saved (joblib) for later prediction.

---

### 3. `validation_MMRF_TARGET.ipynb`
Evaluation of the trained model using independent datasets (MMRF/TARGET, <https://github.com/hjkng/nmdeff>).
It includes:
- loading validation datasets  
- generating NMD-related features  
- comparing predictions with reference measurements  
- producing evaluation plots and performance metrics  

---

### 4. `NMD.ipynb` (WILL BE DEPRECATED / REPLACED)
This notebook originally contained exploratory code for implementing the NMD logic.  
It can be deleted or replaced with a cleaner notebook containing:
- examples of how to use the installed `nmd_scanner` package  
- code snippets to manually inspect transcripts, PTCs, and escape rules  
- debugging notes

---

## Notes
- None of the notebooks are required for end users of the package.  
- They are included for transparency and reproducibility of the development process.
