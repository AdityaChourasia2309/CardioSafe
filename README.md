# Predicting hERG Channel Blockade: A Multi-Modal Framework for Early Cardiac Safety Screening

## Overview

This project predicts whether a drug compound will block the hERG heart channel (a major cause of cardiac side effects and drug failures [1]) based on its chemical structure and its structure variants. The model outputs a blocker probability and a predicted potency value (pIC50) for each compound, and extends further toward assessing binding affinities across hERG variants [2].

## The "Why": Biological & Business Context

**The Stakes:** hERG blockade is a predominant cause of late-stage drug attrition and sudden cardiac death risk [1,4].

**The Biological Mechanism:** The hERG channel controls the potassium current during Phase 3 cardiac repolarization [1,4]. Blockade leads to QT interval prolongation and can result in lethal Torsades de Pointes (TdP) [1].

**Regulatory Context:** ICH S7B regulatory requirements mandate early in silico cardiac screening prior to IND submission (per EMA guidelines).

**The Bottleneck:** Standard ligand-only QSAR models fail when encountering novel chemical space or variation in target conformations/mutants — a gap this project directly investigates.

---

## Repository Structure

```text
herg-blocker-prediction/
├── README.md
├── requirements.txt
├── submission.csv
├── model_codes/
│   ├── 01_fetch_data.ipynb
│   ├── 02_curate.ipynb
│   ├── 03_split_featurize.ipynb
│   ├── 04a_baseline_models.ipynb
│   ├── 04b_baseline_models.ipynb
│   └── 05_applicability_domain.ipynb
├── checkpoints/
│   ├── model_logreg.pkl
│   ├── model_rf_classifier.pkl
│   ├── model_rf_regressor.pkl
│   ├── model_rf_ligand_protein.pkl
│   └── model_xgboost.pkl
├── reports/
│   ├── data_curation_report.md
│   ├── exploration_report.md
│   ├── herg_variant_subset.csv
│   ├── data_summary_plots.png
│   ├── shap_summary_plot.png
│   ├── shap_bar_plot.png
│   └── AD_accuracy_by_similarity.png
├── data/
│   ├── raw/
│   │   └── raw_chembl_herg.csv
│   ├── herg_clean.csv
│   ├── herg_split.csv
│   ├── herg_featurized_meta.csv
│   ├── descriptors.csv
│   ├── herg_protein_features.csv
│   ├── herg_protein_Q12809.fasta
│   ├── test_with_AD_scores.csv
    └── fp_matrix.npy
```

### Notebooks

- **`01_fetch_data.ipynb`** — Pulled real experimental data from the ChEMBL public bioactivity database, specifically hERG target CHEMBL240. Extracted 20,073 raw bioactivity records covering 16,215 unique compounds.

- **`02_curate.ipynb`** — Cleaned the raw data by standardizing molecule representations (SMILES canonicalization) and removing salts/junk. Converted IC50 values into pIC50 and labeled each compound as a blocker/non-blocker using a threshold of pIC50 ≥ 5.0.

- **`03_split_featurize.ipynb`** — Split the data into 80% training and 20% test using a scaffold split to ensure genuinely novel molecule shapes in the test set. Also ran a diagnostic random-split comparison. Converted each molecule into 2048-bit chemical fingerprints (ECFP4) plus 8 physical/chemical properties.

- **`04a_baseline_models.ipynb` & `04b_baseline_models.ipynb`** — Trained four models: Logistic Regression, Random Forest, Random Forest with added protein features, and XGBoost — all evaluated on the same held-out test set. Also trained a separate regression model to predict potency (pIC50). Ran SHAP interpretability analysis to determine why the model makes its predictions (top drivers: TPSA, LogP, hydrogen bond acceptors — matching known hERG pharmacology). Pulled hERG's real protein sequence and added protein-derived features to test ligand+protein integration. Checked the dataset for known hERG mutant variants — found 3 of 4, but only 8 compounds, too few to model separately. Generated `submission.csv` with all predictions.

- **`05_applicability_domain.ipynb`** — Built a confidence-scoring system measuring how chemically similar each test compound is to the training data, allowing the model to flag low-confidence predictions.

### Data & Outputs

- `raw_chembl_herg.csv` — Raw dataset extracted from the database
- `herg_clean.csv` — Curated dataset containing 16,159 unique compounds (9,000 blockers, 7,159 non-blockers)
- `herg_split.csv` — Hold-out set data reflecting the 80/20 train/test split
- `descriptors.csv` — Extracted physicochemical features for the compounds
- `fp_matrix.npy` — Fingerprint matrix (ECFP4) saved as a NumPy array
- `herg_protein_features.csv` — Protein structural and sequence feature data
- `herg_protein_Q12809.fasta` — FASTA file containing the hERG target sequence
- `herg_variant_subset.csv` — Subset data covering varying hERG structural configurations
- `herg_featurized_meta.csv` — Processed metadata tracking the featurization pipeline
- `data_curation_report.md` — Documentation detailing data cleaning and curation steps
- `submission.csv` — Final predicted blockade probabilities and outputs
- `test_with_AD_scores.csv` — Test set predictions appended with Applicability Domain scores

**Note:** the associated assay types in the source data are — B-Binding (73.27%), T-Binding (22.16%), F-Binding (3.48%), and A-ADME (1.10%).

### Models

- `model_logreg.pkl` — Logistic Regression baseline model
- `model_rf_classifier.pkl` — Random Forest classification model
- `model_rf_regressor.pkl` — Random Forest regressor (pIC50)
- `model_rf_ligand_protein.pkl` — Random Forest model utilizing dual-stream ligand-protein features
- `model_xgboost.pkl` — XGBoost classification model

### Visualizations

- `data_summary_plots.png` — Exploratory data analysis distributions
- `shap_summary_plot.png` — Model interpretability summary plot highlighting overall feature importance
- `shap_bar_plot.png` — Bar plot ranking the top SHAP value features
- `AD_accuracy_by_similarity.png` — Model performance mapped against Tanimoto similarity to the training set

---

## Environment Setup

Create and activate a Python 3.10+ environment, then install dependencies:

```bash
pip install -r requirements.txt
```

Core libraries include `rdkit`, `scikit-learn`, `xgboost`, `shap`, `pandas`, `numpy`, and `matplotlib`.

**Known issue:** some environments (particularly fresh Google Colab sessions) hit a numpy/pandas binary incompatibility error on install:

ValueError: numpy.dtype size changed, may indicate binary incompatibility

If this occurs, run:
```bash
pip install --upgrade --force-reinstall numpy pandas
pip install --upgrade --force-reinstall --no-deps rdkit
```
then **restart the Python kernel/runtime** before importing anything further.

---

## How To Run

Run the following notebooks from `model_codes/`, in order, from the repository root. Each notebook depends on outputs from the previous one.

### 1. Fetch data

model_codes/01_fetch_data.ipynb

Requires internet access. Takes ~3–5 minutes.

Expected outputs:
- `data/raw/raw_chembl_herg.csv`

### 2. Curate and clean

model_codes/02_curate.ipynb

Canonicalizes SMILES via RDKit, converts IC50 → pIC50, derives binary blocker labels (pIC50 ≥ 5.0), deduplicates on canonical SMILES.

Expected outputs:
- `data/herg_clean.csv`
- `reports/data_curation_report.md`
- `reports/data_summary_plots.png`

The script also prints explicit counts at every filtering step: raw records in, records dropped for invalid SMILES, records after deduplication, and final class balance — confirming no silent data loss.

### 3. Split and featurize

model_codes/03_split_featurize.ipynb

Creates the primary scaffold-based 80/20 train/test split (Murcko scaffolds), generates ECFP4 fingerprints (2048-bit) and 8 physicochemical descriptors. Also runs a diagnostic random-split comparison to quantify how much easier random splitting makes evaluation. Takes ~2–3 minutes.

Expected outputs:
- `data/herg_split.csv`
- `data/fp_matrix.npy`
- `data/descriptors.csv`
- `data/herg_featurized_meta.csv`

### 4. Train baseline models

model_codes/04a_baseline_models.ipynb
model_codes/04b_baseline_models.ipynb

Trains Logistic Regression, Random Forest (classifier + regressor), and XGBoost. Runs SHAP interpretability analysis. Integrates hERG protein sequence (UniProt Q12809) features alongside ligand features. Checks raw data for known hERG mutant variants. Takes ~10–15 minutes.

Expected outputs:
- `checkpoints/model_logreg.pkl`
- `checkpoints/model_rf_classifier.pkl`
- `checkpoints/model_rf_regressor.pkl`
- `checkpoints/model_rf_ligand_protein.pkl`
- `checkpoints/model_xgboost.pkl`
- `submission.csv`
- `reports/shap_summary_plot.png`
- `reports/shap_bar_plot.png`
- `reports/herg_variant_subset.csv`

### 5. Applicability domain / confidence scoring

model_codes/05_applicability_domain.ipynb

Computes each test compound's maximum Tanimoto similarity to the training set, validates it as a meaningful confidence signal by comparing in-domain vs out-of-domain accuracy. Takes ~2–5 minutes.

Expected outputs:
- `data/test_with_AD_scores.csv`
- `reports/AD_accuracy_by_similarity.png`
- Updated `submission.csv` with `AD_similarity_score` and `high_confidence` columns

### Required Execution Order

model_codes/01_fetch_data.ipynb
model_codes/02_curate.ipynb
model_codes/03_split_featurize.ipynb
model_codes/04a_baseline_models.ipynb
model_codes/04b_baseline_models.ipynb
model_codes/05_applicability_domain.ipynb


**Configuration:** Each notebook defines a working directory via the `SAVE_DIR` variable near the top of the file (e.g. `SAVE_DIR = '/content/drive/MyDrive/herg_hackathon'`), currently configured for the authors' Google Drive environment. Update this single variable to your own path before running each notebook.

Review `reports/exploration_report.md` for the full methodology writeup and discussion of results.

### Skipping Retraining (Inference Only)

To use pre-trained checkpoints directly instead of re-running the full pipeline:
```python
import joblib
model = joblib.load("checkpoints/model_rf_classifier.pkl")
model.predict_proba(X)  # X must be generated the same way as in 03_split_featurize.ipynb
```

---

## Key Metrics & Results

**Dataset Overview:**
- Total unique compounds: 16,159
- Blockers: 9,000 (55.7%)
- Non-blockers: 7,159 (44.3%)
- Train set: 12,866
- Test set (held-out): 3,293

### Classification (Blocker vs. Non-Blocker)
*Evaluated on the held-out scaffold-split test set (n=3,293)*

| Model | Accuracy | Precision | Recall | F1 | ROC-AUC | PR-AUC |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Logistic Regression | 0.734 | 0.759 | 0.782 | 0.770 | 0.810 | 0.847 |
| **Random Forest (best)** | **0.777** | **0.788** | **0.831** | **0.809** | **0.865** | **0.897** |
| RF + Protein Sequence | 0.776 | 0.784 | 0.839 | 0.810 | 0.866 | 0.896 |
| XGBoost | 0.768 | 0.784 | 0.818 | 0.801 | 0.854 | 0.883 |

### Regression (pIC50 Potency Prediction)
**Random Forest Regressor:** RMSE 0.688 | MAE 0.471 | R² 0.486

### Methodology Validation: Scaffold vs. Random Split
- **Random split (easier, inflated):** Accuracy 0.814, ROC-AUC 0.900
- **Scaffold split (honest evaluation):** Accuracy 0.777, ROC-AUC 0.865
- *Diagnostic Conclusion:* the scaffold split is genuinely harder — random splitting inflated accuracy by ~4 points.

### Applicability Domain (Confidence Validation)
- **In-domain (high confidence, n=2,943):** Accuracy 0.790, ROC-AUC 0.876
- **Out-of-domain (low confidence, n=350):** Accuracy 0.660, ROC-AUC 0.736
- *Conclusion:* compounds similar to training data were predicted with 79% accuracy, while structurally novel compounds dropped to 66% accuracy — a validated, actionable confidence signal.

### SHAP Interpretability
- **Top drivers:** TPSA, LogP, and hydrogen bond acceptors
- **Insight:** these drivers match established hERG pharmacology, confirming the model learned real biological structure-activity relationships rather than spurious correlations.

---

## Report Compliance Notes

- **Sanity checks:** `02_curate.ipynb` prints explicit counts at every curation step (raw records, dropped invalid SMILES, post-dedup compound count, final class balance), confirming no silent data loss.
- **Train/test independence:** verified via scaffold split (Murcko scaffolds); a random-split diagnostic in `03_split_featurize.ipynb` confirms the scaffold split is the harder, more honest evaluation.
- **SMILES canonicalization:** performed in `02_curate.ipynb` via RDKit, with salt stripping and largest-fragment retention.
- **Interpretability:** SHAP analysis in `04a_baseline_models.ipynb` identifies TPSA, LogP, and HBA as the top predictive features, consistent with established hERG pharmacology.
- **Variant analysis:** raw ChEMBL data was searched for the 4 classic hERG mutant variants (Y652A, F656A, S624A, N588K); 3 were found, covering only 8 unique compounds — documented in `reports/herg_variant_subset.csv`.

---

## Large File / Archive Note

Due to GitHub's 100MB per-file limit, `data/fp_matrix.npy` and three Random Forest checkpoints are excluded from this repository via `.gitignore`. All are fully reproducible by re-running the corresponding notebook, and are also available in the complete project archive: **https://drive.google.com/drive/folders/1Jk5GhyXGPqdNQ1NYKq01I17xM54iPrwR?usp=sharing**

---

## Future Theoretical Framework

To move past the bottleneck of current ligand-only QSAR models, this project lays the groundwork for a generalized model architecture. The goal is to build an integrated model capable of predicting compound toxicity against varying hERG structural configurations, rather than limiting predictions to a single (wild-type) protein target. Proposed next steps include:

- A GNN with attention-based interpretability, visualizing which atoms/substructures drive predictions
- ProtBERT/ChemBERT embeddings in place of hand-crafted protein composition features
- Expanded variant-resolved data collection to support true multi-variant modeling
- Molecular docking for structural binding-pose evidence
- Conformal prediction to formalize the current applicability-domain flag into calibrated confidence intervals

---

## Acknowledgements & References

**Generative AI Tools Utilized:**
- Google Gemini
- Google AI Mode
- Anthropic's Claude AI

**Bibliography:**
1. Vandenberg, J. I., Perry, M. D., Perrin, M. J., Mann, S. A., Ke, Y., & Hill, A. P. (2012). hERG K(+) channels: structure, function, and clinical significance. Physiological reviews, 92(3), 1393–1478. https://doi.org/10.1152/physrev.00036.2011
2. Syahdi, R. R.; Jasial, S.; Maeda, I.; Miyao, T. Bridging Structure- and Ligand-Based Virtual Screening through Fragmented Interaction Fingerprint. *ACS Omega* 2024, 9 (37), 38957–38969. https://doi.org/10.1021/acsomega.4c05433.
3. Siramshetty, V. B.; Nguyen, D.-T.; Martinez, N. J.; Southall, N. T.; Simeonov, A.; Zakharov, A. V. Critical Assessment of Artificial Intelligence Methods for Prediction of hERG Channel Inhibition in the "Big Data" Era. *J. Chem. Inf. Model.* 2020, 60 (12), 6007–6019. https://doi.org/10.1021/acs.jcim.0c00884.
4. Kalyaanamoorthy, S.; Barakat, K. H. Development of Safe Drugs: The hERG Challenge. *Medicinal Research Reviews* 2018, 38 (2), 525–555. https://doi.org/10.1002/med.21445.
5. Shishir, F. S.; Sarker, B.; Zhong, C.; Shomaji, S. Parameter Efficient Deep Learning Models for Multi-Target Binding Affinity and hERG Cardiotoxicity Prediction. *IEEE Trans. Comput. Biol. Bioinform.* 2026, 1–14. https://doi.org/10.1109/TCBBIO.2026.3705951.
6. Sanches, I. H.; Braga, R. C.; Alves, V. M.; Andrade, C. H. Enhancing hERG Risk Assessment with Interpretable Classificatory and Regression Models. *Chem. Res. Toxicol.* 2024, 37, 910–922. https://doi.org/10.1021/acs.chemrestox.3c00400.
