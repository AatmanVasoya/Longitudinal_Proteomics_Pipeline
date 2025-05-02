# 🧪 Proteomics Preprocessing Reference Code

This repository contains a standardized preprocessing pipeline for longitudinal proteomics data, developed as part of a biomarker discovery project in Alzheimer's disease. The script aligns patient data to disease onset, corrects labels, and generates slope-ready dataframes for downstream modeling.

---

## 📌 What This Script Does

The script:
- Cleans column names for consistency
- Adjusts onset age (`ONSET_AGE`) based on disease progression or DECAGE
- Flags patients who transitioned from Normal to Abnormal
- Corrects misaligned status labels across visits
- Creates piecewise variables for slope modeling
- Returns filtered DataFrames for:
  - All patients with status tracking
  - Transition patients
  - Normal-only and Abnormal-only cohorts

---

## 🗂 Required Columns in Input CSV

Your input file (e.g., `data_with_onset.csv`) **must contain** the following columns:

| Column Name       | Description |
|-------------------|-------------|
| `SUBID`           | Unique patient identifier (repeats across visits) |
| `PROCEDURE_DATE`  | Date of visit/sample collection |
| `BASELINE_AGE`    | Age at study entry (optional but retained) |
| `DECAGE`          | Age at last clinical assessment or visit |
| `PROCEDURE_AGE`   | Age at current visit |
| `SEX`             | Biological sex of the participant |
| `CATEGORY`        | Original diagnosis label (e.g., Normal, MCI, AD) |
| `Status_H`        | Harmonized status (must be either 'Normal' or 'Abnormal') |
| `TRANSITION_AGE`  | Age when disease onset is believed to have occurred |
| `ONSET_AGE`       | Onset age, if available directly |
| `Status`          | (Optional) Original status field |
| `Protein_1`, `Protein_2`, ..., `Protein_n` | Proteomic measurement columns (ideally starting from column 15) |

**⚠️ Note**: Proteomic columns should begin from column 15 onward. These are all numeric columns representing protein expression values.

---

## 💻 How to Use

1. Clone or download this repo
2. Make sure your CSV file follows the format above
3. Import the preprocessing function in your analysis script:

```python
from reference_preprocessing.preprocess_proteomics_data import preprocess_data

# Load and preprocess the data
results = preprocess_data("data_with_onset.csv")

# Access the returned DataFrames
working_df_all = results["working_df_all"]         # All patients
working_df = results["working_df"]                 # Transition patients only
df_normal_only = results["df_normal_only"]         # Patients always Normal
df_abnormal_only = results["df_abnormal_only"]     # Patients always Abnormal
