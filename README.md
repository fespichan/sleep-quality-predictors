# Biological and Anthropometric Determinants of Sleep Quality Across the Adult Life Course

[![R](https://img.shields.io/badge/R-4.3.3-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18902596.svg)](https://doi.org/10.5281/zenodo.18902596)

## Overview

This repository contains the analytical code for the manuscript:

> **Biological and Anthropometric Determinants of Sleep Quality Across the Adult Life Course: A Cross-Sectional Study with Multiple Imputation and Machine Learning**
>
> Espichán F., Carbajal L., Siccha Macassi A.L.
>
> *Submitted to Sleep Medicine*

The study evaluates biological and anthropometric predictors of poor sleep quality (PSQI > 5) in 289 participants aged 17–71 years from a Peruvian university community using logistic regression, random forest, and XGBoost within a multiple imputation framework (m = 100).

## Key Findings

- **Male sex** was the strongest predictor of global sleep quality (OR = 0.08; p = 0.002), driven by reduced sleep medication use (PSQI Component 6)
- **Overweight status** specifically predicted daytime dysfunction (PSQI Component 7; OR = 2.44; p = 0.046)
- Logistic regression outperformed random forest and XGBoost (AUC: 0.722 vs. 0.653 vs. 0.631; Friedman p < 0.001)
- SMOTE did not significantly improve sensitivity (p = 0.091) but reduced specificity and AUC

## Repository Structure

```
├── scripts/
│   ├── 01_Clinical_Audit.R            # Three-stage data quality validation
│   ├── 02_Multiple_Imputation.R       # MICE with random forest (m=100, maxit=50)
│   ├── 03_Reliability_Analysis.R      # Cronbach's alpha & McDonald's omega
│   ├── 04_MFA_Analysis.R             # Multiple Factor Analysis (4 configurations)
│   ├── 05_Table1_Descriptive.R       # Demographic table by sex
│   ├── 06_LR_Pooled_PSQI.R          # Logistic regression with Rubin's rules
│   ├── 07_RF_PSQI.R                  # Random forest (100 imputed datasets)
│   ├── 08_XGBoost_PSQI.R            # XGBoost (100 imputed datasets)
│   ├── 09_Component_Screening.R      # All 7 PSQI components analysis
│   ├── 10_Friedman_Comparison.R      # Model comparison (Friedman + post-hoc)
│   ├── 11_SMOTE_Justification.R      # AUC/Sensitivity/Specificity with vs without SMOTE
│   ├── 12_HL_Calibration.R          # Hosmer-Lemeshow calibration comparison
│   ├── 13_LR_Comp7_Exploratory.R    # Comp.7 (Daytime Dysfunction) LR analysis
│   └── 14_MFA_Comp7.R               # Comp.7 × BMI MFA individual factor map
│
├── figures/                           # Generated figures (main text)
├── supplementary/                     # Supplementary figures and tables
├── data/
│   └── README.md                     # Data description and access information
├── README.md
└── LICENSE
```

## Analytical Pipeline

The analysis follows this workflow:

```
Raw data
  │
  ▼
01_Clinical_Audit.R ──────────► Clean data (219 implausible values flagged)
  │
  ▼
02_Multiple_Imputation.R ─────► 100 imputed datasets (MICE, rf algorithm)
  │
  ▼
03_Reliability_Analysis.R ────► PSQI reliability (ω = 0.70, α = 0.66)
  │
  ▼
04_MFA_Analysis.R ────────────► Variable selection (adiposity + lean mass clusters)
  │
  ▼
05_Table1_Descriptive.R ──────► Table 1 (demographics by sex)
  │
  ├──► 06_LR_Pooled_PSQI.R ──► Pooled OR (Rubin's rules) + 4-panel figure
  ├──► 07_RF_PSQI.R ──────────► RF importance + PDP + ROC + confusion matrix
  └──► 08_XGBoost_PSQI.R ────► XGBoost performance metrics
        │
        ▼
09_Component_Screening.R ────► 7 PSQI components × 6 predictors (42 tests)
  │
  ▼
10_Friedman_Comparison.R ────► Model comparison (LR > RF > XGBoost)
  │
  ▼
11_SMOTE_Justification.R ───► SMOTE impact evaluation (3 metrics)
  │
  ▼
12_HL_Calibration.R ─────────► Calibration: 87% adequate without SMOTE vs 2% with SMOTE
  │
  ▼
13_LR_Comp7_Exploratory.R ──► Comp.7 LR: Overweight OR = 2.44, p = 0.046 (Fig S7)
  │
  ▼
14_MFA_Comp7.R ──────────────► Comp.7 × BMI MFA individual factor map (Fig S8)
```

## Requirements

### R version
- R ≥ 4.3.3

### Required packages

```r
# Data manipulation
install.packages(c("dplyr", "tidyr"))

# Multiple imputation
install.packages("mice")

# Modeling
install.packages(c("caret", "randomForest", "xgboost"))

# SMOTE
install.packages(c("themis", "recipes"))

# Evaluation
install.packages(c("pROC", "ResourceSelection"))

# Multivariate analysis
install.packages(c("FactoMineR", "factoextra"))

# Reliability
install.packages(c("psych"))

# Tables and export
install.packages(c("gtsummary", "writexl"))

# Visualization
install.packages(c("ggplot2", "gridExtra"))
```

## Data and Code Availability

- **Code:** All analytical scripts are archived in Zenodo ([DOI: 10.5281/zenodo.18902596](https://doi.org/10.5281/zenodo.18902596)) and available at this GitHub repository
- **Data:** Raw data are not publicly available due to participant privacy protections approved by the Institutional Ethics Committee of the Universidad Nacional del Callao. Data are available from the corresponding author upon reasonable request

**Contact:** fabio.espichan.j@upch.pe

### Data description

| Variable | Type | Description |
|----------|------|-------------|
| sex | Categorical | Female / Male |
| age | Continuous | Age in years (17–71) |
| height | Continuous | Height in cm |
| weight | Continuous | Weight in kg |
| bmi | Derived | Body mass index (kg/m²) |
| wc | Continuous | Waist circumference (cm) |
| bf | Continuous | Body fat percentage (Deurenberg equation) |
| water | Continuous | Total body water percentage (BIA) |
| mm | Continuous | Skeletal muscle mass (kg, BIA) |
| hb | Continuous | Hemoglobin concentration (g/dL) |
| Comp.1–Comp.7 | Ordinal (0–3) | PSQI component scores |
| psqi.total | Derived | PSQI global score (0–21) |
| psqi | Derived | Good (≤5) vs. Poor (>5) sleep quality |

## Reproducibility Notes

- All scripts use `set.seed(123)` for reproducibility
- Multiple imputation generates 100 datasets with 50 iterations each
- SMOTE is applied exclusively to training partitions
- Test sets preserve original class distribution
- Rubin's rules are used for logistic regression pooling
- Machine learning metrics are averaged across 100 imputed test sets

## Citation

If you use this code or methodology, please cite both the paper and the code:

### Paper
```bibtex
@article{espichan2026sleep,
  title={Biological and Anthropometric Determinants of Sleep Quality 
         Across the Adult Life Course: A Cross-Sectional Study with 
         Multiple Imputation and Machine Learning},
  author={Espich{\'a}n, Fabio and Carbajal, Luz and Siccha Macassi, Ana Lucy},
  journal={Sleep Medicine},
  year={2026},
  note={Submitted}
}
```

### Code
```bibtex
@software{espichan2026sleep_code,
  author       = {Espich{\'a}n, Fabio and Carbajal, Luz and Siccha Macassi, Ana Lucy},
  title        = {Analytical code for: Biological and Anthropometric 
                  Determinants of Sleep Quality Across the Adult Life Course},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18902596},
  url          = {https://doi.org/10.5281/zenodo.18902596}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Archiving and DOI

This repository is archived in [Zenodo](https://zenodo.org/) to ensure long-term preservation and provide a citable DOI.

**Archived version:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18902596.svg)](https://doi.org/10.5281/zenodo.18902596)

### How to cite the code

```bibtex
@software{espichan2026sleep_code,
  author       = {Espich{\'a}n, Fabio and Carbajal, Luz and Siccha Macassi, Ana Lucy},
  title        = {Analytical code for: Biological and Anthropometric 
                  Determinants of Sleep Quality Across the Adult Life Course},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18902596},
  url          = {https://doi.org/10.5281/zenodo.18902596}
}
```

### How to create the Zenodo archive (for authors)

1. Push all scripts to GitHub
2. Link your GitHub account at [zenodo.org](https://zenodo.org)
3. Go to Settings → GitHub and enable this repository
4. In GitHub, create a Release (e.g., v1.0.1) with the tag `v1.0.1`
5. Zenodo will automatically generate a DOI
6. Update the DOI badge in this README and in the manuscript

### Version history

| Version | Date | Description |
|---------|------|-------------|
| v1.0.1 | 2026-03-09 | Initial release (manuscript submission) |

## Acknowledgments

- Universidad Peruana Cayetano Heredia (UPCH)
- Universidad Nacional del Callao
- All study participants
