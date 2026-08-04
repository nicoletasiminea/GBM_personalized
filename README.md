# 🧬 Personalized Drug Discovery for Glioblastoma

![GitHub repo size](https://img.shields.io/github/repo-size/your-username/your-repo)
![GitHub last commit](https://img.shields.io/github/last-commit/your-username/your-repo)
![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Status](https://img.shields.io/badge/status-active-success.svg)

---

## 📖 Overview

This repository focuses on identifying **personalized treatment strategies** for patients with **glioblastoma** using:
- Differential gene expression analysis  
- Network-based modeling  
- Drug-target interaction analysis  

---

## 🗂️ Repository Structure

<summary><b>📁 R_code_data</b></summary>

- Contains datasets for **Differentially Expressed Gene (DEG)** analysis  
- Comparisons:
  - primary tumors vs healthy  
  - recurrent tumors vs healthy  
- The data directory contains:

A table with the RNA-seq analysis results from the TCGA-GBM project: https://portal.gdc.cancer.gov/projects/TCGA-GBM
Two tables defining the comparison groups: one for primary vs. normal and another for recurrent vs. normal.
The file dif_analysis.R, which contains a customized analysis of RNA-seq results at the individual sample level, and test_on_sample_Z.R,
which contains a Z-score-based variant analysis of the RNA-seq results at the individual sample level.

⚡ First step of the pipeline  

---


<summary><b>📁 netcontrol_files</b></summary>

- Scripts for:
  - Network generation  
  - Network analysis  

⚠️ Requires external tool:  
👉 [https://github.com/Vilksar/NetControl4BioMed] 


---


<summary><b>📁 jupyter_code_data</b></summary>

### 📂 code/
- `utils.py` → Shared helper functions  
- Numbered notebooks → Must be executed sequentially  

⚠️ Important:  
Between **Step 1 → Step 2**, run `netcontrol_files` scripts.

---


## 📊 **Step-by-Step Notebooks**
 <summary><b>1️⃣ Data Preparation </b></summary>
&nbsp;
 
 _1_create_sets_for_networks.ipynb_

Generates datasets required for NetControl


<summary><b>2️⃣ Network Statistics</b></summary>
&nbsp;

_2_network_statistics.ipynb_

Describes generated networks


 <summary><b>3️⃣ Centrality Analysis</b></summary>
&nbsp;

_3_centralities.ipynb_ 

Computes centrality metrics

Identifies hub nodes

_3b_centralities_plots.ipynb_

Plot pie charts, histograms, box plots for the results

 <summary><b>4️⃣ Common Graph Analysis</b></summary>
&nbsp;

_4_graph_differences.ipynb_

Builds shared network across cases

Identifies:

 - Common pathways

 - Case-specific nodes

 <summary><b>5️⃣ Modularity & Communities</b></summary>
&nbsp;

_5_modularity.ipynb_

Detects communities

Identifies:

 - Community-specific genes

 - Additional pathways

 <summary><b>6️⃣ Results Processed</b></summary>
&nbsp;

_6_our_results_processed.ipynb_

Extracts:

 - Candidate drugs

 - Target proteins

 - Controllable proteins

 <summary><b>7️⃣ Validation and Enrichment Analysis</b></summary>
&nbsp;

_7_results_statistics_and_validation.ipynb_

Uses Enrichr (L1000 dataset)

Compares common drug effects

 <summary><b>8️⃣ Individual vs Generic Comparison</b></summary>

&nbsp;

_8_individual_vs_general.ipynb_

Compares individual vs generic treatments

Across:

- L1000

 - NIBR
   

_8b_primary_generic_in_individual.ipynb_

Pie chart with our generic results divided according to the number of cases

 <summary><b>9️⃣ Case Comparison</b></summary>
&nbsp;

_9_primary_recurrent_cases.ipynb_


Compares primary vs recurrent tumors:

 - Drugs

 - Target nodes

 - Controlled proteins

 - Pathways



   ## 📊 **Tests**

-_test_1_compare_ensembl.ipynb_


  - Compares the results from the customized analysis to those from the Z-score analysis.
  

 -_test_2_sider_database_prediction.ipynb_

  
  -Finds the SIDER side effects associated with the drugs identified by the network 
  controllability analysis. It relies on the ../data/meddra_all_se.tsv and ../data/drug_names.tsv files, which were downloaded from the SIDER database (https://sideeffects.embl.de/)


 -_test_3_results_on_off_target.ipynb_


  - Finds which drug targets are pharmacologically active among those identified for the selected drugs. It relies on the ../data/ph_active.csv and ../data/drugbank_vocabulary.csv files,
 which were downloaded from DrugBank (https://go.drugbank.com/).


  -_test_4_side_effects_other_targets.ipynb_
  

  - Finds the drugs for which all drug targets are present in the results.


 -_test_5_BBB_permeability.ipynb_

  - Calculates molecular descriptors for each drug identified as a solution in at least one case.
  - Proposes a passive diffusion score that has not been experimentally tested.
  

 -_test_6_interactions_after_exclusion_database.ipynb_
    
  
  - Returns several lists of interactions used in the analysis after excluding a selected database.
  It relies on protein information and direct protein–protein interaction data collected from SIGNOR, STRING, InnateDB, KEGG, and OmniPath.

  
 -_test_7_variants_processed.ipynb_

  
  - Includes sensitivity tests evaluating the effect of varying the gap between proteins, excluding individual databases, and modifying other analysis parameters.
---
Requirements
```
pandas 3.0.3
numpy 1.26.4
matplotlib 3.11.1
seaborn 0.13.2
requests 2.34.2
scipy 1.17.1
statsmodels 0.14.6
scikit-learn 1.9.0
openpyxl 3.1.5
rdkit 2026.3.4
python-louvain 0.16
ipython 9.15.0
