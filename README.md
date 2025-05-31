# Metastatic progression of pheochromocytoma and paraganglioma occurs via parallel evolution
## Abstract
Pheochromocytoma (PCC) and paraganglioma (PGL) are neuroendocrine tumors derived from chromaffin cells of the adrenal medulla and ganglia of the autonomic nervous system. Approximately one-third are causatively associated with germline pathogenic variants. Metastatic disease develops in up to 25% of patients with PCC/PGL, for whom therapeutic options are limited and no targeted treatments exist. Tumor evolution has not been well-delineated in metastatic PCC/PGL. We performed whole exome sequencing of paired specimens from 27 patients with metastatic PCC/PGL to better understand cancer progression. Tumors demonstrated high rates of loss-of-function variants in chromatin remodeling and DNA damage repair genes, suggesting potential therapeutic targets. Low rates of shared somatic variants were observed between primary tumors and metastases, with evidence of independent monoclonal pathogenic variants in metastatic tumors. These findings suggest that PCC/PGL metastases develop through monoclonal seeding and parallel progression.
![graphical_abstract](https://github.com/user-attachments/assets/ec19ce79-bb9b-4eb9-a470-49e5f9018b4c)
## Repository Details
This repository contains all the source code necessary to reproduce the findings of "Metastatic progression of pheochromocytoma occurs via parallel evolution." The subdirectories of the `src` folder are organized in the same order as the Results section of the manuscript. A summary of the folder structures appears below.
### Folder Structure 
| Folder | Purpose |
| --- | --- |
| [01_survival_analysis](src/01_survival_analysis) | Source source code to reproduce the survival curves |
| [02_somatic_variant_analysis](src/02_somatic_variant_analysis) | Contains source code to create the somatic variant oncoplot |
| [03_cnv_analysis](src/03_cnv_analysis) | Contains source code to analyze copy number variants |
| [04_mutation_timing](src/04_mutation_timing) | Contains source code for molecular clock analysis |
| [05_supplementary_analysis](src/05_supplementary_analysis) | Contains source code for supplemental figures |
| [06_archive](src/06_archice) | Contains source code for analyses that were performed but not included in manuscript |
