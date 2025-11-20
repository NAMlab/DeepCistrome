[![Python 3.8](https://img.shields.io/badge/Python-3.8-blue.svg?logo=python&logoColor=white)](https://www.python.org/downloads/release/python-380/)
[![TensorFlow 2.10](https://img.shields.io/badge/TensorFlow-2.10-orange.svg?logo=tensorflow&logoColor=white)](https://www.tensorflow.org/)
[![pyfaidx](https://img.shields.io/badge/pyfaidx-latest-green.svg)](https://github.com/mdshw5/pyfaidx)
[![matplotlib](https://img.shields.io/badge/matplotlib-latest-blue.svg?logo=python&logoColor=white)](https://matplotlib.org/)
[![pandas](https://img.shields.io/badge/pandas-latest-150458.svg?logo=pandas&logoColor=white)](https://pandas.pydata.org/)
[![bedtools](https://img.shields.io/badge/bedtools-latest-success.svg)](https://bedtools.readthedocs.io/)
[![DeepSHAP](https://img.shields.io/badge/DeepSHAP-SHAP-red.svg)](https://github.com/shap/shap)
[![Conda](https://img.shields.io/badge/Conda-latest-44A833.svg?logo=anaconda&logoColor=white)](https://docs.conda.io/)


## Modelling genetic variation effects in plant gene regulatory networks using transfer learning on genomic and transcription factor binding data

### Preparing the necessary data
1. Once you have forked this repository, download the DAP-seq peaks data from http://neomorph.salk.edu/PlantCistromeDB
into the data directory and unzip. The data you need is "Peaks in narrowPeak format (FRiP >=5%)" found using the link
above. The downloaded file is named **dap_download_may2016_peaks.zip**.
2. Run the "fetch_data_and_create_subfolders.sh" script: ``sh fetch_data_and_create_subfolders.sh``
3. Run the "prepare_tf_families.py" script: ``python prepare_tf_families.py``. This script will join peaks from TFs 
belonging to the same families into a single bed file. Therefore, this will produce 46 bed files.
4. Run the "prepare_overlap_matrix.py" script: ``python prepare_overlap_matrix.py``. This script creates an overlap 
matrix, which is a matrix that tells you which of the 46 TF families has experimental binding peaks for each 250 nt
window of the *Arabidopsis thaliana* genome.

### Training CNN models: Multi-label classifiers
6. Once the above steps are completed, you can now train models using the "train.py" script ``python train.py``. This
script will train models, di-nucleotide shuffled models (di-control) and si-nucleotide shuffled models(si-control).

### Model interpretation with SHAP
We have two main steps to interpret our models with SHAP. The first step computes importance scores using the
DeepSHAP/DeepExplain implementation and uses MoDisco to generate motifs. The second step extracts these motifs and 
saves them into an easier to use file for downstream R scripts that use for example MotifStack.
7. Run the "generate_predictive_motifs.py" script: "train.py" script ``python generate_predictive_motifs.py``. This will
generate motifs for each family using SHAP and MoDisco.
8. Run the "extract_motifs.py" script: ``python extract_motifs.py``

### Predictions-Mercartor4 enrichment analysis
To perform the mercator-predictions enrichment analysis, we need the following steps.

9. Using "Arabidopsis_thaliana.TAIR10.pep.all.fa" found within the data/proteome subdirectory, compute Mercator4 
functional annotations using Mercator4 ( https://www.plabipd.de/mercator_main.html ). Download the output and save it
in the data directory as *mercator4_output.txt*.
10. Run the "cluster_gene_by_prediction.py" script: ``python cluster_genes_by_prediction.py``. This script will create
the regulatory modules.
11. Run the "mercator4_cross_promoter_clusters_enrichment_analysis.py": This will perform enrichment of modules in 
Mercator functional groups ``python mercator4_cross_promoter_clusters_enrichment_analysis.py``.

### Effects of SNPs on binding profiles
We used SNPs from the AraGWAs catalog to perform this analysis. We have also uploaded the SNPs downloaded and used.

12. Run the "SNP_effects_on_predictions.py" script:  ``python SNP_effects_on_predictions.py``.

### Transfer learning and heat stress classification analysis in *Zea mays*
We used heat stress MOA-seq data from Liang et al., 2022. This analysis compares peaks which showed a positive fold 
change to those which showed negative fold changes. Precisely, it looks at regions on the genome were recorded MOA-seq
footprints increase or decrease and compares these two using the predicted binding profiles generated using our CNN
models trained on *Arabidopsis thaliana*.

13. Run the "stress_binding_classification.py" script: ``python stress_binding_classification.py``.

