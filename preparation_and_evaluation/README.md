
### Updated Folder Content Summary and Pipeline

The scripts in your folder can be categorized into five distinct phases, forming a cohesive bioinformatics analysis pipeline:

1.  **FASTA Sequence Pre-processing:** A robust set of Python scripts and a master shell script (`prepare_candidate_hubs.sz.1.sh`) are dedicated to preparing high-quality FASTA files for the deepCIS model. This includes tasks like parsing GFF annotations for promoter orientation, ensuring control/variant sequence pairing, filtering for sequence redundancy, and handling reverse complementation. This pipeline has been used to prepare the sequence targets candidates for Plant STARR seq.
2.  **DeepCIS Prediction:** Core Python scripts (`dCIS_predict.*.py`) handle the loading of the neural network model and converting cleaned FASTA sequences into one-hot encoded input for generating the raw TFBS prediction and probability matrices.
3.  **Post-processing & TFBS Event Detection:** Scripts focus on interpreting the raw model output, specifically for large input sequences (e.g., 1500bp) which are processed in chunks. This includes pinpointing TFBS events/gain-loss regions and generating summary statistics on these events.
4.  **Evaluation & Metrics:** A set of Python scripts focuses on calculating quantitative metrics of the deepCIS model's performance (e.g., MCC, F1-Score) and analyzing features like TFBS co-occurrence and IPM predictability using linear and polynomial regression.
5.  **Biological & Variant Analysis:** R and Python scripts perform specialized analyses on the results, such as bootstrap confidence intervals for STARR-seq data, network analysis of TFBS co-occurrence under different genetic variants, and comparative boxplots for biological context.

-----

### Input/Output and Script Order Table (English)
All required Input data are provided with the Supplementary Materials of Peleke et al., 2025 "Modelling genetic variation effects in plant gene regulatory networks using transfer learning on genomic and transcription factor binding data".
The evaluation scripts have beeen designed to start with the overall prediction results across the A. thaliana col-0 TAIR10 genome found in the data/all_actual and data/all_predictions file. This data can be found e.g. supplementary table 1 too.
The results on IPM occurences are proved in suppl. table 4.
The Plant STARR-seq Enrichment Data (CSV) is provided in suppl. table 11.

The following table lists the critical files in their logical order of execution, detailing their primary input and output. Execution of the scripts will generate vizualisations and tabular results which are provided in the study.

| Phase                        | File                                                | Input (Primary)                                              | Output / Result                                                                        |
| :--------------------------: | :-------------------------------------------------: | :----------------------------------------------------------: | :------------------------------------------------------------------------------------: |
| **I. Input Preparation**     | `prepare_candidate_hubs.sz.1.sh`                    | VCF, Ref FASTA, GFF3                                         | Orchestrates all subsequent FASTA creation/cleanup (e.g., \`Athal\_...\_ctrl-var.fa\`) |
|                              | `find_orientation_for_promoter.py`                  | FASTA, GFF3                                                  | FASTA with corrected strand information.                                               |
|                              | `filter_fasta_overlaps.py`                          | FASTA                                                        | FASTA filtered for coordinate redundancy.                                              |
|                              | `keep_pairs.py`                                     | FASTA                                                        | FASTA ensuring control/variant pair integrity.                                         |
| **II. DeepCIS Prediction**   | `split_sequence_chunks_to_250bp.1.py`               | Clean FASTA (\*\_ctrl-var.fa)                                | Chunked FASTA for prediction.                                                          |
|                              | `dCIS_predict.sz.py`                                | Chunked FASTA, Model (\*.h5)                                 | Raw TFBS Prediction CSVs (e.g., \`dCIS\_TFBS\_bound\_candidates.\*.csv\`)              |
| **III. Post-processing**     | `dCIS_interpret_overlap_chunk.15.py`                | Chunked Prediction CSV                                       | Pinpointed TFBS event regions and data for plots.                                      |
|                              | `dCIS_interpret_overlap_chunk_overview_multi.15.py` | Output of \`dCIS\_interpret\_overlap\_chunk.15.py\`          | Summary statistics of TFBS events.                                                     |
|                              | `calculate_coccurences.0.py`                        | all\_predictions.csv                                         | \`average\_cooccurrence\_bin.csv\` (Avg co-occurrence per TF)                          |
| **IV. Evaluation & Metrics** | `evaluate_TFBS_binding_bins.3.py`                   | all\_actual.csv, all\_predictions.csv                        | \`evaluation\_summary\_table.csv\` & Performance Plots.                                |
|                              | `test_TFBS_performance_metrics_lin.3.py`            | TFBS\_performance\_metrics\_complete\_sz25.csv               | Linear Regression Plots/Stats.                                                         |
|                              | `test_TFBS_performance_metrics_poly.4.py`           | TFBS\_performance\_metrics\_complete\_sz25.csv               | Polynomial Regression Plots/Stats.                                                     |
| **V. Biological Analysis**   | `starr_seq_araGwas_analyses_extended.17.py`         | STARR-seq Enrichment Data (CSV)                              | Bootstrap CI, T-test, and Regression results.                                          |
|                              | `dCIS_motif_multi_variant_comparison.sz.0.R`        | \`dCIS\_TFBS\_bound\_candidates.\*.csv\` (Multiple variants) | Consolidated TFBS Co-occurrence Network Graph.                                         |

-----

### Proposed Code Structure for Git Repository

To ensure clarity, reproducibility, and ease of use for other researchers, I recommend adopting a standard project structure. This proposal logically groups your scripts, data, and models.

``` 
dCIS_evaluation_code/
├── README.md               # CRITICAL: Project summary, installation/setup steps, and pipeline overview.
├── LICENSE                 # License file (e.g., MIT, GPL).
├── requirements.txt        # List of Python dependencies (pip install -r requirements.txt).
├── setup_packages.R        # R script to install R dependencies (install.packages(...)).
| 
├── src/                    # All executable source code.
|   ├── 01_preprocessing/   # Scripts for initial data/FASTA preparation (Phase I)
|   |   ├── prepare_candidate_hubs.sh # Renamed for simplicity
|   |   └── ... (all other FASTA-related python scripts)
|   |
|   ├── 02_core_model/      # DeepCIS prediction and core logic (Phase II)
|   |   ├── dCIS_predict.py
|   |   └── dCIS_predict_probs.py
|   |
|   ├── 03_postprocessing/  # Interpretation, Event Detection, Co-occurrence (Phase III)
|   |   ├── dCIS_interpret_overlap_chunk.py
|   |   ├── dCIS_interpret_overlap_chunk_overview.py
|   |   └── calculate_coccurences.py
|   |
|   └── 04_analysis_eval/   # Evaluation, Metrics, and Biological Analysis (Phase IV & V)
|       ├── python/
|       |   ├── evaluate_TFBS_binding_bins.py
|       |   ├── starr_seq_araGwas_analyses.py
|       |   └── test_TFBS_performance_metrics_lin.py (and poly.py)
|       └── R/
|           ├── dCIS_motif_multi_variant_comparison.R
|           └── Kruskal-Wallis_Wilcoxon_Boxplots.R
|
├── data/                   # Input data, often kept separate due to size.
|   ├── raw/                # Large, static input files (GFF3, VCF, Genome FASTA, etc.)
|   └── intermediate/       # Core intermediate files (e.g., all_predictions.csv, TFBS_performance_metrics_complete_sz25.csv)
|
├── models/                 # Neural network model files.
|   └── model_chrom_1_model.h5
|
└── output/                 # Generated plots, final tables, and summaries.

```

This structure clearly separates the different functions of the code and makes the entire pipeline understandable from a high level. You should replace the script names in the proposed structure with the full names you have, or rename the local files to the simplified names provided in the suggestion for cleaner imports and execution.
