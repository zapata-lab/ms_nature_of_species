## Analysis code

The code uses commands specific to the UCLA-HPC. It also uses helper scripts in the `bin` directory. Modify the code according to your HPC and read documentation of helper scripts.

For our genomic analysis, we conducted the following analyses:

1. Data assembly
2. Data filtering
3. Sensitivity Tests
4. Model-based species delimitation: GC Model
5. Model-based species delimitation: AC Model
6. Model-based species delimitation: RI Model
7. Model Selection

In all files, we use data for CladeI as a example. To reproduce all results, change the paths and file names for other clades.

For most analysis, we provide all the input files except when the input takes a large number of files which we cannot share efficiently here (thousands of files). However, we provide the code to generate such files if you want to regenerate all analyses.


Each file includes text and commands describing all our analysis.

* `1_data_assembly.sh` bash script to run iPyrad on the UCLA-HPC

* `2_data_filtering` commands and text to filter data and generate matrices with different amounts of missing data

* `3_sensivity_tests` commands and scripts (`R` and `bash`) to run PCA, infer phylogenies under maximum likelihood, and ADMIXTURE.

* `4_model_based_sp_delimitation_GC_model.R` `R` script to model species under GC model (using `prabclus`)

* `5_model_based_sp_delimitation_AC_model` commands and `R` script to model species under AC model (using `mPTP`)

* `6_model_based_sp_delimitation_RI_model` commands and scripts to  model species under RI model (using `BPP`)

* `7_model_selection` Text file describing how the analyses for model selection were set up for `Beast2`