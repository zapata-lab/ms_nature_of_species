## Overview

This directory includes all the raw data as well as intermediate and final files needed to reproduce all our analyses. The only files not available here are the original sequence reads in fastq format, which are available via [SRA](https://www.ncbi.nlm.nih.gov/sra). Please visit the `analysis` directory for the code and scripts used to generate and use the files included here. If you use any of our data, **please cite our work**. If you would like to carry out a substantial re-analysis or re-interpretation, we would certainly appreciate you contact [us](fzapata@ucla.edu) first. Collecting all these data took us many years and substantial resources!!

Below we provide a description of each file. For clades-specific files, we describe files for CladeI but similar files are provided for all clades.   

* `morphometrics_mspub.csv`  Specimen-based raw morphological measurements and associated metadata.
  
* `speciesTrees/Escallonia_tipsForSppTree_2pertaxon.loci`  original loci file (from iPyrad) for species tree analysis using 2 individuals per species

* `speciesTrees/Escallonia_tipsForSppTree_2pertaxon.vcf`  original vcf file (from iPyrad) used in downstream species tree analyses using 2 individuals per species

* `speciesTrees/Escallonia_tipsForSppTree_2pertaxon_filtered.vcf`  filtered vcf file (after running VCTfools) fused in downstream species tree analyses using 2 individuals per species

* `speciesTrees/Escallonia_tipsForSppTree_2pertaxon_filtered_concat.nex`  nexus file of concatenated loci for species tree analysis using 2 individuals per species. This is the input file for PAUP*; it contains concatenated sequences, loci partitions, taxon partitions, and svdquartets commands 

* `speciesTrees/Escallonia_tipsForSppTree_2pertaxon_filtered_handmade.gphocs`  intermediate gphocs file to generate independent loci files from loci in the filtered dataset for analysis using 2 individuals per species

* `speciesTrees/Escallonia_tipsForSppTree_4pertaxon.loci` original loci file (from iPyrad) for species tree analysis using 4 individuals per species

* `speciesTrees/Escallonia_tipsForSppTree_4pertaxon.vcf` original vcf file (from iPyrad) used in downstream species tree analyses using 4 individuals per species

* `speciesTrees/Escallonia_tipsForSppTree_4pertaxon_filtered.vcf`  filtered vcf file (after running VCTfools) used in downstream species tree analyses using 4 individuals per species

* `speciesTrees/Escallonia_tipsForSppTree_4pertaxon_filtered_concat.nex` nexus file of concatenated loci for species tree analysis using 4 individuals per species. This is the input file for PAUP*; it contains concatenated sequences, loci partitions, taxon partitions, and svdquartets commands

* `speciesTrees/Escallonia_tipsForSppTree_4pertaxon_filtered_handmade.gphocs`  intermediate gphocs file to generate independent loci files from loci in the filtered dataset for analysis using 4 individuals per species

* `speciesTress/params-Escallonia_tipsForSppTree_2pertaxon.txt`  parameters file used in iPyard to generate vcf file for species tree analysis using 2 individuals per species

* `speciesTress/params-Escallonia_tipsForSppTree_4pertaxon.txt`  parameters file used in iPyard to generate vcf file for species tree analysis using 4 individuals per species

* `CladeI/CladeI_BFD_template.xml` template file created in Beauti for SNAPP and BFD* analysis. This file is modified by hand to specify each model to be tested in BFD*

* `CladeI/CladeI_concat.contree` IQ-Tree inferred phylogeny based on concatenated loci in the filtered dataset. Newick format

* `CladeI/CladeI_concat.mldist` pairwise sequence divergence for all samples; calculated during IQTree analysis

* `CladeI/CladeI_concord.cf.stat_loci` gene- and site- concordance factors recorded during IQ-Tree analysis

* `CladeI/CladeI_filtered.str` filtered SNP dataset used in STRUCTURE. The same file is used directly in RMaverick. By changing spaces for commas in this file, the input file for prabclust (GC model) can be generated.

* `CladeI/CladeI_filtered.vcf` filtered vcf file (after running VCTfools) used in downstream analyses

* `CladeI/CladeI_filtered_concat.phy` concatenated loci from the filtered dataset selected for downstream analysis (in all clades this was the smallest data matrix, with the least amount of missing data). Phylip format

* `CladeI/CladeI_filtered_handmade.gphocs` intermediate gphocs file used to create independent loci files from loci in the filtered dataset

* `CladeI/CladeI_filtered_informativeLoci_subBFD.nexus` nexus file in binary format suitable for analysis in BEAST

* `CladeI/CladeI_filtered_informativeLoci_subBFD.vcf` vcf file containing the most informative loci for the subset of samples used for BFD* analysis

* `CladeI/CladeI_imap.txt` mapping file used in BPP analyses linking samples to demes. Note: if STRUCTURE and RMaverick supported different models, there will be a version of this file, one for each analysis; these files will have either `strRun` or `RMavRun` in the filename, respectively

* `CladeI/CladeI_ipyradOut.loci` original loci file output (from iPyrad)

* `CladeI/CladeI_ipyradOut.vcf` original vcf file output (from iPyrad)

* `CladeI/CladeI_loci2keep.txt` the most informative loci from the filtered dataset (determined by gCF/sCF analyses), with appended identifier tags for use in BPP analyses. Note: if STRUCTURE and RMaverick supported different models, there will be a version of this file, one for each analysis; these files will have either `strRun` or `RMavRun` in the filename, respectively

* `CladeI/CladeI_popFile.txt` file linking individuelas to taxonomic species

* `CladeI/CladeI_r1_a.ctl` Control file for BPP runs

* `CladeI/CladeI_sppHull_rangeAdj.csv` min/max traits values extracted from the monograph. File works as input for analysis of hypercubes.

* `CladeI/params-CladeI_cladeSpecificAnalysis.txt` arameters file used in iPyard to generate vcf file for all analyses 
