# pver_epibac
Scripts and important data for "Bacterial inoculation elicits changes to the coral epigenome"

## Methylation analysis
This repository contains the scripts (or references to scripts) and input files used for analyzing the epigenetic response of *Pocillopora verrucosa* to bacteria inoculation and heat stress in the paper "Bacterial inoculation elicits changes to the coral epigenome"

The methylation pipeline used to process samples to this point can be found at https://github.com/lyijin/working_with_dna_meth (which also fuller descriptions + theoretical considerations of the pipeline, and the scripts written to operate on Bismark's output)

### The files in the main repository (and how they were obtained):
- `all.filt.annot.merged.cov.gz`
 - Methylated positions were called, filtered, and merged across all samples (explained at https://github.com/lyijin/working_with_dna_meth)
- `compiled_median_meths.tsv.gz`
 - Median methylation levels were calculated on the gene level for each of the samples and compiled (explained at https://github.com/lyijin/spis_dna_meth/tree/master/bias_density_medians)
