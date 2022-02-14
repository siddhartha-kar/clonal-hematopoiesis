# Causes and consequences of clonal hematopoiesis in the UK Biobank

This is the repository accompanying the manuscript titled, "[Genome-wide analyses of 200,453 individuals yields new insights into the causes and consequences of clonal hematopoiesis](https://www.medrxiv.org/content/10.1101/2022.01.06.22268846v1)".

## Directories

This repository has two directories:

- `variant_calling` which includes the `Nextflow` pipelines used to run Mutect2 and process the subsequent VCF files<sup>**</sup>. We used `Nextflow` version 20.07.1 to run the pipelines.
  
- `Mendelian_randomization` which includes the R scripts and data used for the Mendelian randomization analyses.
  - Use `install.packages("tidyverse")`, `install.packages("vroom")`, and `install.packages("remotes")`, followed by `remotes::install_github("MRCIEU/TwoSampleMR")` to install all the R packages required for the scripts to run.
  - Filenames for the R scripts are in the format ``[CH overall or subtype-specific]_as_exp_[outcome category]_out.R``, where `as_exp` indicates that CH is the exposure and `out` refers to outcomes. Outcome categories include miscellaneous (`misc`), epigenetic (`epigen`), and leukocyte telomere length (`tl`). The outcome datasets can be accessed either through the tab-delimited plain text files ending in `...ids.txt` or can be downloaded via the links provided within the relevant R scripts.

<sup>**</sup>N.B. To apply for access to the exome sequence data please visit the [UK Biobank](https://www.ukbiobank.ac.uk/).
