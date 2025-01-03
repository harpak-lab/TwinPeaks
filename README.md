# No Evidence for Sex-Differential Transcriptomes Driving Genome-Wide Sex-Differential Natural Selection
Matthew J. Ming, Changde Cheng, Mark Kirkpatrick, Arbel Harpak

Here are the code and scripts used to produce the analysis for "No Evidence for Sex-Differential Transcriptomes Driving Genome-Wide Sex-Differential Natural Selection"

## Overview:
1. All code can be found in the main directory. There are no subdirectories
2. The following software were used:
   * [R](https://www.r-project.org/)
   * [RStudio](https://posit.co/download/rstudio-desktop/)
   * [bcftools](https://samtools.github.io/bcftools/bcftools.html)
3. The following R packages were used:
   * library(ggplot2)
   * library(ggpubr)
   * library(vroom)
   * library(ggrepel)
   * library(grid)
   * library(pBrackets)

## Script documentation
Main code is contained in R scripts.
The Slurm script listed as script 2 is a wrapper which run R code over multiple files (e.g., multiple ancestries or tissue types). The slurm file is listed with a number (2) and the R scripts run by this slurm file is listed as number and letter (2a).
Scripts 3 and 4 each have 2 versions, one for all non-gonad tissues and one for gonads. The gonad version is listed as 3a and 4a respectively
All final figures were produced using ```TwinPeaksFigs_Rev_FINAL.R```.

### 1.ChengRepTwinPeaks_Anc.R
Get original Twin Peaks curve from 1000 Genomes data (Fig 2)

### 2. ChengRepRand_GTExONLY.slurm
Get Twin Peaks curves from sex label randomization in GTEx over multiple iterations (Fig 2)

### 2a. ChengRepTwinPeaks_Rand_GTEx.R
R script called by ```ChengRepRand_GTExONLY.slurm``` for getting Twin Peaks sex label randomization (Fig 2)

### 3. NewAnalysis_eQTLRegression_MH.R
Get A-values for updated linear regression method analysis in all tissues (except gonads) (Fig S1)

### 3a. NewAnalysis_eQTLRegression_MH_Gonads.R
Same as ```NewAnalysis_eQTLRegression_MH.R```, but modified slightly for the gonads specifically (Fig S1)

### 4. NewAnalysis_eQTLRegression_MH_Revision.R
Get A-values using an improved weighted Standard Major Axis regression method (Fig 3)

### 4a. NewAnalysis_eQTLRegression_MH_Gonads_Revision.R
Same as ```NewAnalysis_eQTLRegression_MH_Revision.R```, but modified slightly for the gonads specifically (Fig 3)

### 5. TwinPeaksFigs_Rev_FINAL.R
Final script for composing both main text and supplemental figures (Figs 1-3, S1). Depends on output from all previous scripts
