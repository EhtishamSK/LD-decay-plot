# LD Decay Plot

Quick R script to visualize **linkage disequilibrium (LD) decay** using TASSEL LD output.

## What it does

- Reads LD results from TASSEL
- Fits a nonlinear regression model to LD vs. distance
- Estimates recombination parameter (rho = 4Nr)
- Calculates **half-decay distance**
- Creates LD decay plot (scatter + fitted curve)

## Requirements

- R
- Packages: ggplot2 (and possibly others like dplyr – check script)

## How to use

1. Have your TASSEL LD output file ready (usually with columns like distance, r², etc.)

2. Open `script.R` and adjust:
   - Input file name/path
   - Column names if needed

3. Run the script

Output:
- LD decay plot (saved as image)
- Key values: half-decay distance, rho estimate

## Model

Nonlinear fit to model LD decay over physical/genetic distance.

## Notes

- Designed for TASSEL-generated LD tables
- Useful in population genetics, GWAS, haplotype block analysis
- Clean and filter your LD data first if necessary
