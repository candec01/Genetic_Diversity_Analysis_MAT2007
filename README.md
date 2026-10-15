# Genetic Diversity vs. Extinction Risk (MAT2007)

Tests whether threatened animal species have lower COI nucleotide diversity (π)
than non-threatened species, using a Mann-Whitney U test implemented manually,
both globally and per taxonomic class.

# Dataset

Source: Petit-Marty, N., Vázquez-Luis, M., & Hendriks, I. E. (2020). *Use of the
nucleotide diversity in COI mitochondrial gene as an early diagnostic of
conservation status of animal species*. Conservation Letters.
https://doi.org/10.1111/conl.12756
(Supplementary DataSet1, PDF table of 4,363 species.)

The PDF table was converted to CSV; the result is included in `Data_set/`.
**No manual cleaning is needed**: the Threatened / Not Threatened classification
and the exclusion of species without a valid category (6 species) are done
automatically by the script.

# Requirements

- Python 3.x (developed with 3.14)
- pandas
- matplotlib >= 3.9 (the script uses `tick_labels`)

    pip3 install pandas "matplotlib>=3.9"

# How to run

From the project root folder (keep the `Data_set/` folder next to the script):

    python3 genetic_diversity_analysis.py

The script prints the cleaning summary, the test results (U, z, p-value, effect
size) and a per-class table, and opens several figures one after another.
Close each figure window to continue to the next step.
