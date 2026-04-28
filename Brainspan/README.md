# BrainSpan canonical analysis

Reproduces the canonical BrainSpan rerun used for the Somatic5 developmental-expression analysis.

## What it does

Given the full BrainSpan RNA-seq download:
- verifies row mapping using `row_num` from `rows_metadata.csv` and column 0 from `expression_matrix.csv`
- assigns each sample to one of 11 developmental periods
- filters to protein-coding genes using valid `entrez_id`
- computes peak developmental period for every background gene
- extracts the 78 Somatic5 target genes and computes peak periods
- computes prenatal vs postnatal summary metrics for the target genes
- runs the brain-gene-background permutation test
- writes the heatmap input matrix used for the 78-gene period-level Z-score plot

## Input

Either:
- a `.zip` file containing `expression_matrix.csv`, `rows_metadata.csv`, `columns_metadata.csv`
- or a directory containing those three files

## Usage

```bash
python brainspan_canonical_analysis.py \
  --input /path/to/genes_matrix_csv.zip \
  --output-dir brainspan_canonical_rerun \
  --permutations 10000 \
  --seed 42
```

## Output files

- `all_genes_peak_periods.csv`
- `somatic78_canonical_peak_periods.csv`
- `somatic78_canonical_detailed.csv`
- `somatic78_zscored_by_period.csv`
- `permutation_null_distribution.csv`
- `brainspan_canonical_permutation_histogram.png`
- `brainspan_canonical_summary.csv`
- `brainspan_canonical_console_summary.txt`
- `run_metadata.json`

## Notes

- Row mapping is done by exact `row_num` match. The script does **not** use positional indexing.
- The first column of `expression_matrix.csv` is treated as the row identifier and is used as the DataFrame index.
- Developmental periods 1–4 are classified as prenatal.
- The background set is defined as genes with a non-empty `gene_symbol` and valid non-zero `entrez_id`, with duplicate symbols removed by keeping the first occurrence.
