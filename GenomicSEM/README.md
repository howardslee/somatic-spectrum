# GenomicSEM Analysis Pipeline

This directory contains the GenomicSEM workflow used to estimate shared genetic architecture across the somatic-spectrum traits, fit latent factor models, run factor GWAS, and perform reviewer-requested diagnostics.

Large input data, reference panels, and bulk result directories are intentionally not tracked in Git. They are expected to be present locally under `data/` and `results/`.

## Environment

Run commands from the `GenomicSEM/` directory. The original analyses were run with `Rscript` available from `/usr/local/bin/`.

```bash
cd GenomicSEM
export PATH="/usr/local/bin:$PATH"
```

If using the local Framework R installation that has `GenomicSEM` installed:

```bash
/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/bin/Rscript --version
```

## Inputs

The pipeline expects:

- `config/manifest.tsv`
- original GWAS summary statistics under `data/original/`
- harmonised/processed outputs under `data/processed/hg19/`
- LDSC reference files under `data/reference/`

The scripts write model outputs under `results/`.

## Core Preprocessing

Run these steps first:

```bash
Rscript scripts/01_convert_and_harmonise_sumstats.R
Rscript scripts/01b_munge_sumstats_genomicsem.R
Rscript scripts/02_genomicsem_ldsc_universe.R
Rscript scripts/02b_export_ldsc_matrices.R
```

## Primary Somatic5 Model

Somatic5 includes chronic pain, fibromyalgia, ME/CFS, IBS, and migraine.

```bash
Rscript scripts/03_factor_model.R \
  --core "Chronic_Pain,Fibromyalgia,MECFS,IBS,Migraine" \
  --outdir results/models/somatic5/

Rscript scripts/04_factor_GWAS.R \
  --core "Chronic_Pain,Fibromyalgia,MECFS,IBS,Migraine" \
  --factor-rds results/models/somatic5/factor_model.rds \
  --outdir results/models/somatic5 \
  --cores 40
```

## Sensitivity Models

### Somatic4

Somatic4 excludes chronic pain.

```bash
Rscript scripts/03_factor_model.R \
  --core "Fibromyalgia,MECFS,IBS,Migraine" \
  --outdir results/models/somatic4/

Rscript scripts/04_factor_GWAS.R \
  --core "Fibromyalgia,MECFS,IBS,Migraine" \
  --factor-rds results/models/somatic4/factor_model.rds \
  --outdir results/models/somatic4 \
  --cores 40
```

### Somatic4_CP

Somatic4_CP excludes fibromyalgia.

```bash
Rscript scripts/03_factor_model.R \
  --core "Chronic_Pain,MECFS,IBS,Migraine" \
  --outdir results/models/somatic4_CP/

Rscript scripts/04_factor_GWAS.R \
  --core "Chronic_Pain,MECFS,IBS,Migraine" \
  --factor-rds results/models/somatic4_CP/factor_model.rds \
  --outdir results/models/somatic4_CP \
  --cores 40
```

## Additional-Trait Models

These models add one external trait to the Somatic5 core and fit the factor model only.

```bash
Rscript scripts/03_factor_model.R \
  --core "Chronic_Pain,Fibromyalgia,MECFS,IBS,Migraine,Depression" \
  --outdir results/models/somatic5_plus_Depression/

Rscript scripts/03_factor_model.R \
  --core "Chronic_Pain,Fibromyalgia,MECFS,IBS,Migraine,ADHD" \
  --outdir results/models/somatic5_plus_ADHD/

Rscript scripts/03_factor_model.R \
  --core "Chronic_Pain,Fibromyalgia,MECFS,IBS,Migraine,PTSD" \
  --outdir results/models/somatic5_plus_PTSD/

Rscript scripts/03_factor_model.R \
  --core "Chronic_Pain,Fibromyalgia,MECFS,IBS,Migraine,SCZ" \
  --outdir results/models/somatic5_plus_SCZ/

Rscript scripts/03_factor_model.R \
  --core "Chronic_Pain,Fibromyalgia,MECFS,IBS,Migraine,Tinnitus" \
  --outdir results/models/somatic5_plus_Tinnitus/
```

## Reviewer Diagnostics

### Factor LDSC Intercept

The reviewer-requested LDSC intercept/ratio for the Somatic5 factor GWAS is produced from the saved factor GWAS output:

```bash
/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/bin/Rscript \
  scripts/05_factor_ldsc_intercept.R \
  --model-dir results/models/somatic5 \
  --ld data/reference/eur_w_ld_chr \
  --outdir results/models/somatic5/ldsc_factor \
  --trait-name Somatic5_factor
```

The summary output is:

- `results/models/somatic5/ldsc_factor/Somatic5_factor_ldsc_summary.tsv`
- `results/models/somatic5/ldsc_factor/Somatic5_factor_ldsc.log`

Current Somatic5 factor LDSC result:

```text
Mean chi2: 1.7559
Lambda GC: 1.5694
Intercept: 1.0014 (0.0100)
Ratio: 0.0018 (0.0133)
Observed-scale h2: 0.1057 (0.0036)
```

### Formal Residual Covariance Tests

To formally test whether pairwise residual covariances remain beyond the Somatic5 factor, run all-pairs residual-coupling tests:

```bash
/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/bin/Rscript \
  scripts/03_factor_model.R \
  --only-coupling \
  --ldsc-rds results/ldsc_universe/ldsc/ldsc_universe.rds \
  --outdir results/models/somatic5 \
  --test-all-pairs
```

This writes:

- `results/models/somatic5/residual_coupling_tests.tsv`
- `results/models/somatic5/model_summary.txt`

Current result: all 10 possible pairwise residual covariances were tested. Several were nominally significant, but none survived Benjamini-Hochberg FDR correction.

## Notes

- BrainSpan exploratory analysis files are kept local and are not tracked in this repository.
- Bulk GWAS inputs, reference panels, and generated full result directories are excluded from Git to avoid committing large data artifacts.
