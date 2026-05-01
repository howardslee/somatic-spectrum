#!/usr/bin/env bash

# Somatic4 Fibromyalgia,Migraine,MECFS,IBS
#Rscript scripts/03_factor_model.R --core "Fibromyalgia,Migraine,MECFS,IBS" --outdir results/models/somatic4/
#Rscript scripts/04_factor_GWAS.R --core "Fibromyalgia,Migraine,MECFS,IBS"  --factor-rds results/models/somatic4/factor_model.rds  --outdir results/models/somatic4 --cores 40

# Somatic4_CP Chronic_Pain,Migraine,MECFS,IBS (replaces Fibro with Chronic Pain)
#Rscript scripts/03_factor_model.R --core "Chronic_Pain,Migraine,MECFS,IBS" --outdir results/models/somatic4_CP/
#Rscript scripts/04_factor_GWAS.R --core "Chronic_Pain,Migraine,MECFS,IBS"  --factor-rds results/models/somatic4_CP/factor_model.rds  --outdir results/models/somatic4_CP --cores 40

# Somatic5_Both Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS (replaces Fibro with Chronic Pain)
Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS" --outdir results/models/somatic5_Both/
Rscript scripts/04_factor_GWAS.R --core "Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS"  --factor-rds results/models/somatic5_Both/factor_model.rds  --outdir results/models/somatic5_Both --cores 40

# Somatic5_Depression Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,Depression 
Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,Depression" --outdir results/models/Somatic5_plus_Dep/
Rscript scripts/04_factor_GWAS.R --core "Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,Depression"  --factor-rds results/models/Somatic5_plus_Dep/factor_model.rds  --outdir results/models/Somatic5_plus_Dep --cores 40

# Somatic5_PTSD Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,PTSD 
Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,PTSD" --outdir results/models/Somatic5_plus_PTSD/
Rscript scripts/04_factor_GWAS.R --core "Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,PTSD"  --factor-rds results/models/Somatic5_plus_PTSD/factor_model.rds  --outdir results/models/Somatic5_plus_PTSD --cores 40

# Somatic5_Tinnitus Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,Tinnitus 
Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,Tinnitus" --outdir results/models/Somatic5_plus_Tinnitus/
Rscript scripts/04_factor_GWAS.R --core "Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,Tinnitus"  --factor-rds results/models/Somatic5_plus_Tinnitus/factor_model.rds  --outdir results/models/Somatic5_plus_Tinnitus --cores 40

# Somatic5_ADHD Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,ADHD
Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,ADHD" --outdir results/models/Somatic5_plus_ADHD/
Rscript scripts/04_factor_GWAS.R --core "Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,ADHD"  --factor-rds results/models/Somatic5_plus_ADHD/factor_model.rds  --outdir results/models/Somatic5_plus_ADHD --cores 40

# Somatic5_SCZ Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,SCZ
Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,SCZ" --outdir results/models/Somatic5_plus_SCZ/
Rscript scripts/04_factor_GWAS.R --core "Fibromyalgia,Chronic_Pain,Migraine,MECFS,IBS,SCZ"  --factor-rds results/models/Somatic5_plus_SCZ/factor_model.rds  --outdir results/models/Somatic5_plus_SCZ --cores 40
