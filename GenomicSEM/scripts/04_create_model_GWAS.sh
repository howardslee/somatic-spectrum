#!/usr/bin/env bash

# Somatic4 Fibromyalgia,MigAura,MECFS,IBS
#Rscript scripts/03_factor_model.R --core "Fibromyalgia,MigAura,MECFS,IBS" --outdir results/models/somatic4/
#Rscript scripts/04_factor_GWAS.R --core "Fibromyalgia,MigAura,MECFS,IBS"  --factor-rds results/models/somatic4/factor_model.rds  --outdir results/models/somatic4 --cores 40

# Somatic4_CP Chronic_Pain,MigAura,MECFS,IBS (replaces Fibro with Chronic Pain)
#Rscript scripts/03_factor_model.R --core "Chronic_Pain,MigAura,MECFS,IBS" --outdir results/models/somatic4_CP/
#Rscript scripts/04_factor_GWAS.R --core "Chronic_Pain,MigAura,MECFS,IBS"  --factor-rds results/models/somatic4_CP/factor_model.rds  --outdir results/models/somatic4_CP --cores 40

# Somatic5_Both Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS (replaces Fibro with Chronic Pain)
Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS" --outdir results/models/somatic5_Both/
Rscript scripts/04_factor_GWAS.R --core "Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS"  --factor-rds results/models/somatic5_Both/factor_model.rds  --outdir results/models/somatic5_Both --cores 40

# Somatic5_Depression Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS,Depression 
Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS,Depression" --outdir results/models/Somatic5_plus_Dep/
Rscript scripts/04_factor_GWAS.R --core "Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS,Depression"  --factor-rds results/models/Somatic5_plus_Dep/factor_model.rds  --outdir results/models/Somatic5_plus_Dep --cores 40

# Somatic5_PTSD Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS,PTSD 
Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS,PTSD" --outdir results/models/Somatic5_plus_PTSD/
Rscript scripts/04_factor_GWAS.R --core "Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS,PTSD"  --factor-rds results/models/Somatic5_plus_PTSD/factor_model.rds  --outdir results/models/Somatic5_plus_PTSD --cores 40

# Somatic5_Tinnitus Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS,Tinnitus 
Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS,Tinnitus" --outdir results/models/Somatic5_plus_Tinnitus/
Rscript scripts/04_factor_GWAS.R --core "Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS,Tinnitus"  --factor-rds results/models/Somatic5_plus_Tinnitus/factor_model.rds  --outdir results/models/Somatic5_plus_Tinnitus --cores 40


Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS,ADHD" --outdir results/models/Somatic5_plus_Tinnitus/


Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS,SCZ" --outdir results/models/Somatic5_plus_SCZ/


Rscript scripts/03_factor_model.R --core "Fibromyalgia,Chronic_Pain,MigAura,MECFS,IBS,ADHD" --outdir results/models/Somatic5_plus_ADHD/








