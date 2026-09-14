# somatic-spectrum

Reproducible Genomic Structural Equation Modelling (GenomicSEM) analysis of a shared genetic dimension across fibromyalgia, multisite chronic pain (MCP), broad migraine, myalgic encephalomyelitis/chronic fatigue syndrome (ME/CFS), and irritable bowel syndrome (IBS).

The primary Somatic5 model uses the FinnGen R12 broad-migraine endpoint `G6_MIGRAINE` (`Migraine` in the analysis manifest; N = 357,295, comprising 28,504 cases and 328,791 controls). The separate `MigAura` manifest row represents the FinnGen migraine-with-aura endpoint and is retained as an auxiliary input; it is not part of the Somatic5 core model.

Downstream gene, tissue, and cell-type analyses characterise the shared factor's biological annotations. Brain-tissue and prenatal and adult neuronal-cell enrichments are interpreted as convergent associations, not as evidence of prenatal causality or cerebellar specificity.

The repository documentation is aligned with the final HGG submission snapshot `Lee_Somatic5_v23_HGG_submission_ready.zip`. No analysis code or generated results were changed for this documentation update.

See [GenomicSEM/README.md](GenomicSEM/README.md) for the reproducible analysis command sequence and model definitions.
