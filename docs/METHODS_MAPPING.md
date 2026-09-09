# Methods-to-code mapping

| Manuscript analysis | Script |
|---|---|
| RNAAgeCalc blood-specific transcriptomic age | `R/01_rnaagecalc.R` |
| quanTIseq immune-cell deconvolution through immunedeconv | `R/02_quantiseq_deconvolution.R` |
| Base and leukocyte-adjusted disease models | `R/03_primary_models.R` |
| Composition-adjusted gene-level age associations with limma | `R/04_gene_level_limma.R` |
| Hallmark preranked GSEA with fgsea | `R/05_fgsea_hallmark.R` |
| R/package environment snapshot | `R/99_sessionInfo.R` |

## Definitions

Transcriptomic age acceleration is defined as:

`delta = RNAAge - chronological age`

The primary base model is:

`delta ~ class + age + sex`

The leukocyte-composition-adjusted model is:

`delta ~ class + age + sex + Neutrophil + Monocyte + Lymphocytes`

The gene-level model is:

`expression ~ age + sex + Neutrophil + Monocyte + Lymphocytes`

`Lymphocytes` is the arithmetic sum of quanTIseq B-cell, CD4 non-regulatory T-cell,
Treg and CD8 T-cell fractions. No post-hoc renormalization is performed in the curated script.
