# Transcriptomic aging during severe infection

Reproducibility code for the study **“Immune-cell redistribution explains transcriptomic aging signatures during severe infection.”**

The project reanalyzes two public RNA-seq datasets:

- **GSE152641**: severe COVID-19 and healthy controls
- **GSE196399**: severe non-COVID-19 pneumonia and healthy controls

The workflow estimates transcriptomic age with the blood-specific **RNAAgeCalc** model,
inferes leukocyte composition using **quanTIseq** through **immunedeconv**, fits disease-association
and composition-adjusted models, performs gene-level modeling with **limma**, and carries out
Hallmark pathway enrichment with **fgsea**.

## Repository structure

```text
R/
  00_setup.R                    install required packages
  01_rnaagecalc.R               raw counts -> RNAAgeCalc -> delta
  02_quantiseq_deconvolution.R  expression -> quanTIseq fractions
  03_primary_models.R           disease models before/after cell adjustment
  04_gene_level_limma.R         composition-adjusted age coefficients
  05_fgsea_hallmark.R           Hallmark GSEA + leading-edge genes
  99_sessionInfo.R              frozen R/package information
archive/                        recovered original project scripts, unedited
data/                           local inputs; large matrices are not committed
results/                        generated outputs; not committed by default
docs/METHODS_MAPPING.md         manuscript-to-code map
run_all.R                       sequential runner
```

## Reproducibility and provenance

The scripts in `archive/` are the recovered project scripts and are retained **unchanged** for
provenance. The numbered scripts in `R/` are cleaned, modular versions intended for public
release and reviewer inspection. They avoid machine-specific Desktop paths, use explicit input
files, export full regression coefficients/standard errors/P values, and export fgsea leading-edge
genes.

The corrected age- and sex-adjusted disease analysis yields positive disease coefficients in both
cohorts: approximately **+4.92 years** for severe COVID-19 and **+7.94 years** for severe pneumonia.

### Important deconvolution provenance note

The original analysis supplied DESeq2 variance-stabilized expression matrices to quanTIseq through the immunedeconv interface. This preprocessing choice is retained in the present repository to ensure that the deposited code accurately reflects the analysis underlying the reported results. The transcriptomic-age analysis was performed separately using raw gene-level counts as input to RNAAgeCalc, with the package performing its required internal normalization and transformation.

## Inputs

The raw RNA-seq data are public and are **not redistributed** here. Download the corresponding
count/expression files from GEO and prepare the filenames defined in `R/config.R`.

Curated metadata must contain at least:

- `sample`
- `class`
- `age`
- `sex`

Disease labels expected by `R/03_primary_models.R` are `Severe COVID-19` and `Severe pneumonia`,
with `Healthy` as the reference. Edit these labels if your metadata uses different strings.

## Running the workflow

From an R session opened at the repository root:

```r
source("R/00_setup.R")
source("run_all.R")
```

Or run each numbered script separately. This is preferable during manuscript revision because
it allows each stage to be inspected before downstream analyses are regenerated.

## Outputs

Key outputs include:

- sample-level RNAAge and delta
- quanTIseq neutrophil, monocyte and aggregated lymphocyte fractions
- complete regression coefficients, standard errors and P values
- composition-adjusted limma results for each cohort
- cross-cohort directional concordance table
- Hallmark NES, P value, FDR and leading-edge genes
- `sessionInfo.txt`

## Software

Core R packages: RNAAgeCalc, immunedeconv, limma, fgsea, msigdbr, DESeq2, dplyr, tibble.
Run `R/99_sessionInfo.R` after analysis to freeze the exact versions used in the final revision.

## Citation

Please cite the associated manuscript and the archived GitHub/Zenodo release. Update
Zenodo DOI before release.

## License

MIT.
