# Derived inputs

This directory is for derived matrices and curated metadata used by the analysis.
It is excluded from git by default because the files may be large.

Expected inputs are configured in `R/config.R`:

- VST expression matrices for limma analyses
- normalized expression matrices used for quanTIseq deconvolution
- curated sample metadata with sample ID, disease class, chronological age and sex

Important provenance note: archived project scripts recovered for this repository passed VST matrices to `immunedeconv::deconvolute(method = "quantiseq")`. The curated workflow preserves this fact transparently through `DECONV_INPUT_SCALE` in `R/config.R`. If the final manuscript analysis is rerun using a different input scale (for example TPM), change the configuration and regenerate all downstream leukocyte-adjusted results before publication.
