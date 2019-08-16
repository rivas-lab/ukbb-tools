This code is designed to replicate the pipeline used in the UK Biobank marker paper (Bycroft et al.), just applied to populations other than the self-identified White British. It runs an initial PCA of self-identified individuals, detects outliers in this PCA, and then runs a second PCA to project everyone for use as covariates in analyses downstream. It also runs KING to filter out only non-outlier related individuals.