[README.md](https://github.com/user-attachments/files/27264874/README.md)
# Juglans Population Genomics — Climate Adaptation & Assisted Gene Flow

R and shell scripts for population genomic analyses of *Juglans hindsii* (Northern California black walnut) and *J. californica* (Southern California black walnut), with a focus on identifying climate-adaptive genetic variation and evaluating its relevance to commercial walnut croplands under future climate scenarios.

---

## Repository contents

| File | Description |
|------|-------------|
| `JuglansAnalyses.R` | Master analysis script (both species combined): VCF filtering, PCA, ADMIXTURE, RDA, LFMM, Gradient Forest, genomic offset mapping, moving-window diversity, gene annotation |
| `JuglansAnalyses_JhinONLY.R` | Same workflow as above restricted to *J. hindsii* samples only |
| `ForwardAdaptedness_Jcal2Cropland.R` | Forward genomic offset — how well-adapted is *J. californica* to future climate at current walnut orchard locations |
| `ForwardAdaptedness_Jhin2Cropland.R` | Forward genomic offset — same analysis for *J. hindsii* |
| `ReverseAdaptedness_Jcal2Walnuts.R` | Reverse genomic offset — which wild *J. californica* individuals are most pre-adapted to future orchard conditions |
| `ReverseAdaptedness_Jhin2Walnuts.R` | Reverse genomic offset — same analysis for *J. hindsii* |

---

## Analysis overview

### 1 · VCF filtering & LD pruning
Raw SNP data is filtered for mean depth, individual genotype depth, MAF (≥ 0.01), and missingness (≤ 10%) using **vcftools**. **PLINK** is then used to prune for linkage disequilibrium (window 50 SNPs, step 10, r² < 0.1) and output `.bed/.bim/.fam` and `.raw` genotype files. Multiple filtering rounds remove outlier individuals iteratively.

### 2 · PCA
Principal components analysis is performed from PLINK eigenvectors in R (**ggplot2**), with points colored by latitude or longitude to visualize geographic structure.

### 3 · ADMIXTURE
Ancestry proportions are estimated for K = 1–10 using **ADMIXTURE** with cross-validation; optimal K is selected from CV error.

### 4 · Genotype–environment association (GEA)

**RDA + rdadapt** (`vegan`, `robust`, `qvalue`)  
Redundancy analysis with climate variables as predictors and population structure (PC1, latitude, longitude, lat × lon interaction) partialled out as a Condition. Outlier SNPs are identified via a Mahalanobis-distance test across RDA axes (rdadapt). The top 10% of outliers by p-value are retained for downstream analyses.

**LFMM** (`lfmm`, `qvalue`)  
Latent factor mixed-model ridge regression using climate PCs as predictors; K = 1 latent factor controls for population structure. Outliers from LFMM and RDA are merged into a combined adaptive SNP set.

### 5 · Gradient Forest (GF)
`gradientForest` models cumulative allele turnover along climate gradients using the adaptive outlier SNP set. The trained model is used to transform climate space so that distances reflect genotypic change.

### 6 · Genomic offset
For each future climate scenario — **CNRM-CM5** and **HadGEM2-ES** under **RCP 4.5** and **RCP 8.5**, at **2040–2069** and **2070–2099** — genomic offset is computed as the Euclidean distance in GF-transformed climate space between current and projected conditions. Results are written as GeoTIFFs.

### 7 · Forward & reverse adaptedness (cropland scripts)
Using the GF model trained on wild populations:

- **Forward offset** (recipient perspective): for each walnut cropland pixel (from the USDA Cropland Data Layer), calculate how different future climate will be from what current wild populations are adapted to. Low offset → cropland is likely to remain in climatic space that wild trees are adapted to.
- **Reverse offset** (donor perspective): for each wild individual, calculate how different its current climate is from future cropland conditions. Low offset → that wild individual carries alleles pre-adapted to future orchard environments and may be a useful seed/pollen donor.

Both directions are computed across all eight future scenarios; results are written as GeoTIFF rasters.

### 8 · Moving-window genetic diversity
`wingen` computes observed heterozygosity (Ho) and nucleotide diversity (π) across a spatial raster using a moving window. `krig_gd` then kriges these estimates to produce a smooth diversity surface. Analyses are run range-wide and within individual preserves.

### 9 · Isolation by distance / Mantel tests
Pairwise genetic distance matrices are correlated with geographic distance matrices using Mantel tests to assess IBD patterns.

### 10 · Gene annotation
Adaptive outlier SNPs are intersected with gene models using **bedtools closest**. Gene models are first lifted from the *J. regia* reference annotation onto the *J. hindsii* genome assembly using **Liftoff**.

### 11 · Environmental coverage
Density plots compare the distribution of climate variables across the full species range versus sampled locations, confirming that the sample set captures the species' climatic breadth.

---

## Software & R packages

**Command-line tools:** PLINK, vcftools, bcftools, ADMIXTURE, bedtools, Liftoff, conda

**R packages:** `terra`, `gradientForest`, `extendedForest`, `vegan`, `lfmm`, `qvalue`, `wingen`, `algatr`, `vcfR`, `adegenet`, `parallel`, `doParallel`, `foreach`, `fields`, `gstat`, `sp`, `raster`, `ggplot2`, `dplyr`, `tidyr`, `readr`, `readxl`, `psych`, `robust`, `Rcpp`

**Climate data:** Rose et al. current (1981–2010) and future (CNRM-CM5 / HadGEM2-ES, RCP 4.5 & 8.5) climate variables

**Crop mask:** USDA Cropland Data Layer (CDL) — walnut pixels

---

## Notes

- Scripts were developed for a SLURM/SGE HPC environment (`qrsh` header lines). File paths beginning with `/u/home/r/rcbuck/project-vlsork/` will need to be updated for your system.
- The `JuglansAnalyses.R` script is a combined analysis log and may include exploratory or deprecated code blocks (commented out). Species-specific downstream analyses live in the dedicated scripts.
- Chromosome accession IDs (CM084183.1–CM084198.1) in PLINK `.bim` files are renamed to integers 1–16 via `sed` before running ADMIXTURE.
