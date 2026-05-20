# Postprandial profiling of the duodenal microbiome

Custom code accompanying:

> Warren, F.J., Petropoulou, K., Harris, H.C., Barbas-Bernardos, C., Kasapi, M.,
> Garcia, A., Holmes, E., Domoney, C., Wist, J., Garcia-Perez, I. and Frost, G.
> "Postprandial profiling of the duodenal microbiome reveals the impact of food
> structure and association with luminal metabolite and gut hormone responses."

Source repository: <https://github.com/fredoowarren/Duodenal_microbiome_CRESTAR> &nbsp;·&nbsp;
Preprint: <https://www.biorxiv.org/content/10.64898/2026.05.06.723166v1>

---

## Contents

```
.
├── README.md                            this file
├── LICENSE                              CC0 1.0 Universal public-domain dedication
├── duodenal_microbiome_analysis.R       all statistical analyses + Figures 1–6
└── data/
    ├── physeq_hormones.rds              phyloseq object: MetaPhlAn4 species
    │                                    abundances + metadata + taxonomy
    ├── metabolites_hormones.xlsx        quantified NMR metabolites,
    │                                    bile acids, GIP, GLP-1
    ├── metabolites_plot.xlsx            long-format per-metabolite data
    │                                    used by the trajectory plots
    └── amino_acids.xlsx                 chiral LC–MS D/L amino acid
                                         concentrations (0 and 60 min)
```

The upstream shotgun metagenomics pipeline (KneadData → MetaPhlAn 4 →
HUMAnN 3) is not included here. Those steps are wrappers around standard,
cited bioinformatics tools and their exact parameters are documented in
the Methods section of the manuscript.

A small **demo dataset** (`data/demo/`) containing four randomly selected
participants is bundled with the repository so reviewers can run the analysis
end-to-end in a few minutes. The raw sequencing reads are deposited in NCBI
SRA under **PRJNA1425766**.

---

## 1. System requirements

### Operating systems

The R analysis (`duodenal_microbiome_analysis.R`) runs on any platform
supported by R 4.5.x: Linux, macOS 12+, Windows 10/11.

The R analysis has been tested on:

- Ubuntu 22.04 LTS, R 4.5.1
- macOS 14 (Sonoma), R 4.5.1
- Windows 11, R 4.5.1

### Hardware

No non-standard hardware is required. Tested on a 2020 laptop with a 4-core
CPU and 16 GB RAM. The Figure 5 correlation step is the heaviest in-memory
operation; 8 GB will be tight for the full dataset but is comfortable for
the bundled demo.

### Software dependencies

R ≥ 4.5.1 with the following packages (versions shown are the ones tested):

| Package        | Version  | Source       |
|----------------|----------|--------------|
| phyloseq       | 1.50.0   | Bioconductor |
| vegan          | 2.6-8    | CRAN         |
| lme4           | 1.1-35   | CRAN         |
| lmerTest       | 3.1-3    | CRAN         |
| compositions   | 2.0-8    | CRAN         |
| Hmisc          | 5.1-3    | CRAN         |
| igraph         | 2.0.3    | CRAN         |
| ggraph         | 2.2.1    | CRAN         |
| ggforce        | 0.4.2    | CRAN         |
| ggplot2        | 3.5.1    | CRAN         |
| patchwork      | 1.2.0    | CRAN         |
| ggrepel        | 0.9.5    | CRAN         |
| ggtext         | 0.1.2    | CRAN         |
| pheatmap       | 1.0.12   | CRAN         |
| tidyverse      | 2.0.0    | CRAN         |
| readxl         | 1.4.3    | CRAN         |
| RColorBrewer   | 1.1-3    | CRAN         |
| viridis        | 0.6.5    | CRAN         |
| scales         | 1.3.0    | CRAN         |

> Replace these with the exact versions returned by `sessionInfo()` on the
> machine that produced the published figures before submission.

---

## 2. Installation guide

Install R 4.5.x from <https://cran.r-project.org/>, then from an R session:

```r
# CRAN packages
install.packages(c(
  "tidyverse", "readxl", "vegan", "lme4", "lmerTest", "compositions",
  "Hmisc", "igraph", "ggraph", "ggforce", "ggplot2", "patchwork",
  "ggrepel", "ggtext", "pheatmap", "RColorBrewer", "viridis", "scales"
))

# Bioconductor (phyloseq)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
```

**Typical install time on a normal desktop computer with broadband internet:
30 – 45 minutes** (most of which is compiling C/C++ dependencies of `Hmisc`,
`lme4`, and `phyloseq`; on Windows and macOS, binary builds reduce this to
~10 minutes).

---

## 3. Demo

A small demo dataset (`data/demo/`) is bundled with the repository. It
contains four randomly chosen participants (one per meal arm) so that
reviewers can verify the analysis end-to-end without needing the complete
dataset.

### Running the demo

From the repository root, in R:

```r
# Point the script at the demo directory rather than the full data
Sys.setenv(DUODENAL_DATA = "data/demo")

# Run end-to-end
source("duodenal_microbiome_analysis.R")
```

(The script's `read*` calls read from `data/` by default; if you prefer not
to use the environment variable, simply copy the four files in
`data/demo/` over the equivalents in `data/`.)

### Expected output

A `results/` directory is created in the working directory containing:

```
results/
├── Figure1A_baseline_composition.pdf
├── Figure1B_baseline_pcoa.pdf
├── Figure2_alpha_diversity.pdf
├── Figure3A_composition_areaplot.pdf
├── Figure3B_pcoa_by_meal.pdf
├── Figure4A_diff_abundance_dotplot.pdf
├── Figure4B_diff_abundance_trajectories.pdf
├── Figure5_network.pdf
├── Figure6_DL_amino_acid_ratios.pdf
├── Sup_metabolite_*.pdf                (one per metabolite)
├── differential_abundance.csv
├── pairwise_meal_permanova.csv
├── time_stratified_permanova.csv
├── network_edges.csv
├── network_nodes.csv
└── DL_ratio_statistics.csv
```

PDFs of all six main figures should look like reduced-power versions of the
manuscript figures (fewer participants → noisier LOESS bands, broader
confidence ellipses, fewer significant network edges) but the structure
and the analysis flow are identical.

The console will also print:

- A PERMANOVA table for `Genotype × Food_Type + Participant` (Section 5)
- Type III ANOVA tables from each LME model (Sections 4 and 6)
- A summary of the network: number of nodes, edges, Louvain communities,
  and modularity (Section 7)

### Expected run time for the demo

**~2 – 3 minutes total** on a 2020-era 4-core laptop, broken down roughly as:

- Setup, data loading, alpha and beta diversity: ~15 s
- Differential abundance loop (top 50 taxa × LME): ~40 s
- CLR transform + Spearman + network construction: ~30 s
- All figure rendering (PDFs): ~30 s
- Supplementary metabolite panels: ~20 s

### Expected run time on the full dataset

**~10 – 20 minutes** on the same machine. The Hmisc::rcorr step in the
network section is the rate-limiting operation.

---

## 4. Instructions for use

### Reproducing the published figures

```bash
# Clone the repository
git clone https://github.com/fredoowarren/Duodenal_microbiome_CRESTAR.git
cd Duodenal_microbiome_CRESTAR

# Run the analysis (uses the data/ files included with the repo)
Rscript duodenal_microbiome_analysis.R
```

Outputs land in `results/` exactly as for the demo. All seeds are fixed
(`set.seed(123)`) so PERMANOVA p-values, Louvain community labels, and
network layouts are byte-for-byte reproducible across runs on the same
machine. (Layout coordinates are not guaranteed to match exactly across
operating systems because of differences in BLAS implementations, but
network topology and community membership will.)

### Re-generating the species-level abundance tables from raw FASTQ

Raw paired-end FASTQ files are available from NCBI SRA project
**PRJNA1425766**. The shotgun metagenomics pipeline used to derive the
species-level abundance tables (KneadData v0.10.0 → MetaPhlAn 4 → HUMAnN 3)
consists of standard, cited bioinformatics tools invoked with the parameters
documented in the Methods section of the manuscript. The pipeline is not
distributed as part of this repository because it contains no custom code,
only wrappers around the published tools.

### Applying the analysis to your own data

The R script expects four input files in `data/`:

1. **`physeq_hormones.rds`** — a phyloseq object with:
   - `otu_table`: species × sample relative abundances (percentage, not
     proportions), produced by MetaPhlAn4
   - `tax_table`: taxonomy with `Kingdom`, `Phylum`, `Class`, `Order`,
     `Family`, `Genus`, `Species` columns
   - `sample_data`: must contain columns `PPT` (participant), `Time`
     (minutes post-meal, numeric), and `SampleID` (one of `A`, `B`, `C`,
     `D` encoding RR-peas, rr-peas, RR-flour, rr-flour)

2. **`metabolites_hormones.xlsx`** — one row per sample, indexed by a
   column named `sample` whose values match the phyloseq sample names.
   Columns: `Valine`, `Alanine`, `Acetate`, `Methionine`, `Succinate`,
   `Lysine`, `Tyrosine`, `Phenylalanine`, `beta_Glucose`, `beta_Maltose`,
   `Stachyose`, `Trigonelline`, `Fumarate`, `Tryptophan`, `Formate`,
   `Propionate`, `Lactate`, `Threonine`, `BA1`–`BA7`, `GLP`, `GIP`.

3. **`metabolites_plot.xlsx`** — long-format `Participant`, `Time`,
   `Group`, then one column per metabolite.

4. **`amino_acids.xlsx`** — one row per (Participant, Timepoint,
   Randomization) with paired `D <aa>` and `L <aa>` concentration columns
   for each amino acid of interest.

To adapt to a different study design (e.g. different number of meal arms),
edit the `Meal_Full` mapping near the top of Section 2 of the R script.

---

## License

This code is released into the public domain under the
[**CC0 1.0 Universal**](https://creativecommons.org/publicdomain/zero/1.0/)
dedication (see `LICENSE`). You are free to copy, modify, and redistribute
without restriction.

---

## Citation

If you use this code, please cite:

> Warren, F.J. *et al.* (2026) Postprandial profiling of the duodenal
> microbiome reveals the impact of food structure and association with
> luminal metabolite and gut hormone responses. *Journal*, vol(issue),
> pages. DOI: `<paper DOI>`.

A `CITATION.cff` file is included for machine-readable citation metadata.

---

## Contact

Frederick J. Warren — fred.warren@quadram.ac.uk
Quadram Institute Bioscience, Norwich Research Park, Norwich NR4 7UQ, UK

Issues and pull requests are welcome at
<https://github.com/fredoowarren/Duodenal_microbiome_CRESTAR/issues>.
