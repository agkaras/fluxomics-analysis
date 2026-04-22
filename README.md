# Metabolomics Flux Analysis — ¹³C₆-Glucose Tracing

R analysis pipeline for isotope tracing metabolomics data from murine thoracic aorta.

**Experimental design:** 2 h ¹³C₆-glucose labelling · 24 h total incubation · aorta from young and old mice, CTR / IL-1β

---

## Data source

Representative results from this experimental dataset are published in:

> Karaś A, Buczek E, Pyka J, Annunciato I, Kutryb-Zając B, Jędrzejewska A, Stawarska K, Bar A, Kurpińska A, Nieminen AI, Szabo C, Kaczara P, Chłopicki S.
> **Age-dependent reprogramming of vascular metabolism compromises endothelial resilience to inflammation-induced endothelial dysfunction.**
> *Free Radical Biology and Medicine.* 2026;245:433–446.
> https://doi.org/10.1016/j.freeradbiomed.2025.12.046

This script was developed independently for visualisation and statistical analysis of the underlying flux data. The underlying raw data are not included in this repository.

---

## What this script does

- Loads and cleans raw metabolomics data from the FIMM platform (Excel, TIC-normalised)
- Reshapes data from wide to long format for statistical analysis and plotting
- Generates boxplots with individual data points overlaid for each metabolite
- Runs one-way ANOVA (4 groups) with Šidák post-hoc for biologically defined pairwise contrasts
- Saves individual per-metabolite PNGs, a combined multi-panel figure, a significance heatmap, and a clustered z-score heatmap of the top 50 most significantly changed metabolites

---

## Input data format

The script expects an `.xlsx` file with:

| Requirement | Detail |
|---|---|
| Sheet name | `Results normalized TIC` |
| Header rows to skip | 14 (data starts at row 15) |
| Required columns | `Pathway`, `Compounds`, `FLUX` + one column per sample |
| Sample column names | Must match the vectors defined in Section 3 |

Run `colnames(raw)` after loading to verify column names match your file.

---

## Quick start

```r
# 1. Set your file path at the top of the script
DATA_FILE <- "your_file.xlsx"

# 2. Update sample group vectors (Section 3) to match your column names
#    The example names currently in the script are defined basef on the published dataset 

# 3. Update metabolites_to_plot (Section 8) with compounds to show in the panel:
unique(data_clean$Compound)   # run this to see all available names

# 4. Source the script
source("fluxomics_analysis.R")
```

---

## Dependencies

```r
install.packages(c(
  "readxl",      # Excel import
  "tidyverse",   # data wrangling + ggplot2
  "ggbeeswarm",  # quasirandom point jitter
  "rstatix",     # statistical tests
  "ggpubr",      # significance brackets on plots
  "patchwork",   # multi-panel figure composition
  "emmeans",     # post-hoc contrasts (one-way ANOVA)
  "pheatmap"     # clustered heatmap
))
```

---

## Output

| File | Description |
|---|---|
| `metabolomics_flux_panel.png` | Multi-panel boxplot figure (300 dpi) |
| `plots_individual/*.png` | One boxplot PNG per metabolite, all compounds (300 dpi) |
| `heatmap_anova.png` | Significance heatmap: −log₁₀(FDR) for ANOVA + 4 contrasts |
| `heatmap_top50.png` | Clustered z-score heatmap, top 50 metabolites by ANOVA FDR |
| `anova_results.csv` | Full results table: F, p, FDR (omnibus) + FDR for each contrast |

---

## Statistical approach

One-way ANOVA (`Intensity ~ Group`, 4 groups) per metabolite. Post-hoc: three *a priori* biological contrasts with Šidák correction within each metabolite:

1. young CTR vs young IL-1β
2. old CTR vs old IL-1β
3. young CTR vs old CTR

FDR (Benjamini-Hochberg) applied across metabolites separately for each test.


---

## Plot style

- White background, no gridlines, single panel per metabolite
- X-axis: `young CTR` / `young IL-1β` / `old CTR` / `old IL-1β`
- Boxplot: median, IQR, whiskers with end-caps; semi-transparent fill
- Individual data points overlaid (quasirandom jitter), coloured by group
- Significance brackets from ANOVA Šidák post-hoc (* p < 0.05, ** p < 0.01, *** p < 0.001; Šidák post-hoc for 3 contrasts)

