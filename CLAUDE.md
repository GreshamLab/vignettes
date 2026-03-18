# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is the [Gresham Lab](https://greshamlab.bio.nyu.edu/) vignettes repository — a collection of self-contained RMarkdown notebooks documenting common wet-lab and bioinformatics workflows. Each vignette knits to an HTML document published at https://greshamlab.bio.nyu.edu/vignettes/. These are tutorials, not an R package.

## Rendering Notebooks

Knit an individual notebook in R:
```r
rmarkdown::render("path/to/notebook.Rmd")
# or for html_notebook output:
rmarkdown::render("path/to/notebook.Rmd", output_format = "html_notebook")
```

Most notebooks have `eval = FALSE` in the setup chunk so they render without executing code (since they depend on HPC environments or large data files not in the repo).

## Repository Structure

Each vignette is either a root-level `.Rmd` or a subdirectory containing the `.Rmd` plus supporting files (images, CSVs, FASTAs):

| File/Directory | Topic |
|---|---|
| `vignette_SimpleFlow.Rmd` | Flow cytometry analysis with CytoExploreR (Cytek Aurora) |
| `growth_curves_gresham_lab.Rmd` | Growth curve analysis from Tecan plate reader using Growthcurver |
| `ggplot_notebook.Rmd` | ggplot2 tutorial |
| `tidyverse_notebook.Rmd` | tidyverse tutorial |
| `Reform/reform_genomes.Rmd` | Custom genome + GFF generation with the `reform` tool |
| `Common_Nanopore_Tools/nanopore.Rmd` | Nanopore basecalling (Guppy) and alignment (sbatch scripts) |
| `cDNA_Isosform_identifier/Nanopore_cDNA_isoform.Rmd` | Nanopore cDNA isoform calling via Snakemake pipeline |
| `windchime/Windchime_pipeline.rmd` | Paired-end RNAseq with UMIs: trim → align (STAR) → dedup (UMI-tools) → counts (bedtools) |

## Biological Context

All workflows are built around *S. cerevisiae* (and sometimes *C. albicans*) experiments in the Gresham Lab. Common themes:
- GAP1 locus, mCitrine fluorescent reporters, copy number variation
- Growth in YPD and glutamine-limited media
- HPC (NYU Greene) via SLURM `sbatch` for Nanopore/RNAseq pipelines
- Illumina SNV analysis is handled separately in the `vivaldi` R package (`~/projects/vivaldi/`)

## Adding a New Vignette

New vignettes follow the existing pattern: one `.Rmd` (using `html_notebook` output), with a setup chunk setting `eval = FALSE` if the code requires data/tools not in the repo. Supporting data files go in the same subdirectory as the `.Rmd`.
