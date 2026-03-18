# Gresham Lab Vignettes

Notebooks documenting common experimental and bioinformatics workflows used in the [Gresham Lab](https://greshamlab.bio.nyu.edu/) at NYU.

Rendered notebooks are available at **https://greshamlab.bio.nyu.edu/vignettes/**. Each `.nb.html` file embeds the source `.Rmd`, which can be downloaded directly from the website.

## Vignettes

| Directory | Topic |
|---|---|
| [`flow_cytometry/`](flow_cytometry/) | Flow cytometry analysis with CytoExploreR (Cytek Aurora) |
| [`growth_curves/`](growth_curves/) | Growth curve analysis from Tecan plate reader data using Growthcurver |
| [`reform/`](reform/) | Building custom genome and annotation files with [reform](https://github.com/gencorefacility/reform) |
| [`nanopore_tools/`](nanopore_tools/) | Nanopore basecalling (Guppy) and alignment (minimap2) on the HPC |
| [`nanopore_cdna_isoforms/`](nanopore_cdna_isoforms/) | cDNA isoform calling from Nanopore reads using ONT's Snakemake pipeline |
| [`windchime/`](windchime/) | Paired-end RNAseq with UMIs: trimming, STAR alignment, UMI-tools deduplication, and feature counts for DESeq2 |
| [`ggplot/`](ggplot/) | Introduction to `ggplot2` for data visualization |
| [`tidyverse/`](tidyverse/) | Introduction to tidy data and `dplyr` |

## Repository Structure

All vignettes are self-contained RMarkdown (`.Rmd`) notebooks. Supporting files (images, example data, config templates) are kept in the same subdirectory as their notebook. Reference files shared across vignettes are in [`data/`](data/).

## Contributing

Each vignette should be an `html_notebook` RMarkdown document. To render locally:

```r
rmarkdown::render("vignette_dir/notebook.Rmd")
```

Rendered `.nb.html` files are excluded from the repository via `.gitignore`.
