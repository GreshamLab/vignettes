---
title: "Windchime pipeline"
author: "Pieter Spealman"
date: "11/04/2021"
output: html_notebook
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

This 

```
 conda create -n umi_tools python=3.8
 conda activate umi_tools
 conda install -c bioconda umi_tools
```

```
python windchime.py -a -im 500 -i ctrl_file.tab -o run_windchime_star.sh
```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.

conda install -c bioconda umi_tools
