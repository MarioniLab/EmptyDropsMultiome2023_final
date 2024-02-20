# Distinguishing empty and nuclei-containing droplets

## Overview

This repository contains analysis scripts and simulation code for the manuscript **EmptyDropsMultiome discriminates real cells from background in single-cell multiomics assays**
by [Megas _et al._ (2024)](link).

## Downloading files

Download the following data:

- 10k PBMC multiome data tobe used from the 10x website. (https://www.10xgenomics.com/welcome?closeUrl=%2Fdatasets&lastTouchOfferName=PBMC%20from%20a%20Healthy%20Donor%20-%20Granulocytes%20Removed%20Through%20Cell%20Sorting%20%2810k%29&lastTouchOfferType=Dataset&product=chromium&redirectUrl=%2Fdatasets%2Fpbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0)
- gonadal development multiome data from [the Garcia-Alonso _et al._ study] (https://www.nature.com/articles/s41586-022-04918-4)

## Simulations

Enter `simulations` and run:

- `simrun.R`, which will perform simulations based on the real datasets to assess cell detection methods.
- `plotsim.R`, to recreate the plots in the manuscript for the simulation results.

## Real data

Enter `real` and run:

- `realrun.R`, which will apply cell detection methods to the real datasets.
- `negcheck.R`, which will examine the p-value distribution reported for low-count barcodes.

Each subdirectory of `analysis` contains self-contained analysis files for each data set.
Run:

- `analysis.Rmd`, which will perform the analysis of each dataset.
- `plot_maker.R`, which will create the plots used in the manuscript.

