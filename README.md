# Distinguishing empty and nuclei-containing droplets

## Overview

This repository contains analysis scripts and simulation code for the manuscript **EmptyDropsMultiome discriminates real cells from background in single-cell multiomics assays**
by [Megas _et al._ (2024)](link).

## Downloading files

Download the following data:

- 10k PBMC multiome data tobe used from the 10x website. (https://www.10xgenomics.com/welcome?closeUrl=%2Fdatasets&lastTouchOfferName=PBMC%20from%20a%20Healthy%20Donor%20-%20Granulocytes%20Removed%20Through%20Cell%20Sorting%20%2810k%29&lastTouchOfferType=Dataset&product=chromium&redirectUrl=%2Fdatasets%2Fpbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0)
- gonadal development multiome data from [the Garcia-Alonso _et al._ study] (https://www.nature.com/articles/s41586-022-04918-4)

Place them under `data/input`

## Simulations

Run:

- `bash simulations/bsub_simrun_multiomics.sh`, which will perform simulations based on the PBMC dataset above to assess cell detection methods. This creates several new small cell types.
- `bash simulations/bsub_simrun_multiomics_mono.sh`, which will perform simulations based on the PBMC dataset above to assess cell detection methods. This creates only one new small cell type by downsamplind and scrambling monocytes.

After the simulations are complete, create the figures for the simulations by running:

- `bash simulations/bsub_figure2_monocytes.sh`
- `bash simulations/bsub_figure2.sh`

## Real data using EmptyDropsMultiome

Run:

- `bash real/bsub_real.sh`, which will apply cell detection methods to the gonadal development datasets.
- `bash real/bsub_downstream.sh`, which will apply QC on the detected droplets.

To create the figures run:

- `bash real/bsub_figure1.sh`
- `bash real/bsub_figure3.sh`
- `bash real/bsub_figure4.sh`


## Real data using only EmptyDrops

To perform analysis using EmptyDrops only on RNA, run:
- `bash real/bsub_real_rna.sh`, which will apply EmptyDrops and then QC controls.

To create the figures run:
- `bash real/bsub_figure_rna.sh`

