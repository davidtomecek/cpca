# Complex Principal Component Analysis (cPCA) for fMRI Data

This repository contains a collection of R scripts designed to perform Complex Principal Component Analysis (cPCA) on fMRI data projected onto cortical surfaces.

## Overview

The primary goal of this codebase is to estimate the loading matrix $W$ using advanced stochastic power iteration methods, specifically focused on handling complex-valued signals derived from fMRI data. The implementation includes features for variance reduction and efficient handling of large datasets.

## Project Structure

The scripts are organized into versioned directories under `code/fmri_on_surface_scripts/`, representing different iterations of the algorithm and utility functions:

- `code/fmri_on_surface_scripts/20241122/`: Initial implementations for cPCA estimation on specific datasets.
- `code/fmri_on_surface_scripts/20250606/`: Updated implementations featuring the **PIMVR** (Principal Iterative Method for Variance Reduction) algorithm and improved utility functions.

## Core Algorithms

- **cPCA Estimation**: Utilizes stochastic power iteration or PIMVR approaches to estimate the loading matrix.
- **PIMVR**: An accelerated version of stochastic power iteration designed for improved variance reduction.

## Data Pipeline

The codebase implements a complete preprocessing and estimation pipeline:
1. **Data Loading**: Reading NIfTI data and FreeSurfer surface/annotation formats.
2. **Preprocessing**: 
    - Regressing out motion parameters (using `.mcdat` files).
    - Bandpass filtering and resampling.
    - **Motion Analysis**: Calculating Framewise Displacement (FD) metrics to assess head motion.
3. **Signal Transformation**: Applying the Hilbert transform to convert real-valued signals into complex-valued signals.
4. **Estimation**: Running the cPCA/PIMVR algorithms to find the principal components.

## Motion Analysis (FD Metrics)

The repository includes two complementary pipelines for head motion analysis in `code/calculate-fd/`:

### MATLAB/DPARSFA Pipeline (Voxel-wise FD)
- `dparfsa_voxel_based_fd.m`: Full preprocessing pipeline using SPM12 + DPABI
  - Slice timing correction, realignment, voxel-specific FD calculation
  - Outputs: 4D FD maps, directional motion components (X/Y/Z), mean FD/TD maps
  - **Requirements**: MATLAB 2024a, SPM12, DPABI V9.0+

### Jupyter Notebook Pipeline
- `calculate_fd.ipynb`: Calculates FD from motion parameter files (.txt/.mat)
- `export_metrics.ipynb`: Aggregates FD data across subjects, identifies motion outliers
- `plot_agg_fd.ipynb`: Visualizes average FD distribution across subjects

### Language & Core Libraries
- **Language**: R
- **Signal Processing**: `gsignal`, `pracma`, `complexlm`
- **Parallel Computing**: `doMC`, `doFuture`, `future`, `foreach`
- **Data Handling**: `bigmemory` (for large-scale data), `freesurferformats`

### Running the Scripts
Most scripts are executable via `Rscript`. For example:
```bash
Rscript code/fmri_on_surface_scripts/20241122/cpca4set_estimation.R
```
*Note: Some scripts may contain hardcoded paths that require adjustment to your local environment.*

## Implementation Notes
- **Memory Management**: The project leverages `bigmemory` and chunked iterations to manage the high RAM requirements of large fMRI datasets.
- **Parallelization**: Heavy use of `foreach` backends. For debugging purposes, it is recommended to set the number of workers to 1.
- **Complexity**: The algorithms specifically handle complex-valued matrices, necessitating careful management of real and imaginary components.
