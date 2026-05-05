# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview
This repository contains R scripts for performing Complex Principal Component Analysis (cPCA) on fMRI data projected onto cortical surfaces. The codebase is divided into two main versioned directories under `code/fmri_on_surface_scripts/`, representing different iterations of the implementation.

## Project Structure
- `code/fmri_on_surface_scripts/20241122/`: Contains scripts for cPCA estimation for specific datasets (e.g., `cpca4set_estimation.R`).
- `code/fmri_on_surface_scripts/20250606/`: Contains updated implementations, including more advanced PIMVR (Principal Iterative Method for Variance Reduction) algorithms and utility functions (`utils_pimvr.R`).

## Development Environment
- **Language**: R
- **Key Libraries**:
  - `gsignal`, `pracma`, `complexlm`: Signal processing and complex number manipulation.
  - `doMC`, `doFuture`, `future`: Parallel and distributed computing.
  - `freesurferformats`: Reading FreeSurfer surface and annotation formats.
  - `bigmemory`: Handling large datasets efficiently.
- **Dependencies**: Ensure R packages listed in the scripts (e.g., `RhpcBLASctl`, `freesurferformats`, `gsignal`) are installed.

## Common Tasks

### Running Scripts
Most scripts are designed to be run via `Rscript`.
- To run a specific estimation script:
  ```bash
  Rscript code/fmri_on_surface_scripts/20241122/cpca4set_estimation.R
  ```
- Note that many scripts contain hardcoded paths (e.g., `/home/anyzjiri/...`). These may need to be updated to match your local environment or passed as arguments if modified.

### Data Preprocessing
The scripts implement a pipeline including:
1. Loading NIfTI data.
2. Regressing out motion parameters (from `.mcdat` files).
3. Bandpass filtering and resampling.
4. Applying Hilbert transform to create complex-valued signals.

### Core Algorithms
- **cPCA Estimation**: Uses stochastic power iteration or PIMVR approaches to estimate the loading matrix $W$.
- **PIMVR**: An accelerated stochastic power iteration method for variance reduction.

## Architecture Notes
- **Parallelization**: The codebase heavily relies on `foreach` with `doMC` or `doFuture` backends. When debugging, you might want to set workers to 1 to trace execution more easily.
- **Memory Management**: Large datasets are processed using `bigmemory` or chunked iterations to manage RAM usage.
- **Complexity**: The algorithms deal with complex-valued matrices, requiring careful handling of real and imaginary parts.
