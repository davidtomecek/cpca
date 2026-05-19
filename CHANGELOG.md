# Changelog

## [2026-05-19]

### Added
- **MATLAB-based voxel-wise FD pipeline** using DPARSFA/SPM12 framework:
  - `dparfsa_voxel_based_fd.m`: Headless pipeline for voxel-specific head motion calculation
  - Generates per-voxel FD maps (`FDvox_4DVolume.nii`, `MeanFDvox.nii`)
  - Outputs directional motion components (X, Y, Z) as 4D NIfTI volumes
  - Requires MATLAB 2024a + SPM12 + DPABI toolboxes

## [2026-05-11]

### Added
- New Jupyter notebook pipeline for calculating, exporting, and visualizing Framewise Displacement (FD) metrics in the `code/calculate-fd/` directory.
  - `calculate_fd.ipynb`: Calculates FD from motion parameter files.
  - `export_metrics.ipynb`: Aggregates FD data and identifies motion outliers.
  - `plot_agg_fd.ipynb`: Visualizes the distribution of average FD across subjects.
