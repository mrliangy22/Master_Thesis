# Master's Thesis Code

This repository contains code for my MSc thesis:

Applying Cointegration Analysis to Model Brain Connections in iEEG Data: Possibilities and Limitations

The project applies Vector Error Correction Models to high-dimensional intracranial EEG time series. The analysis includes diagnostic checks, preprocessing, cointegration rank estimation, sparse VECM estimation, ROI-level interaction analysis, cross-subject comparison, and visualization.

## Repository structure

- scripts/00_diagnostics: diagnostic checks, including stationarity tests, lag selection, Johansen eigenvalues, trace statistics, condition numbers, and sensitivity analysis
- scripts/01_preprocessing: preprocessing workflow for iEEG time series
- scripts/02_rank_estimation: cointegration rank estimation
- scripts/03_vecm_lasso: sparse VECM/LASSO estimation
- scripts/04_plotting: plotting and visualization scripts
- scripts/05_cross_subject_tests: cross-subject tests and final comparison analyses
- data: placeholder only; original data are not included
- docs: additional documents, notes, or PDFs

## Data availability

The original intracranial EEG data are not included due to data privacy and access restrictions.

## Main tools

R, time series analysis, VECM, cointegration, LASSO, and iEEG data analysis.