# Direct Gaussian Process Predictive Regressions with Mixed Frequency Data

## Description

> These codes come without technical support of any kind. The code is free to use, provided that the paper is cited properly.

This repository contains code and replication files for the paper "Direct Gaussian Process Predictive Regressions with Mixed Frequency Data" by Niko Hauzenberger, Massimiliano Marcellino, Michael Pfarrhofer and Anna Stelzer (working paper available at https://arxiv.org/abs/2402.10574).

## Contents

| File | Description |
|---|---|
| `example.R` | Worked example: nowcast and one-quarter-ahead forecast of US real GDP growth and GDP deflator inflation, with data pulled from ALFRED |
| `gpmf_func.R` | `npmf()`: direct mixed-frequency predictive regressions with a Gaussian process, BART or linear (horseshoe) conditional mean, several MIDAS compression schemes and optional stochastic volatility |
| `mfvar_func.R` | `mfvar()`: mixed-frequency VAR for monthly and quarterly variables, used as a benchmark |
| `aux_func.R` | Helper functions (horseshoe prior, MIDAS weighting polynomials, kernels, quantile-weighted CRPS) |
| `facpc_func.R` | Principal components and EM factor estimation, and selection of the number of factors |
| `replication/` | Script and intermediate data that reproduce the figures and tables of the main text; see `replication/README.md` |

## Example

Run `example.R` with the root of this repository as working directory:

```
Rscript example.R
```

The script downloads the monthly predictors and the quarterly targets as available on the vintage date (today by default) from ALFRED, so it needs an internet connection. It then estimates the model for the nowcast and the one-quarter-ahead forecast, adds a backcast if the previous quarter has not been released yet, prints a summary of the predictive distributions and plots them next to the latest released values.

The settings at the top of the script select the vintage date, the MIDAS compression (`run.type`), the conditional mean (`run.model`), stochastic volatility (`run.sv`) and the number of MCMC draws. The defaults use a shorter chain than the paper.

The data differ from those used in the paper (FRED-MD/QD) in two respects, which are also noted in the script: the S&P 500 is not included, and the sample starts in 1967.

Required R packages: `alfred`, `zoo`, `lubridate`, `ggplot2`, `forecast`, `stochvol`, `dbarts`, `MASS`, `spam`, `fields`, `abind` and `pracma`.

`mfvar()` is not used in `example.R` but features in the paper. It needs the packages `Matrix` and `stochvol`. The function has no forecast horizon argument; the periods to nowcast or forecast are set by the rows of the input matrix. The requirements on the input are listed at the top of `mfvar_func.R`.

## Replication

`replication/out_script.R` reproduces the figures and tables of the main text from two CSV files with the evaluated forecast losses and model confidence set results. The underlying collected forecast evaluation files are large and available upon request. Details are in `replication/README.md`.

## Contact

* Michael Pfarrhofer (michael.pfarrhofer@wu.ac.at, mpfarrho@gmail.com)
