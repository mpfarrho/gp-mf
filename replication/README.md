# Replication of the main-text figures

This folder reproduces the eight figure files used in the main text of the GP-MF paper with a single script, `out_script.R`.

## What is reproduced

All files are written to `figures/`, in the same subfolders and with the same names as in the paper.

| File | Paper item |
|---|---|
| `expl_lvl/mean/expl_mean_GDPC1_RMSE_Full_outl.pdf` | Figure 1, panel (a) `GDPC1`, left |
| `expl_lvl/mean/expl_mean_GDPC1_CRPS_Full_outl.pdf` | Figure 1, panel (a) `GDPC1`, right |
| `expl_lvl/mean/expl_mean_GDPCTPI_RMSE_Full_outl.pdf` | Figure 1, panel (b) `GDPCTPI`, left |
| `expl_lvl/mean/expl_mean_GDPCTPI_CRPS_Full_outl.pdf` | Figure 1, panel (b) `GDPCTPI`, right |
| `mcs/mcs-shares_25_full_CRPS_all.pdf` | Figure 2, panel (a), full sample |
| `mcs/mcs-shares_25_full-outl_CRPS_all.pdf` | Figure 2, panel (b), full sample without outliers |
| `tabs/bestmods_GDPC1_CRPS_full.pdf` | Table 1 (`GDPC1`) |
| `tabs/bestmods_GDPCTPI_CRPS_full.pdf` | Table 2 (`GDPCTPI`) |

## Folder contents

```
replication/
  out_script.R          # the script
  raw/                  # source files for raw mode
    eval_out.rda        # collected forecast evaluation (about 660 MB; not in the repository)
    eval_out.md         # note on eval_out.rda
    mcs_out.rda         # collected model confidence set results (not in the repository)
    mcs_out.md          # note on mcs_out.rda
  data/                 # intermediate CSVs
    losses_avg.csv      # subsample-averaged losses (CRPS, RMSE) and benchmark losses
    mcs.csv             # model confidence set inclusion for CRPS
  figures/              # output
```

## How to run

Use `replication/` as the working directory:

```
Rscript out_script.R
```

The switch `use.raw` at the top of the script selects the mode:

- `use.raw <- NA` (default): raw mode if both `raw/eval_out.rda` and `raw/mcs_out.rda` exist, CSV mode otherwise.
- `use.raw <- TRUE`: raw mode; stops if a raw file is missing.
- `use.raw <- FALSE`: CSV mode, even if the raw files are present.

**Raw mode** needs `raw/eval_out.rda` and `raw/mcs_out.rda`. It computes the losses from the source files and (over)writes `data/losses_avg.csv` and `data/mcs.csv`. Loading `eval_out.rda` takes about a minute and needs about 4 GB of memory.

**CSV mode** needs only `data/losses_avg.csv` and `data/mcs.csv`.

In both modes, the figures are produced from the two CSVs, so both modes run the same plotting code. The paths to the raw files can be changed in `raw.eval` and `raw.mcs`.

## Software

The output was checked against the figures of the paper with:

- R 4.6.1
- dplyr 1.2.1, tidyr 1.3.2, reshape2 1.4.5, lubridate 1.9.5
- ggplot2 4.0.3, ggh4x 0.3.1, scales 1.4.0, cowplot 1.2.0

In both modes, all eight files have the same text layer and render to the same image (Ghostscript, 150 dpi) as the versions in the paper.

## Data and upstream code

The raw source files are large, and thus only available upon request from the corresponding author.
