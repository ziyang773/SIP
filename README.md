# Autocorrelation-Test-under-Frequent-Mean-Shifts

Software accompanying the paper  
**"Autocorrelation Test under Frequent Mean Shifts"**  
by Liu, Z., Hao, N., Niu, Y. S., Xiao, H., & Ding, H. (2025).  [arXiv:2510.21047]

---

## Overview

This repository provides R code implementing the **Shift-Immune Portmanteau (SIP) Test** proposed in the above paper.  

Two main functions are provided:

- **`SIP.test`** — performs the autocorrelation test for a given order of dependence.  
- **`SIP.acf`** — computes and plots autocovariance or autocorrelation with confidence intervals for a given sequence.

---

## File Structure

| File | Description |
|------|--------------|
| `main.R` | Contains main functions `SIP.test` and `SIP.acf`. |
| `helper.R` | Contains functions for estimating $W(\theta)$, autocovariance, and confidence intervals. |
| `simulation_helpers.R` | Simulation setup for Type I error and power analysis under different noise structures. |
| `numerical_analysis.R` | Simulation and real-data analysis (nanopore sequencing data). |
| `AJO242_4fc438f9_0.fast5` | Example nanopore sequencing dataset (requires the `rhdf5` package to read). |

---

## Installation and Usage

```r
# 1. Clone or download this repository
#    (in terminal or Git Bash)
# git clone https://github.com/ziyang773/Autocorrelation-Test-under-Frequent-Mean-Shifts.git
# setwd("Autocorrelation-Test-under-Frequent-Mean-Shifts")

# 2. Open R and source the main script
source("main.R")

# 3. Read the documentation inside main.R for parameter descriptions and usage notes.

# 4. Example usage
set.seed(111)
x <- rnorm(100) + c(rep(1, 50), rep(0, 50))

# Perform SIP test for dependence order 4
SIP.test(x, m = 4)

# Plot autocorrelation with common CI width
SIP.acf(x, lag.max = 4, type = "correlation", estimation = FALSE, plot = TRUE)
```
---

## Notes

- All helper scripts (`helper.R`, `simulation_helpers.R`, `numerical_analysis.R`) must be in the same working directory.  
- For nanopore data analysis, install **`rhdf5`** from Bioconductor:
- The dataset AJO242_4fc438f9_0.fast5 is used as an example of raw nanopore sequencing input.

  ```r
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("rhdf5")


For questions, comments, or collaborations, please contact: ziyang773@arizona.edu
