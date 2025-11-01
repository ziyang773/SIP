# Autocorrelation-Test-under-Frequent-Mean-Shifts

Software accompanying the paper  
**"Autocorrelation Test under Frequent Mean Shifts"**  
by Liu, Z., Hao, N., Niu, Y. S., Xiao, H., & Ding, H. (2025).  [arXiv:2510.21047](https://arxiv.org/abs/2510.21047)

---

## Overview

This repository provides R code implementing the **Shift-Immune Portmanteau (SIP) Test** proposed in the above paper.  

Two main functions are provided:

- **`SIP.test`** — performs the autocorrelation test for a given order of dependence.  
- **`SIP.acf`** — computes and plots autocovariance or autocorrelation with confidence intervals for a given sequence.
- **`mean_sim`** — simulates sequences with frequently changed means.  
- **`acf_sim`** — simulates noise under different structures and different distributions.
---

## File Structure

| File | Description |
|------|--------------|
| `R/sip.R` | Contains main functions `SIP.test` and `SIP.acf`. |
| `R/helper.R` | Contains functions for estimating $W(\theta)$, autocovariance, and confidence intervals. |
| `R/simulation_helpers.R` | Simulation setup for Type I error and power analysis under different noise structures. |
| `Numerical/AJO242_4fc438f9_0.fast5` | Example nanopore sequencing dataset (requires the `rhdf5` package to read). |
| `Numerical/numerical_analysis.R` | Simulation and real-data analysis (nanopore sequencing data). |
  

---

## Installation and Usage

```r
# 1 install the package
devtools::install_github("ziyang773/SIP")
library(SIP)

# 2. Example usage
## 2.1 IID case
set.seed(1)
y <- rnorm(1000) + rep(c(rep(1, 50), rep(0, 50)), 10)

SIP.test(y)
Box.test(y) # for comparison with Box-Pierce test 

SIP.acf(y)  # defaults: lag.max = 4, type = "correlation", estimation = FALSE
acf(y) # for comparison with default acf

## 2.2 MA(1) alternative with mean shifts (theta = 0.3)
set.seed(1)
e <- rnorm(1001)
x <- e[-1] + 0.3 * e[-1001]
y <- x + rep(c(rep(1, 50), rep(0, 50)), 10)

SIP.test(y)
Box.test(y) # for comparison with Box-Pierce test 

SIP.acf(y)  # defaults: lag.max = 4, type = "correlation", estimation = FALSE
acf(y) # for comparison with default acf
```
---

## Notes

- For nanopore data analysis, install **`rhdf5`** from Bioconductor:
- The dataset AJO242_4fc438f9_0.fast5 is used as an example of raw nanopore sequencing input.

  ```r
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("rhdf5")


For questions, comments, or collaborations, please contact: ziyang773@arizona.edu
