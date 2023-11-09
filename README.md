
<!-- README.md is generated from README.Rmd. Please edit that file -->

# partialLDSC

<!--<img src="inst/Figures/logo.png" align="right" height=180/> -->

:information\_source: `partialLDSC` is still under active development.  
<!-- Check the [NEWS](NEWS.md) to learn more about what has been modified\! -->

## Overview

`partialLDSC` is an R-package to estimate partial genetic
correlations.  
Relies on cross-trait LD-score regression (LDSC)  
This package builds up on the
[`GenomicSEM`](https://github.com/GenomicSEM/GenomicSEM/) R-package:
modification of the LDSC implementation + should use their munge
function to process the data

There are three functions available:

-   **`ldsc_partial()`**  
    main function that performs LDSC to estimate unadjusted and partial
    genetic correlations+ heritability (observed scale only)

-   **`heatmap()`**  
    …

-   **`forest_plot()`**  
    …

More details about its usage can be found in the
[manual](doc/partialLDSC-manual.pdf).

## Installation

You can install the current version of `partialLDSC` with:

``` r
# Directly install the package from github
# install.packages("remotes")
remotes::install_github("GEMINI-multimorbidity/partialLDSC")
library(partialLDSC)
```

## Usage

To run the analysis with `MRlap` different inputs are needed:

#### 1. The munged GWAS summary statistics (`conditions` & `confounder`):

#### 2. The input files for LDSC (`ld`):

These are needed by the
[`GenomicSEM`](https://github.com/GenomicSEM/GenomicSEM/) R-package.

-   ld:

> Expects LD scores formated as required by the original LD score
> regression software. Weights for the european population can be
> obtained by downloading the eur\_w\_ld\_chr folder in the link below
> (Note that these are the same weights provided by the original
> developers of LDSC):
> <https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v>

### Analysis

Before running the examples, please make sure to have downloaded the
input files for LDSC. You may also need to modify the `ld` & `hm3`
parameters to indicate the correct paths.

-   **Example A**

-   **Example B**

### Results

**`ldsc_partial()`** returns a named list containing the following
results:

-   **Example A**

-   **Example B**

## Runtime

Example A \~

Example B \~

The runtime can be influenced by the number of traits.

<!-- <font color="grey"><small> Results from analyses performed on a MacBook Pro (2020) - Processor : 2 GHz Quad-Core Intel Core i5 - Memory : 16 GB 3733 MHz LPDDR4X.</font> </small>    -->

## Contributors

## Citation

If you use the `partialLDSC` package, please cite:

## Contact

<!-- <mounier.ninon@gmail.com> -->
