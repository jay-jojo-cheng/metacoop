  <h3 align="center">metacoop: Metaanalysis of Individualized Treatment Rules via Sign-Coherency</h3>

  <p align="center">
    Functions for performing metaanalysis of individualized treatment rules and other causal inference utilities. 
    <br />
    <br />
    <a href="https://github.com/jay-jojo-cheng/metacoop/issues">Report Bug</a>
    ·
    <a href="https://github.com/jay-jojo-cheng/metacoop/issues">Request Feature</a>
  </p>
</p>


### Description

This package implements efficient procedures for estimating a solution path for individual level meta-analysis of individualized treatment rules. It penalizes sign-inconsistent gradients within groups of similar variables. The base learners implemented are A learning and weighted learning. The package also contains utilities for crossfitting propensity scores and semiparametric efficiency augmentation functions, an adaptive information criterion for value and concordance based model selection, and functions for reproducing simulations from Cheng, Huling, and Chen 2022.


<!-- GETTING STARTED -->
### Installation

You can install this package by running the following commands:

``` r
if (!require(devtools)) install.packages('devtools')
devtools::install_github("jay-jojo-cheng/metacoop")
```


### References
\[1\] Cheng JJ, Huling J, Chen G (2022). Meta-analysis of Individualized Treatment Rules via Sign-Coherency
Proceedings of Machine Learning for Health (ML4H), 2022. In press.
