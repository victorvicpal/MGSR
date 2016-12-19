## Multivariate Gaussian Subspatial Regression

This repository contains the required functions to apply the **MGSR**. 


----------
# Table of Contents
1. [Summary](#summary)
2. [Algorithm](#algorithm)
3. [Installation](#installation)
4. [Reference](#reference)
5. [Example](#example)


-------------
## Summary

Our proposal is a mixture between Factorial Techniques and Gaussian Processes.

Factorial techniques ([Biplot](https://www.wikiwand.com/en/Biplot), [PCA](https://www.wikiwand.com/en/Principal_component_analysis), [CA](https://www.wikiwand.com/en/Correspondence_analysis) ...) are useful to represent our observations in terms of unobserved variables called factors. 
All these techniques provide a set of coordinates linked to the observations, which display information of our analyzed variables. 
These kind of procedures are merely descriptive and have a low predictive power.

On the other hand, [Gaussian Processes](https://www.wikiwand.com/en/Gaussian_process) are statistical methods where observations occur in a continuous domain (mainly time or space). 
Furthermore, variables have a multivariate normal distribution. Gaussian Processes use similarity between points to predict the value in an unobserved point.

The **MGSR** applies a Gaussian Process Regression to the created subspatial of the Factorial Technique.
We make use of this subspace to simulate a continuous domain that permit the application of Gaussian Processes, such as Cokriging.

## Algorithm

1. Factorial Technique
2. [Cross-variogram](https://github.com/victorvicpal/MGSR/blob/master/crossvariogram.R)
    * [CV plot (cross variogram and LMC)](https://github.com/victorvicpal/MGSR/blob/master/plot.crossvariogram.R)
3. [Linear Model of Coregionalization](https://github.com/victorvicpal/MGSR/blob/master/lmc.R)
4. [Subspatial Grid](https://github.com/victorvicpal/MGSR/blob/master/grid.R)
5. [Cokriging](https://github.com/victorvicpal/MGSR/blob/master/cokrig.R)

> **Note:**

> - Cross-variograms normally follows a "Power Distribution" due to the small scale of Factorial Techniques.
> - Unlike geostatistical analyses, we don't have a real field where boundaries restrict our study. On the other hand, this aspect is more positive than negative because we can portray a more simple grid without losing information.

## Installation
```
install.packages('devtools')
library(devtools)
install_github("victorvicpal/MGSR")
library(MGSR)
```

## Reference
[![DOI](https://zenodo.org/badge/76850843.svg)](https://zenodo.org/badge/latestdoi/76850843)

## Example


-------------
