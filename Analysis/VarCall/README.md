# VarCall

This repository provides the VarCall model tool used in the analysis of high throughput CRISPR based saturated genome editing (SGE) data, one type of multiplexed assays of variant effects (MAVEs). It also includes related bioinformatics processing and statistical analysis code. The VarCall model mainly uses R and the R package rjags to perform the Bayesian two Gaussian components modeling.

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Test and Run](#Test-and-run)
- [License](#license)

# Overview
``VarCall`` aims to provide a comprehensive statistical model for SGE data analysis, including the modeling of batch effects, and the prediction of variant effects. The package utilizes a Bayesian hierachical two Gaussian components modeling. The package should be able to run on all major platforms (e.g. BSD, GNU/Linux, OS X, Windows) as long as R and related R Packages are installed.

# System Requirements
## Hardware requirements
`VarCall` package requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
This package is supported for *Linux*, but in principle should be able to run on all other platforms such as OS X and Windows. The package has been tested on the following systems:
+ Linux: Red Hat Enterprise Linux 8.8
+ Windows: Windows 10 Enterprise

### Package Dependencies
`VarCall` mainly depends on JAGS and the following R packages.

```
knitr, rjags, R2WinBUGS, R2jags, mgcv
```

# Installation Guide:
### Download or clone from github
```
git clone https://github.com/najiemayo/Couch_SGE_RAD51D/Analysis/VarCall/
```

### Install packages
First install JAGS in the system following the link here: https://mcmc-jags.sourceforge.io/
Then within R, type the following:
```
install.packages(c("knitr", "rjags", "R2WinBUGS", "R2jags", "mgcv"))
```
This should be done within one miniute.

# Run:
- To run the analysis:
  - `cd VarCall` make sure you are in this directory, and make sure the subfolder cache and figs are empty
  - Start R by type R in the command line
  - Within R, type the following `library(knitr); knit("RAD51Dmave25.noPosNorm.ldaErrVarModel.ClinVarExclude.Rtex")`

- Generate the pdf file with code and explanations:
  - In Linux, type `pdflatex RAD51Dmave25.noPosNorm.ldaErrVarModel.ClinVarExclude.tex`.
  - In Windows, use an installed Tex systeme, such as TexStudio to open `RAD51Dmave25.noPosNorm.ldaErrVarModel.ClinVarExclude.tex` and build the pdf file `RAD51Dmave25.noPosNorm.ldaErrVarModel.ClinVarExclude.pdf`
    
- Specified prior:
  - Currently the prior is set to a mean value of 0.2 using a beta distribution Beta(2, 8). The change the prior, modifiy the line `beta.a<-2.0` and `beta.b<-8.0`. 

- Expected output and running time.
The output files will be `RAD51Dmave25.noPosNorm.ldaErrVarModel.ClinVarExclude.tex`, `MAVEpostProbs.csv` and several pdf plots. File `RAD51Dmave25.noPosNorm.ldaErrVarModel.ClinVarExclude.tex` and the plots can be further compiled into a pdf file `RAD51Dmave25.noPosNorm.ldaErrVarModel.ClinVarExclude.pdf` with code, running result and code explanations. File `MAVEpostProbs.csv` is the main output file for predicted probabilities of being pathogenic based on the training labels. The following shows the main columns of output:
  - PrDel: the probability of being pathogenic
  - lPostOdds: the log posterior odds
  - logBF: the log Bayes factor
  - eta: the estimated effect size
  - eta.ll, eta.ul: 95% lower bound and upper bound of eta
The running time for the full data set takes ~12 hours (see R sessionInfo at the end of this page). Due to the stochastic nature of the model, there might be some minor differences in the results, but it should converge when the iteration is long enough.

# License

## R sessionInfo()

> sessionInfo()
R version 4.5.1 (2025-06-13)
Platform: x86_64-pc-linux-gnu
Running under: Rocky Linux 8.10 (Green Obsidian)

Matrix products: default
BLAS:   /usr/lib64/libblas.so.3.8.0 
LAPACK: /usr/lib64/liblapack.so.3.8.0  LAPACK version 3.8.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=C              
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/Chicago
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] mgcv_1.9-4       nlme_3.1-169     R2jags_0.8-9     R2WinBUGS_2.1-24 boot_1.3-32      rjags_4-17      
[7] coda_0.19-4.1    knitr_1.51      

loaded via a namespace (and not attached):
 [1] vctrs_0.7.2       cli_3.6.6         rlang_1.2.0       xfun_0.57         stringi_1.8.7     otel_0.2.0       
 [7] highr_0.12        glue_1.8.1        grid_4.5.1        evaluate_1.0.5    abind_1.4-8       lifecycle_1.0.5  
[13] stringr_1.6.0     compiler_4.5.1    codetools_0.2-20  rstudioapi_0.18.0 lattice_0.22-9    digest_0.6.39    
[19] pillar_1.11.1     parallel_4.5.1    splines_4.5.1     magrittr_2.0.5    Matrix_1.7-5      tools_4.5.1 


This project is covered under the **MIT License**.


