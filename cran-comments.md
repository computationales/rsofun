Dear CRAN team,

This is the re-submission of the {rsofun} package. We have removed non-default lower bound specification in derived type, in order to fix the compilation issues that appeared in the Debian check.

We were not able to reproduce the error shown by Prof. Ligges but have tried to solve the issue anyways, and our checks pass cleanly. If other issues arise or this one is still present, we would need to be provided with the configuration of the machine used to run the checks on CRAN's end. This way, we could reproduce the error and see if our changes get rid of it.

---

The {rsofun} package provides the implementation of a modelling framework for site-scale simulations of ecosystem processes, with low level routines in Fortran 90. It contains the following models:
- P-model for leaf-level acclimation of photosynthesis from Stocker et al. (2019).
- SPLASH for bioclimatic variables, including the surface radiation budget and the soil water balance from Davis et al. (2017).
- BiomeE for comprehensive simulations of ecosystem carbon and water cycling, tree growth, and tree cohort-explicit forest dynamics following the Perfect Plasticity Approximation, from Weng et al. (2015).

This package is an extension of {rpmodel} in the sense that it expands the P-model implementation and provides functions for multiple-site simulations and model parameter calibration. 

The full documentation can be found at the github repository link: https://geco-bern.github.io/rsofun

Code coverage sits at ~76%, with remaining uncovered code pertaining to minor input data format checks of the main functions. The underlying P-model implementation is based on the {rpmodel} and the parameter calibration routines use packages {GenSA} and {BayesianTools}.

I hope this package is useful for other earth system scientists and the larger CRAN community. Kind regards, Josefa Arán.

---

I have read and agree to the CRAN policies enumerated here: https://cran.r-project.org/web/packages/policies.html

## Local, github actions and r-hub checks

- Pop!_OS 22.04 install on R 4.3

- Ubuntu 22.04, MacOS and Windows on github actions (devel / release)

- rhub::check_on_cran() with only notes for latex elements

- codecove.io code coverage at ~76%

## Github actions R CMD check results

0 errors | 0 warnings | 0 notes
