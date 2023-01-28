# Addressing non-normality in multivariate analysis

[![R](https://img.shields.io/badge/Made%20with-R%20under%20development-success)](https://cran.r-project.org/)
[![MVT](https://img.shields.io/badge/MVT-0.3-orange)](http://mvt.mat.utfsm.cl/)
[![fastmatrix](https://img.shields.io/badge/fastmatrix-0.3--8196-orange)](https://faosorios.github.io/fastmatrix/)
[![fMultivar](https://img.shields.io/badge/fMultivar-4021.83-orange)](https://cran.r-project.org/package=fMultivar)
[![sn](https://img.shields.io/badge/sn-2.1.0-orange)](https://cran.r-project.org/package=sn)
[![plot3D](https://img.shields.io/badge/plot3D-1.4-orange)](https://cran.r-project.org/package=plot3D)
[![DOI](https://img.shields.io/badge/DOI-10.1007/s10182--022--00468--2-blue)](http://doi.org/10.1007/s10182-022-00468-2)

Supplementary material to **Addressing non-normality in multivariate analysis using the t-distribution** by Felipe Osorio, Manuel Galea, Claudio Henriquez and Reinaldo Arellano-Valle (AStA Advances in Statistical Analysis, DOI: [10.1007/s10182-022-00468-2](https://doi.org/10.1007/s10182-022-00468-2)).

Code written by: Felipe Osorio

Correspondence author: Felipe Osorio, Email: felipe.osorios@usm.cl

Code tested on:
- R under development (2018-02-21 r74285), running Linux Mint 18.3 (64 bits)
- R version 3.3.0, running OS X 10.13.4 (64 bits)
- R version 3.4.3, running Windows 10 (64 bits)

Attached packages: MVT 0.3, fastmatrix 0.3-8196, fMultivar 4021.83, sn 2.1.0, plot3D 1.4

CONTENTS:
- case_studies/Ex1_WindSpeed.R: R commands for the analysis of wind speed dataset (analyzed at Sections 5.2.1 from manuscript).
- case_studies/Ex2_PSG.R: R commands for the analysis of transient sleep disorder dataset (analyzed at Sections 5.2.2 from manuscript).
- code/boots.R: R functions for bootstrap and Negentropy.
- code/envelope.R: computation of QQ-plot with simulation envelopes.
- code/stats.R: R functions to obtain summary statistics.
- code/student.influence.R: R functions to compute the conformal curvature under the t-perturbation.
- app_B/App_B.R and mardia.R: R commands to construct the Mardia's skewness coefficient plot.
- README.md: this file.

DATASETS:
- PSG: clinical trial on transient sleep disorder dataset.
- WindSpeed: Wind speed in the Pacific North-West of the United States collected at three meteorological towers approximately located on a line and ordered from west to east: Goodnoe Hills (gh), Kennewick (kw), and Vansycle (vs).

PSG and WindSpeed now are provided by MVT package.
