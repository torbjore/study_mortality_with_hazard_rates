# study_mortality_with_hazard_rates

Author: Torbjørn Ergon
Date: 9 April 2018

These files accompany the paper «Study mortality with hazard rates, not probabilities» - https://doi.org/10.1101/216739


* The directory ‘supporting_information’ contains the R Markdown file used to generate the Supporting Information containing 3 appendices.


* The directory ‘modified_marked_files’ contains modified files from the ‘marked’ R package by Laake, Johnson & Conn (2013). Files are modified to enable modelling of cause-specific mortality hazard rates:

- ‘multistate_csmhr.tpl’ (original ‘multistate.tpl’) is the modified model template file for ADMB (Fournier et al. 2012).

- ‘mscjs_csmhr.R’ (original mscjs.R) calls the executable generated from ‘multistate_csmhr.tpl’ by ADMB.

- crm_csmhr.R (original crm.R) is the model fitting function (calling ‘mscjs_csmhr.R’) to be used.

See Appendix S3 in the Supporting Information (the ‘supporting_information’ directory) for documentation and a worked example using simulated data.


REFERENCE:

Fournier, D.A., Skaug, H.J., Ancheta, J., Ianelli, J., Magnusson, A., Maunder, M.N., Nielsen, A. \& Sibert, J. (2012) AD Model Builder: using automatic differentiation for statistical inference of highly parameterized complex nonlinear models. Optim. Methods Softw., 27, 233-249.

Laake, J.L., Johnson, D.S. & Conn, P.B. (2013) marked: An R package for maximum-likelihood and MCMC analysis of capture-recapture data. Methods in Ecology and Evolution, 4, 885-890
