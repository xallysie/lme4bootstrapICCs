# lme4bootstrapICCs
R script &amp; tutorial to boostrap 95% confidence intervals around intraclass correlation coefficients (ICCs) and variance components in linear mixed-effects models (based on lme4). Researchers may wish to compare ICCs to one another. This function modifies the lmeresample() package in R (which modifies the boot() function) to bootstrap 95% confidence intervals around ICC estimates and variance components. The 95% CIs of each ICC can then be examined for overlap to determine if they are different.
