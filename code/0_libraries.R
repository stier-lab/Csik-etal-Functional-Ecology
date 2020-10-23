##############################
# libraries
##############################

library(bbmle) # functional response analysis (in Bart's JAGS setup)
# library(boot) # bootstrapping
# library(broom) # used alonside nls.multstart for tidying up model outputs
library(car) # model predictions & general stats
library(chron) # working with dates 
library(corrplot) # correlation matrix
library(cowplot) # organizing fig panels & printing plots
# library(effsize) # calculate effect sizes 
# library(frair) # fitting FR
library(ggpubr) # (in Bart's JAGS setup)
library(ggridges) # ridgeline plots for displaying LTER temps 
library(ggExtra) # ggMarginal for correlation plots 
library(here) # because we love reproducibility 
# library(Hmisc) # plotting CI using stat_summary
# library(kableExtra) # knitting tables
library(lme4) # max consumption analysis (in Bart's JAGS setup)
library(lmerTest) # max consumption analysis (ended up using this rather than lme4)
library(lubridate) # tidyverse friendly dates 
library(MCMCvis) # monte carlo 
# library(minpack.lm) # nonlinear least squares for fitting FR
library(naniar) # dealing with missing data in LTER time series 
# library(nlme) # nonlinear mixed effects models for consumption data
# library(nls.multstart) # fitting FR
# library(nlstools) # used alongside nls.multstart for estimating CI of non-linear regression
# library(nls2) # nonlinear least squares with brute force for fitting FR
# library(patchwork) # arranging plots
library(R2jags) # JAGS analysis (in Bart's JAGS setup)
library(rjags) # JAGS analysis (in Bart's JAGS setup)
library(sjPlot) # create tables from model outputs
library(tidyverse) # the true hero - data wrangling & visualization (in Bart's JAGS setup)
library(tidybayes) # get credible intervals for functional responses

