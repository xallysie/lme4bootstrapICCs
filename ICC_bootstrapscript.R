# Bootstrapping Confidence Intervals around ICCs [osf.io/anwx2/]
# Updated 17/01/2020 for lme4 version 1.1-21 and boot version 1.3-24
# updated 12/02/2020 to be able to accommodate any number of clusters and their interactions
# updated 11/17/2021 - the latest update (01/05/2021) to the lmeresampler package basically makes this function obsolete, hooray!
#                    - it is now a convenience function to make the output easier to read

# library
library(boot)     # bootstrapping
library(psych)    # descriptives
library(lme4)     # fitting mixed-effects models
library(lmeresampler) # resampling mixed-effects models
library(tidyr)    # data hygiene

## uncomment to install the lmeresampler package from github
## if(!require(devtools)) install.packages("devtools",dependencies=TRUE, repos='http://cran.rstudio.com/')
## devtools::install_github("aloy/lmeresampler")

# import sample dataset 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df <- read.csv("sample_df.csv")
df <- df[df$WhatWasRating=="agg",] # look at ratings of aggressiveness

# This script will estimate intraclass correlation coefficients (ICC) for each cluster in a model with k clusters. 
# It can handle models that are cross-classified, as well as interactions between clusters.
# All ICC formulas based on Raudenbush & Bryk (2002).

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# Step 1. Fit your lmer model ####

# In this example the DV is cross-classified across Participants (ParticipantID) and Targets (StimName).
# DV = Rating
# cluster 1 variable     = ParticipantID
# cluster 2 variable     = StimName
# interaction component  = ParticipantID:StimName

# model <- lmer(Rating ~ 1 + (1 | ParticipantID), data=df)                  # sample model with 1 cluster
# model <- lmer(Rating ~ 1 + (1 | ParticipantID) + (1 | StimName), data=df) # sample model with 2 clusters (but no interaction)
# model <- lmer(Rating ~ 1 + (1 | ParticipantID) + (1 | StimName) + (1 | ParticipantID:StimName), data=df) # sample model with 2 clusters & an interaction

# For this demonstration, we will build a model with 2 clusters (but no interaction)
model <- lmer(Rating ~ 1 + (1 | ParticipantID) + (1 | StimName), data=df)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# Step 2. Load the bootICC function (requires lmeresampler package)

bootstrapICC <- function(model,             # Lmer model object
                         iterations = 5000, # Number of iterations or bootstrap resamples
                         boot.type,         # A character string indicating the type of bootstrap that is being requested,
                         ## e.g., "parametric", "residual", "case", "cgr", "reb" (see ?lmeresampler::bootstrap for more info)
                         boot.ci.type,      # Enter a character string indicating the type of interval that is being requested,
                         ## e.g., "all", "norm", "basic", "stud", "perc", "bca" (see ?lmeresampler::confint for more info)
                         ...){             
  m <- model
  ranefs = names(getME(m,"cnms")) # get cluster names
  ranef_k = length(ranefs)   # number of clusters
  require(lmeresampler)
  print(paste("Starting lmer resampling at",iterations,"iterations using",boot.type,"bootstrap and requesting",boot.ci.type,"intervals",sep=" "))
  
  mysumm   <- function(m){  #:::::::::::::# function returning the outputs/statistics of interest
    sig    <- getME(m, "sigma")           # sigma
    R_ij   <- sig^2; names(R_ij) <- "Residual" 
    taus   <- ((getME(m,"theta")*sig)^2)  # variance components for k clusters
    names(taus) <- gsub("\\s*\\.\\([^\\)]+\\)\\s*$","",names(taus)) # clean up column names
    taus   <- taus[c(ranefs)]             # re-order according to original cluster names
    vars   <- c(taus,R_ij)
    totvar <- sum(vars)                   # Total Variance = Sum of all variance components including residual variance (i.e., sigma^2)
    ICCs   <- sapply(c(vars),
                     function(x){x/totvar},   # Calculate ICC for each cluster, including the residual ICC
                     simplify=T, USE.NAMES=T) # (Raudenbush & Bryk, 2002)
    names(sig) <- "(Sigma)"
    names(vars)<- paste0(names(vars), " Variance")
    names(ICCs)<- paste0(names(ICCs), " ICC")
    c(sig,vars,ICCs)
  }
  # store summary object to call later on
  pointsummary = mysumm(m)
  
  # bootstrap your lmer object
  bootie <- lmeresampler::bootstrap(model = m, .f=mysumm, type=boot.type, B=iterations)
  
  # store bootstrapped samples in dataframe
  #dfbootie <- as.data.frame(bootie[['replicates']])
  #names(dfbootie) <- names(pointsummary)

  # what to extract from boot.ci (don't touch this unless boot.ci got updated and broke things)
  #boot.ci.index <- as.numeric(c(1:(2*ranef_k+3)))
  CIs <- confint(bootie, level = 0.95, type = boot.ci.type)
  
  # store 95% CIs
  #CIs     <- lapply(boot.ci.index, function(i){boot.ci(bootie, index=i, type=boot.ci.type)})
  #CIs.out <- sapply(CIs, function(j){tail(as.numeric(j[[4]]),2)}) 
  #CIs.out <- as.data.frame(t(CIs.out))
  #names(CIs.out) <- c("l-95% CI", "u-95% CI")
  #row.names(CIs.out) <- c(names(pointsummary[1:length(pointsummary)]))
  #pointsummary <- as.data.frame(pointsummary); names(pointsummary) = "Estimate"
  #CIs.out <- cbind(CIs.out, pointsummary)
  
  # split variances and ICCs into diff columns
  #CIs.ICC <- CIs.out[grepl("ICC", rownames(CIs.out))==TRUE,]
  #CIs.var <- CIs.out[grepl("ICC", rownames(CIs.out))==FALSE,]
  
  # print CI output
  #print(CIs.var,digits=5); print(CIs.ICC,digits=5)
  
  # return list containing boot object, boot.ci output, and the table of 95% CIs with their Point Estimates
  #merp <- list(bootie, dfbootie, CIs, CIs.var, CIs.ICC)
  #names(merp) <- c("boot_object", "df_bootsamples", "bootCI_object", "CI_table_Variance", "CI_table_ICC")
  
  print(CIs, digits=5)
  merp = list(bootie,CIs)
  names(merp) = c("boot_object","CIs")
  print("Done! Have a great day.")
  return(merp)
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# Step 3. Adjust parameters below and run the bootstrapICC function

# model         = Lmer model object
# iterations    = Number of iterations or bootstrap resamples (default:5000)
# boot.type     = A character string indicating the type of bootstrap that is being requested (e.g., "parametric", "residual", "case", "cgr", "reb")
#                 (see ?lmeresampler::bootstrap for more info)
# boot.ci.type  = A character string indicating the type of interval that is being requested (e.g., "norm", "basic", "stud", "perc", "bca")
#                 (see ?boot.ci for more info)

booty <- bootstrapICC(model, iterations=100, boot.type="parametric", boot.ci.type="norm")

# The output is a list containing:
View(booty[["boot_object"]]) # the original boot object from lmeresampler::bootstrap()
booty[["CIs"]]               # dataframe of 95% Interval Estimates of Variance Components and intraclass correlation coefficients 
 
