## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=6, 
  fig.height=4
)
# Legge denne i YAML på toppen for å skrive ut til tex
#output: 
#  pdf_document: 
#    keep_tex: true
# Original:
#  rmarkdown::html_vignette:
#    toc: true

## ----setup--------------------------------------------------------------------
# Start the HDANOVA R package
library(HDANOVA)

## -----------------------------------------------------------------------------
# Load Candy data
data(candies)

# Fit ASCA model
mod <- hdanova(assessment ~ candy*assessor, data=candies)
summary(mod)

## -----------------------------------------------------------------------------
# Fit ASCA model
mod <- asca(assessment ~ candy*assessor, data=candies)
summary(mod)

## -----------------------------------------------------------------------------
# Fit ASCA model step by step
mod_hd <- hdanova(assessment ~ candy*assessor, data=candies)
mod_asca <- sca(mod_hd)

## -----------------------------------------------------------------------------
# Permutation testing (default = 1000 permutations, if not specified)
mod <- asca(assessment ~ candy*assessor, data=candies, permute=TRUE)
summary(mod)

