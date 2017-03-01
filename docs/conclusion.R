## ----setup, include=FALSE, echo=F, warning= F, message=F-----------------
knitr::opts_chunk$set(message = FALSE, 
                      warning = FALSE, 
                      error = FALSE, 
                      tidy = FALSE,
                      results = 'asis',
                      fig.align="center", 
                      dpi = 300, 
                      cache = FALSE,
                      progress=FALSE, 
                      quite = TRUE)

## ---- echo=FALSE, results="hide", warning=FALSE,message=FALSE------------
suppressPackageStartupMessages({
  if("TCGAWorkflow" %in% installed.packages()[,1]) {
    library(TCGAWorkflow)
  } else {
    devtools::load_all(".")
  }
})

## ------------------------------------------------------------------------
pander::pander(sessionInfo(), compact = FALSE)

