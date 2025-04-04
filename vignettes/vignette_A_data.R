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
# Find directory extdata from the multiblock package
mbdir <- system.file('extdata/', package = "multiblock")

# Comma separated values, row names in first column
meta_data <- read.csv(paste0(mbdir, "/meta_data.csv"), row.names = 1)
# If working directory matches file location:
# meta_data <- read.csv('meta_data.csv', row.names = 1)
meta_data

# Semi-colon separated values (locales where the decimal point is comma),
# no row names
proteins <- read.csv2(paste0(mbdir, "/proteins.csv"))
proteins

# Blank space separated data without labels
genes <- read.table(paste0(mbdir, "/genes.dat"))
genes

## -----------------------------------------------------------------------------
# Column-centring
genes_centred <- scale(genes, scale=FALSE)
colMeans(genes_centred) # Check mean values

# Autoscaling
genes_scaled <- scale(genes)
apply(genes_scaled, 2, sd) # Check standard deviations

# SNV (transpose, autoscale, re-transpose)
genes_snv <- t(scale(t(genes)))
apply(genes_snv, 1, sd) # Check standard deviations

## -----------------------------------------------------------------------------
# Default is sum-coding
dummycode(meta_data$colour)

# Treatment coding
dummycode(meta_data$colour, "contr.treatment")

# Full dummy-coding (rank deficient)
dummycode(meta_data$colour, drop = FALSE)

# Replace categorical with dummy-coded, use I() to index by common name
meta_data2 <- meta_data
meta_data2$colour <- I(dummycode(meta_data$colour, drop = FALSE))
meta_data2
meta_data2$colour

## -----------------------------------------------------------------------------
# Direct approach
blocks1 <- list(meta = meta_data2, proteins = proteins, genes = genes)

# Two-step approach
blocks2 <- list(meta_data2, proteins, genes)
names(blocks2) <- c('meta', 'proteins', 'genes')

# Same result
identical(blocks1, blocks2)

# Access by name or number
blocks1[['meta']]
blocks2[[1]]

## -----------------------------------------------------------------------------
# Construct block data.frame from list
df1 <- block.data.frame(blocks1)

# Construct block data.frame from data.frame:
# First merge blocks into data.frame
my_data <- cbind(meta_data2, proteins, genes)
# Then construct block data.frame using named 
# list of indexes
df2 <- block.data.frame(my_data, block_inds = 
        list(meta = 1:2, proteins = 3:5, genes = 6:8))

# Same result
identical(df1,df2)

# Access by name or number
df1[[2]]
df2[['proteins']]
df1[c(1,3)]
df1[-2]
df2[c('proteins','genes')]

# Use with formula interface (see other vignettes)
# sopls(meta ~ proteins + genes, data = df1)

# Use with single list interface (see other vignettes)
# mfa(df1[c(1,3)], ncomp = 3)

