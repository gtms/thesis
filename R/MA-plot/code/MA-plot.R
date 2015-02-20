## This script generates an image with two MA plots: one prior to normalization
## and another post-normalization

## * Preamble
library (affyPLM)
library (affydata)

## sets up work directory
work.dir <- "~/Dropbox/thesis/R/MA-plot/"
setwd (work.dir)

## * Plots
## sets display
par.bak <- par ()
par (mfrow = c (2, 1))

data(Dilution)
rma.dilution <- rma(Dilution)

svg ("img/ma-plot.svg",
     width = 4, height = 9)
par (mfrow = c (2, 1))
MAplot (Dilution, which = 1, plot.method = "smoothScatter", cex = 2)
MAplot (rma.dilution, which = 1, plot.method = "smoothScatter")
dev.off ()
