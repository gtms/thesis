## Computes differentiation signatures based on bodyMap ranked FPKM's

library (plyr)

work.dir <- "~/Dropbox/thesis/R/differentiation-signatures"
setwd (work.dir)

source ("code/extract-bodyMap-sig.R")

load ("raw-data/bodyMapRnk.Rda")
tissues <- colnames (bodyMapRnk.mtx)

differentiationSigs <- sapply (tissues,
                               extractBodyMapSig,
                               tol = 0,
                               prjDir = ".")

sigSize.dfr <- ldply (differentiationSigs, length)
names (sigSize.dfr) <- c ("tissue", "number of tissue-specific genes")

write.csv2 (sigSize.dfr, "out-data/sig-size.csv")
