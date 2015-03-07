## Extracts Illumina signatures for range of parameters

## * preamble
prj.dir <- "/Users/giltomas/projects/illumina"

library (plyr)
library (doMC)
library (CellMix)

source (file.path (prj.dir, "code/extract-illumina-sig.R"))
source (file.path (prj.dir, "code/select-probesets-ttest.R"))

load (file.path (prj.dir, "raw-data/illumina-rank.Rda"))

## registers parallel backend with the 'foreach' package 
registerDoMC (16)

## * generates sigs
get.sigs <- function(rk, h,l,t) {
  sigs <- sapply(1:16, function(x) get.sig.illumina(rk, x, h,l,t))
  names(sigs) <- colnames(rk)
  sigs <- sigs[0!=sapply(sigs, length)]

  sigs.u133v2 <- sapply(names(sigs), function(x) get.u133v2.probes(sigs[x], x))
  sigs.u133v2 <- sigs.u133v2[0!=sapply(sigs.u133v2, length)]

  sigs.u95av2 <- sapply(names(sigs), function(x) get.u95av2.probes(sigs[x], x))
  sigs.u95av2 <- sigs.u95av2[0!=sapply(sigs.u95av2, length)]

  list(genes=sigs, u133v2=sigs.u133v2, u95av2=sigs.u95av2)
}

getSigs <- function (top.tss, top.other, tol) {
  sigs <- llply (colnames (rnk.mtx), function (tss) {
    extractIlluminaSig2 (tss = tss,
                         top.tss = top.tss,
                         top.other = top.other,
                         tol = tol)
  },
                 .progress = "text",
                 .parallel = TRUE)
  names (sigs) <- colnames (rnk.mtx)
  sigs <- sigs[sapply (sigs, length) != 0]
  
  sigs.u133v2 <- sapply (names (sigs), function (tss) {
    getU133v2Probesets (unlist (sigs[tss]), tss)
  })
  sigs.u133v2 <- sigs.u133v2[sapply (sigs.u133v2, length) != 0]

  sigs.u95av2 <- sapply (names (sigs), function (tss) {
    getU95av2Probesets (unlist (sigs[tss]), tss)
  })
  sigs.u95av2 <- sigs.u95av2[sapply (sigs.u95av2, length) != 0]
  
  list (genes = sigs, u133v2 = sigs.u133v2, u95av2 = sigs.u95av2)
}

sigs.top1000.other5000.tol0 <- getSigs (1000, 5000, 0)
sigs.top1000.other5000.tol1 <- getSigs (1000, 5000, 1)
sigs.top5000.other7000.tol1 <- getSigs (5000, 7000, 1)

save (sigs.top1000.other5000.tol0,
      sigs.top1000.other5000.tol1,
      sigs.top5000.other7000.tol1,
      file = file.path (prj.dir, "raw-data/ref-sigs.Rda"),
      compress = TRUE)
