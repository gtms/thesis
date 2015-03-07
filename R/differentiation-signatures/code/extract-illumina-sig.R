## extractIlluminaSig
extractIlluminaSig <- function (tss,
                                top.tss = 1000,
                                top.other = 5000,
                                tol = 3) {
  prj.dir <- "/Users/giltomas/projects/illumina"
  load (file.path (prj.dir, "raw-data/illumina-list.Rda"))
  
  all.tss <- names (illumina)
  all.other.tss <- all.tss[all.tss != tss]
  
  illumina.coding <- lapply (names (illumina), function (tss) {
    subset (illumina[[tss]], source == "protein_coding" & chr != "chrM")
  })
  illumina.coding <- setNames (illumina.coding, names (illumina))
  
  tss.sig <- illumina.coding[[tss]]$gene.name[1:top.tss]
  
  tss.sig[sapply (tss.sig, function (gene) {
    sum (sapply (all.other.tss, function (tissue) {
      gene %in% illumina.coding[[tissue]]$gene.name[1:top.other]
    })) <= tol
  })]
}

extractIlluminaSig2 <- function (tss,
                                 top.tss = 1000,
                                 top.other = 5000,
                                 tol = 0) {
  prj.dir <- "/Users/giltomas/projects/illumina"
  load (file.path (prj.dir, "raw-data/illumina-rank.Rda"))

  all.tss <- colnames (rnk.mtx)
  all.other.tss <- all.tss[all.tss != tss]

  x <- rnk.mtx[, tss] <= top.tss
  y <- apply (rnk.mtx[, all.other.tss], 1, function (z) {
    sum (z < top.other) <= tol
  })

  rownames (rnk.mtx)[x & y]
}

extractRankDiff <- function (tss,
                             tol = 1) {
  prj.dir <- "/Users/giltomas/projects/illumina"
  load (file.path (prj.dir, "raw-data/ranked-illumina.Rda"))
  source (file.path (prj.dir, "code/utility-functions.R"))
  
  all.tss <- colnames (ranked.illumina)
  all.other.tss <- all.tss[all.tss != tss]
  
  tss.mtx <- ranked.illumina[, tss, drop = FALSE]
  all.other.mtx <- t (ranked.illumina[, all.other.tss])
  
  norm.mtx <- sapply (colnames (all.other.mtx), function (gene) tss.mtx[gene, 1] / all.other.mtx[, gene])
  ttest.stat <- sapply (as.data.frame (norm.mtx), function (gene.norm) {
    t.test (gene.norm,
            mu = 1,
            alternative = "greater")$statistic
  })
}
