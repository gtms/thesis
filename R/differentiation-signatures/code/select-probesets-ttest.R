## * import datasets
library (GenePattern)

roth.avg.gct <- read.gct ("/Users/giltomas/projects/differentiation-indices/raw-data/GSE3526/out-data/GSE3526-collapsed-frma.gct")
roth.navg.gct <- read.gct ("/Users/giltomas/projects/differentiation-indices/raw-data/GSE3526/out-data/GSE3526-frma.gct")
su.navg.gct <- read.gct ("/Users/giltomas/projects/illumina/raw-data/GSE96.gct")

renameColnamesRoth <- function (gct) {
  colnames (gct$data) <- gsub ("CNS", "brain", colnames (gct$data))
  colnames (gct$data) <- gsub ("lymph.nodes", "lymph", colnames (gct$data))
  colnames (gct$data) <- gsub ("colon.cecum", "colon", colnames (gct$data))
  colnames (gct$data) <- gsub ("adipose.+", "adipose", colnames (gct$data))
  gct
}

roth.navg.gct <- renameColnamesRoth (roth.navg.gct) # all tissue samples
roth.navg.noBrain.gct <- roth.navg.gct
roth.navg.noBrain.gct$data <- roth.navg.noBrain.gct$data[, colnames (roth.navg.noBrain.gct$data) != "brain"]

## * selectProbesetsTTest
selectProbesetsTTest <- function (sig,
                                  tissue,
                                  platform = "u133v2", # one of
                                        # "u133v2" or "u95av2"
                                  probs = .9) {
                                        # requires
                                        # roth.navg.noBrain.gct and
                                        # su.navg.gct
  stopifnot (platform %in% c ("u133v2", "u95av2"))
  ifelse (platform == "u133v2",
          gct <- roth.navg.noBrain.gct,
          gct <- su.navg.gct)
  
  getTopProbesets <- function (.df, tissue) {
    ttest <- sapply (.df$probeset, function (probeset) {
      ifelse (.df[probeset, "is.diagnostic"],
              t.test (gct$data[probeset, grepl (tissue, colnames (gct$data))],
                      gct$data[probeset, !grepl (tissue, colnames (gct$data))])$statistic,
              NA)
    })
    
    topProbeset <- function (ttest, .df) {
      sapply (tapply (ttest, .df$gene.name, function (x) x,
                      simplify = FALSE),
              function (y) {
                ifelse (sum (is.na (y)) == length (y),
                        NA,
                        names (y)[which (y == max (y, na.rm = TRUE))])
              })
    }
    
    data.frame (gene = levels (.df$gene.name),
                top.probeset = topProbeset (ttest, .df),
                stringsAsFactors = FALSE)
  }
  
  probesets <- rownames (gct$data)[gct$row.descriptions %in% sig]
  
  .df <- data.frame (probeset = probesets,
                     gene.name = factor (sapply (probesets, function (probeset) {
                       with (gct, row.descriptions[rownames (data) == probeset])
                     })),
                     is.diagnostic = sapply (probesets, function (probeset) {
                       median (gct$data[probeset, grepl (tissue, colnames (gct$data))],
                               na.rm = TRUE) > quantile (gct$data[probeset, ], probs = probs)
                     }),
                     stringsAsFactors = FALSE)
  
  list (all.probesets = .df[order (.df$gene.name, .df$probeset), ],
        top.probesets = getTopProbesets (.df, tissue))
}

## * getU133v2Probesets
getU133v2Probesets <- function (sig, tissue, probs = .9) {
  probesets.lst <- selectProbesetsTTest (sig = sig,
                                         tissue = tissue,
                                         platform = "u133v2",
                                         probs = probs)
  na.omit (probesets.lst$top.probesets$top.probeset)
}

## * getU95av2Probesets
getU95av2Probesets <- function (sig, tissue, probs = .9) {
  probesets.lst <- selectProbesetsTTest (sig = sig,
                                         tissue = tissue,
                                         platform = "u95av2",
                                         probs = probs)
  na.omit (probesets.lst$top.probesets$top.probeset)
}

## * selectProbesetsTTestRoth
selectProbesetsTTestRoth <- function (sig,
                                      tissue,
                                      probs = .9) {
                                        # requires
                                        # roth.navg.noBrain.gct
  getTopProbesets <- function (.df, tissue) {
    ttest.roth <- sapply (.df$probeset, function (probeset) {
      ifelse (.df[probeset, "is.diagnostic.roth"],
              t.test (roth.navg.noBrain.gct$data[probeset, grepl (tissue, colnames (roth.navg.noBrain.gct$data))],
                      roth.navg.noBrain.gct$data[probeset,
                                                 !grepl (tissue,
                                                         colnames (roth.navg.noBrain.gct$data))])$statistic,
              NA)
    })
    
    topProbeset <- function (ttest, .df) {
      sapply (tapply (ttest, .df$gene.name, function (x) x,
                      simplify = FALSE),
              function (y) {
                ifelse (sum (is.na (y)) == length (y),
                        NA,
                        names (y)[which (y == max (y, na.rm = TRUE))])
              })
    }
    
    data.frame (gene = levels (.df$gene.name),
                top.probeset.roth = topProbeset (ttest.roth, .df),
                stringsAsFactors = FALSE)
  }
  
  probesets <- rownames (roth.navg.noBrain.gct$data)[roth.navg.noBrain.gct$row.descriptions %in% sig]
  
  .df <- data.frame (probeset = probesets,
                     gene.name = factor (sapply (probesets, function (probeset) {
                       with (roth.navg.noBrain.gct, row.descriptions[rownames (data) == probeset])
                     })),
                     is.diagnostic.roth = sapply (probesets, function (probeset) {
                       median (roth.navg.noBrain.gct$data[probeset,
                                                          grepl (tissue,
                                                                 colnames (roth.navg.noBrain.gct$data))],
                               na.rm = TRUE) > quantile (roth.navg.noBrain.gct$data[probeset, ], probs = probs)
                     }),
                     stringsAsFactors = FALSE)
  
  list (all.probesets = .df[order (.df$gene.name, .df$probeset), ],
        top.probesets = getTopProbesets (.df, tissue))
}

## * selectProbesetsTTestSu95
selectProbesetsTTestSu95 <- function (sig,
                                      tissue,
                                      probs = .9) {
                                        # requires
                                        # non-tissue averaged and su.navg.noBrain.gct
  getTopProbesets <- function (.df, tissue) {
    ttest.su <- sapply (.df$probeset, function (probeset) {
      ifelse (.df[probeset, "is.diagnostic.su95"],
              t.test (su.navg.noBrain.gct$data[probeset, grepl (tissue, colnames (su.navg.noBrain.gct$data))],
                      su.navg.noBrain.gct$data[probeset,
                                               !grepl (tissue,
                                                       colnames (su.navg.noBrain.gct$data))])$statistic,
              NA)
    })
    
    topProbeset <- function (ttest, .df) {
      sapply (tapply (ttest, .df$gene.name, function (x) x,
                      simplify = FALSE),
              function (y) {
                ifelse (sum (is.na (y)) == length (y),
                        NA,
                        names (y)[which (y == max (y, na.rm = TRUE))])
              })
    }
    
    data.frame (gene = levels (.df$gene.name),
                top.probeset.su95 = topProbeset (ttest.su, .df),
                stringsAsFactors = FALSE)
  }
  
  probesets <- rownames (su.navg.noBrain.gct$data)[su.navg.noBrain.gct$row.descriptions %in% sig]
  
  .df <- data.frame (probeset = probesets,
                     gene.name = factor (sapply (probesets, function (probeset) {
                       with (su.navg.noBrain.gct, row.descriptions[rownames (data) == probeset])
                     })),
                     is.diagnostic.su95 = sapply (probesets, function (probeset) {
                       median (su.navg.noBrain.gct$data[probeset,
                                                        grepl (tissue,
                                                               colnames (su.navg.noBrain.gct$data))],
                               na.rm = TRUE) > quantile (su.navg.noBrain.gct$data[probeset, ], probs = probs)
                     }),
                     stringsAsFactors = FALSE)
  
  list (all.probesets = .df[order (.df$gene.name, .df$probeset), ],
               top.probesets = getTopProbesets (.df, tissue))
}
