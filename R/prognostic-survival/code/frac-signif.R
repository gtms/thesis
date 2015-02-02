## Produces image with fraction of significant tests predicting outcome across
## all cancer datasets in the study

## libraries
library (ggplot2)
library (data.table)
library (plyr)
library (extrafont)
library (grid)

## setup work directory
work.dir <- "~/Dropbox/thesis/R/prognostic-survival"
setwd (work.dir)

## loads data
load ("raw-data/fracGenesSignif.dfr.RData")
load ("raw-data/fracSignif.dfr.RData")
studies.dfr <- read.csv2 ("raw-data/studies.csv",
                          stringsAsFactors = FALSE)

## preprocess data
## defines function to choose event, should there be more than one available
selectEvent <- function (dfr,
                         eventOrder = c ("dss", "os", "dmfs", "dfs")) {
    idx <- tapply (1:nrow (dfr),
                   dfr$study,
                   function (x) {
                       events <- dfr$event[x]
                       x[events == eventOrder[eventOrder %in% events][1]]
                   })
    dfr[idx, ]
}

## renames certain cancers
## atypical-teratoid-rhabdoid-tumors to teratoid-rhabdoid
studies.dfr$cancer[studies.dfr$cancer == "atypical-teratoid-rhabdoid-tumors"] <- "teratoid-rhabdoid"
studies.dfr$cancer[studies.dfr$cancer == "cervical-cancer"] <- "cervical"

## defines vars2keep
vars2keep <- c ("study", "cancer", "os", "dfs", "dss", "dmfs", "nGenes")

## produces fracSigs.dfr
fracSigs.dfr <- merge (fracSignif.dfr,
                       studies.dfr[, vars2keep])
fracSigs.dfr <- selectEvent (fracSigs.dfr)
fracSigs.dfr$eventSigs.rdx <- ifelse (fracSigs.dfr$event %in% c ("os", "dss"),
                                  "death", "relapse")

## produces fracGenes.dfr
fracGenes.dfr <- merge (fracGenesSignif.dfr[, c ("fracSignif.pval", "study", "event")],
                        studies.dfr[, vars2keep])
names (fracGenes.dfr)[names (fracGenes.dfr) == "fracSignif.pval"] <- "fracSignifGenes.pval"
fracGenes.dfr <- selectEvent (fracGenes.dfr)
fracGenes.dfr$eventGenes.rdx <- ifelse (fracGenes.dfr$event %in% c ("os", "dss"),
                                   "death", "relapse")

## merges fracSigs.dfr and fracGenes.dfr
## merged.dfr <- merge (fracSigs.dfr[, !names (fracSigs.dfr) %in% "event"],
merged.dfr <- merge (fracSigs.dfr,
                     fracGenes.dfr[, c ("study", "fracSignifGenes.pval", "eventGenes.rdx")],
                     by = "study")

## redefines nSamples according to the chosen event
merged.dfr$nSamples <- sapply (1:nrow (merged.dfr), function (n) {
    merged.dfr[n, merged.dfr[n, "event"]]
})

## reorders study according to fraction of MSigDBc2 signatures associated with outcome
merged.dfr$study <- with (merged.dfr,
                          ordered (study, study[order (random, decreasing = TRUE)]))
rownames (merged.dfr) <- merged.dfr$study
merged.dfr <- merged.dfr[levels (merged.dfr$study), ]

## defines gse variable
merged.dfr$gse <- with (merged.dfr,
                        ifelse (grepl ("GSE", study),
                                gsub ("(.*)-(.*)$", "\\1", study), as.character (study)))

## discretizes by endpoint
merged.dfr$msigdb.death <- ifelse (merged.dfr$eventSigs.rdx == "death", merged.dfr$msigdb, NA)
merged.dfr$msigdb.relapse <- ifelse (merged.dfr$eventSigs.rdx == "relapse", merged.dfr$msigdb, NA)
merged.dfr$random.death <- ifelse (merged.dfr$eventSigs.rdx == "death", merged.dfr$random, NA)
merged.dfr$random.relapse <- ifelse (merged.dfr$eventSigs.rdx == "relapse", merged.dfr$random, NA)
merged.dfr$genes.death <- ifelse (merged.dfr$eventGenes.rdx == "death", merged.dfr$fracSignifGenes.pval, NA)
merged.dfr$genes.relapse <- ifelse (merged.dfr$eventGenes.rdx == "relapse", merged.dfr$fracSignifGenes.pval, NA)

## defines nStudies variable
nStudies <- nrow (merged.dfr)

## defines color palette (The Economist)
## as taken from the object ggthemes_data$economist$stata$fg
## from the ggthemes package
colors <- c (blue_gray = "#6794a7",
             blue_dark = "#014d64",
             green_light = "#76c0c1",
             blue_mid = "#01a2d9",
             blue_light = "#7ad2f6",
             green_dark = "#00887d",
             gray = "#adadad",
             blue_light = "#7bd3f6",
             red_dark = "#7c260b",
             red_light = "#ee8f71",
             green_light = "#76c0c1",
             brown = "#a18376")

## colors in rgb
colorsRGB.mtx <- sapply (colors, col2rgb)

## defines layout of device; initializes pdf device
## opens device
svg (file = "img/fraction-significant-tests-base-R-graphics.svg",
     width = 12,
     height = 10)

## pdf (file = "img/fraction-significant-tests-base-R-graphics.pdf",
##      width = 12,
##      height = 10)

## sets layout
layout (mat = matrix (1:4, ncol = 4),
        widths = c (1.5, 1, .5, 1.5))

## defines the setup the numbers of columns and the margins of the device
## par (mar = c (4, 0, 5, .5))

## panel #1
## study names
## creates the empty plot canvas to accomodate output of text command
plot (x = NA,
      y = NA,
      xlim = c (0, 1), # defines 1 column
      ylim = c (0, nStudies),
      xaxt = "n", # no x axis
      yaxt = "n", # no y axis
      xlab = "", # no x label
      ylab = "", # no y label
      main = "Study", # title of study
      frame.plot = FALSE) # whether to give the plot a frame or not

## lays out the text on the canvas created by the plot command
text (x = .5,
      y = nStudies:1,
      labels = merged.dfr$gse,
      cex = .8, # sets the size of the labels
      ## pos = 2, # Values of '1', '2', '3' and '4', respectively indicate
                                        # positions below, to the left of, above
                                        # and to the right of the specified
                                        # coordinates.
      adj = 0,
      font = 1)

## panel #2
## cancer type
## par (mar = c (9.9, .1, 9.9, .1))
## creates the empty plot canvas to accomodate output of text command
plot (x = NA,
      y = NA,
      xlim = c (0, 1), # defines 1 column
      ylim = c (0, nStudies),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      main = "Cancer",
      frame.plot = FALSE)

## lays out the text on the canvas created by the plot command
text (x = .5,
      y = nStudies:1,
      labels = merged.dfr$cancer,
      cex = .8, # sets the size of the labels
      ## pos = 1, # Values of '1', '2', '3' and '4', respectively indicate
                                        # positions below, to the left of, above
                                        # and to the right of the specified
                                        # coordinates.
      adj = 0,
      font = 1)

## panel #3
## number of samples
## par (mar = c (9.9, .1, 9.9, .1))
## creates the empty plot canvas to accomodate output of text command
plot (x = NA,
      y = NA,
      xlim = c (0, 1), # defines 1 column
      ylim = c (0, nStudies),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      main = "n",
      frame.plot = FALSE)

## lays out the text on the canvas created by the plot command
text (x = .5,
      y = nStudies:1,
      labels = merged.dfr$nSamples,
      cex = .8, # sets the size of the labels
      ## pos = 2, # Values of '1', '2', '3' and '4', respectively indicate
                                        # positions below, to the left of, above
                                        # and to the right of the specified
                                        # coordinates.
      adj = 1,
      font = 1)

## panel #4
## This panel plots the fraction of significant tests
## for studies with DEATH as endpoint
## colours to be used: blue_mid and blue_dark
## par (mar = c (9.9, .1, 9.9, .1))
## sets the emply plot device
plot (x = NA,
      y = NA,
      ## col = colors["blue_mid"],
      pch = 19,
      main = "Fraction associated with outcome",
      xlim = c (0, 1),
      ylim = c (0, nStudies),
      yaxt = "n", ## no x axis
      xlab = "", ## no x label
      ylab = "", ## no y label
      frame.plot = FALSE) # whether to give the plot a frame or not

## adds horizontal faint guidelines
abline (h = 1:nStudies,
        lty = 3,
        col = "grey")

## adds red line at x = .05
abline (v = .05,
        lty = 3,
        col = "darkred")

## adds results with randomized signatures
points (x = rev (merged.dfr$random.death),
        y = 1:nStudies,
        col = colors["blue_mid"],
        pch = 19)

## adds results with regular signatures
points (x = rev (merged.dfr$msigdb.death),
        y = 1:nStudies,
        col = colors["blue_dark"],
        pch = 19)

## adds results with fraction of genes associated with outcome
points (x = rev (merged.dfr$genes.death),
        y = 1:nStudies,
        col = colors["brown"],
        pch = 19)

## adds results with randomized signatures
points (x = rev (merged.dfr$random.relapse),
        y = 1:nStudies,
        col = colors["blue_mid"],
        pch = 1)

## adds results with regular signatures
points (x = rev (merged.dfr$msigdb.relapse),
        y = 1:nStudies,
        col = colors["blue_dark"],
        pch = 1)

## adds results with fraction of genes associated with outcome
points (x = rev (merged.dfr$genes.relapse),
        y = 1:nStudies,
        col = colors["brown"],
        pch = 1)

## panel #5
## This panel plots the fraction of significant tests
## for studies with RELAPSE as endpoint
## colours to be used: blue_mid and blue_dark
## par (mar = c (9.9, .1, 9.9, .1))
## sets the emply plot device
## plot (x = NA,
##       y = NA,
##       ## col = colors["blue_mid"],
##       pch = 19,
##       main = "Relapse",
##       xlim = c (0, 1),
##       ylim = c (0, nStudies),
##       yaxt = "n", ## no x axis
##       xlab = "", ## no x label
##       ylab = "", ## no y label
##       frame.plot = FALSE) # whether to give the plot a frame or not

## adds red line at x = .05
## abline (v = .05,
##         lty = 3,
##         col = "darkred")

## adds results with randomized signatures
## points (x = rev (merged.dfr$random.relapse),
##         y = 1:nStudies,
##         col = colors["blue_mid"],
##         pch = 19)

## adds results with regular signatures
## points (x = rev (merged.dfr$msigdb.relapse),
##         y = 1:nStudies,
##         col = colors["blue_dark"],
##         pch = 19)

## adds results with fraction of genes associated with outcome
## points (x = rev (merged.dfr$genes.relapse),
##         y = 1:nStudies,
##         col = colors["red_light"],
##         pch = 19)

## adds legend text at the bottom right corner
## regular
## text (x = 1,
##       y = 4,
##       labels = "MSigDB C2",
##       pos = 2,
##       cex = 1.5,
##       col = colors["blue_dark"])

## random
## text (x = 1,
##       y = 1,
##       labels = "random",
##       pos = 2,
##       cex = 1.5,
##       col = colors["blue_mid"])

## genes
## text (x = 1,
##       y = 7,
##       labels = "singluar genes",
##       pos = 2,
##       cex = 1.5,
##       col = colors["brown"])

legend ("bottomright",
        legend = c ("death endpoint",
            "relapse endpoint",
            "single genes",
            "MSigDB C2 signatures",
            "random signatures"),
        bty = "n",
        lwd = 2,
        cex = 1.5,
        text.col = c (rep ("black", 2),
            colors["brown"],
            colors["blue_dark"],
            colors["blue_mid"]),
        lty = rep (NA, 5),
        pch = c (19, 1, rep (NA, 3)))

## closes device
dev.off ()
