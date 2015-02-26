## plots censorhip ilustrative graphic
## refer to file survival-analysis/org/censorhip.org
## for reference data being plotted

## sets up work directory
work.dir <- "~/Dropbox/thesis/R/survival-analysis"
setwd (work.dir)

## loads fonts
library (extrafont)
loadfonts ()

## plots
svg ("img/censorship.svg",
     family = "OfficinaSansITC",
     height = 4,
     width = 4)
## patient A
plot (4, 5,
      pch = 19,
      xlab = "Time (years)",
      ylab = "Patient",
      yaxt = "n",
      xlim = c (0, 12),
      ylim = c (.5, 6.5))
segments (0, 5, 4, 5)
## patient B
points (11, 4, pch = 19)
segments (0, 4, 11, 4)
## patient C
segments (0, 3, 2, 3)
points (2, 3, pch = 19, col = "white")
points (2, 3, pch = 21)
## patient D
points (9, 2, pch = 19)
segments (0, 2, 9, 2)
## patient E
segments (0, 1, 6, 1)
points (6, 1, pch = 19, col = "white")
points (6, 1, pch = 21)
## redefines y-axis values
axis (2, at = 1:5, labels = LETTERS[5:1])
## adds grey vertical line at x = 0
abline (v = 0,
        lty = "dotted",
        col = "grey60")
## adds red vertical line at x = 10
abline (v = 10,
        lty = "dotted",
        col = "darkred")
dev.off ()
