## sets up work directory
work.dir <- "~/Dropbox/thesis/R/mds"
setwd (work.dir)

## loads fonts
library (extrafont)
loadfonts ()

euro.mds.2 <- cmdscale (eurodist, eig = T)
euro.mds.2

## items to remove
cities2rm <- c ("Hook of Holland",
                "Calais",
                "Cologne",
                "Cherbourg",
                "Lyons")

items2keep <- !rownames (euro.mds.2$points) %in% cities2rm

mds2plot <- euro.mds.2
mds2plot$points <- mds2plot$points[items2keep, ]
mds2plot$eig <- mds2plot$eig[items2keep]

## plot (euro.mds.2$points[, 1],
##       euro.mds.2$points[, 2],
##       type = "n",
##       xlab = "Coordinate 1",
##       ylab = "Coordinate 2",
##       xlim = c (-2500, 2500),
##       ylim = c (-2500, 2500))

## text (euro.mds.2$points[, 1],
##       euro.mds.2$points[, 2],
##       labels = labels (eurodist))

## euro.abb <- abbreviate (labels (eurodist))

## plot (euro.mds.2$points[, 1],
##       euro.mds.2$points[, 2],
##       type = "n",
##       xlab = "Coordinate 1",
##       ylab = "Coordinate 2",
##       xlim = c (-2500, 2500),
##       ylim = c (-2500,2500) )
## text (euro.mds.2$points[, 1],
##       euro.mds.2$points[, 2],
##       labels = euro.abb)

# Maybe rotating across the x-axis would produce a better reflection of reality!

## plots
svg ("img/mds.svg",
     family = "OfficinaSansITC",
     height = 5,
     width = 5)
plot (mds2plot$points[, 1],
      -mds2plot$points[, 2],
      type = "n",
      xlab = "",
      ylab = "",
      xlim = c (-2500, 2500),
      ylim = c (-2500, 2500))
text (mds2plot$points[, 1],
      -mds2plot$points[, 2],
      labels = labels (eurodist)[items2keep],
      cex = .6)
dev.off ()
