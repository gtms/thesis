## Analysis of publication trends for MeSH (Medical Subject Headings) keywords
## "microarray" and "cancer"

## libraries
library (ggplot2)
library (ggthemes)
library (plyr)
library (extrafont)
library (grid)

## setup work directory
work.dir <- "~/Dropbox/thesis/R/publication-trends/"
setwd (work.dir)

## loads fonts
## font_import ()
loadfonts ()

## loads data
dataSources <- c ("gopubmed", "ligercat")

frames.lst <- setNames (llply (dataSources, function (source) {
    file.name <- sprintf ("raw-data/%s.csv", source)
    dfr <- read.csv2 (file.name,
                      comment.char = "#",
                      stringsAsFactors = FALSE)
}),
                        dataSources)

## ggplots
pubtrends.plt <- ggplot (frames.lst$ligercat,
                         aes (x = year, y = count)) +
                             geom_line (size = 2, colour = "#00526D") +
                                 scale_x_continuous (breaks = seq (1996, 2012, 1))

pubtrends.plt +
    theme_bw () +
        theme (text = element_text (family = "OfficinaSansITC-Book", size = 20),
               axis.title.x = element_blank (),
               axis.title.y = element_blank (),
               axis.ticks = element_blank (),
               panel.grid.major = element_line (colour = "grey60"),
               panel.grid.major.x = element_blank (),
               panel.grid.minor = element_blank (),
               legend.position = "top",
               legend.title = element_blank (),
               panel.border = element_blank (),
               panel.background = element_blank (),
               plot.margin = unit (c (-0.5, 1, 0, 0.5), "lines"))

ggsave (width = 13,
        height = 5,
        "img/pubtrends.pdf")
