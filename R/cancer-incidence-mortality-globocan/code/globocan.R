## Analysis of global estimated cancer mortality and incidence data by sex
## GLOBOCAN 2012

## libraries
library (ggplot2)
library (data.table)
library (plyr)
library (extrafont)
library (grid)

## setup work directory
work.dir <- "~/Dropbox/thesis/R/cancer-incidence-mortality-globocan/"
setwd (work.dir)

## loads fonts
## font_import ()
loadfonts ()

## loads data
vars.dfr <- expand.grid (data.type = c ("incidence", "mortality"),
                         sex = c ("males", "females"),
                         stringsAsFactors = FALSE)

frames.lst <- alply (vars.dfr, 1, function (x) {
    var <- sprintf ("%s-%s", x[1], x[2])
    file.name <- sprintf ("raw-data/globocan-cancer2012-%s.csv", var)
    dfr <- read.csv2 (file.name,
                      comment.char = "#",
                      stringsAsFactors = FALSE)
    dfr$data.type = unlist (x[1])
    dfr$sex = unlist (x[2])
    dfr
})

## cleans data
mlt.dfr <- do.call (rbind,
                    lapply (frames.lst, function (dfr) {
                        dfr[-1,
                            c ("Cancer", "ASR..W.", "data.type", "sex")]
                    }))
colnames (mlt.dfr) <- c ("cancer", "asr", "data.type", "sex")
mlt.dfr$asr <- as.numeric (mlt.dfr$asr)
mlt.dtb <- data.table (mlt.dfr)
mlt.dtb[, asr.sex := ifelse (sex == "males", -asr, asr)]

## select only most frequent types of cancer
top.cancers <- c ("Lung",
                  "Liver",
                  "Stomach",
                  "Colorectum",
                  "Breast",
                  "Bladder",
                  "Oesophagus",
                  "Kidney",
                  "Non-Hodgkin lymphoma",
                  "Leukaemia",
                  "Pancreas",
                  "Brain, nervous system",
                  "Melanoma",
                  "Thyroid",
                  "Prostate")
topCancers.dtb <- mlt.dtb[cancer %in% top.cancers]

## sort top cancers by their mortality
## top.cancers.srt <- topCancers.dtb[data.type == "mortality",
##                                   list (incidenceByCancer = asr[which.max (asr)]),
##                                   by = cancer][, cancer[order (incidenceByCancer)]]

top.cancers.srt <- topCancers.dtb[data.type == "mortality",
                                  list (mortalityByCancer = sum (asr)),
                                  by = cancer][, cancer[order (mortalityByCancer)]]

topCancers.dtb$cancer <- factor (topCancers.dtb$cancer,
                                 levels = top.cancers.srt,
                                 ordered = TRUE)
topCancers.dtb$y.pos <- as.numeric (topCancers.dtb$cancer) - .5

## computePolygonVertices
computePolygonVertices <- function (dtb = topCancers.dtb,
                                    c, d, s) { # respectively, cancer, data.type and sex
    ifelse (d == "incidence", y.move <- .35, y.move <- -.35)
    lst <- dtb[cancer == c & data.type == d & sex == s,
               list (x.pos = asr.sex,
                     y.pos = y.pos)]
    if (length (lst$x.pos) == 0) {
        dfr <- data.frame (x = rep (0, 4),
                           y = rep (0, 4))
    } else {
        dfr <- data.frame (x = c (0, lst$x.pos,
                               lst$x.pos, 0),
                           y = c (lst$y.pos,
                               lst$y.pos,
                               lst$y.pos + y.move,
                               lst$y.pos + y.move))
    }
    dfr
}

## builds polygon.dfr
loopVars.dfr <- expand.grid (c = unique (topCancers.dtb$cancer),
                             d = c ("mortality", "incidence"),
                             s = c ("females", "males"))
polygon.dfr <- mdply (loopVars.dfr, computePolygonVertices)
polygon.dfr$colour <- with (polygon.dfr, sprintf ("%s-%s", d, s))

## ggplots
key.colours <- c ("incidence-females" = "orange",
                  "mortality-females" = "orangered",
                  "incidence-males" = "skyblue",
                  "mortality-males" = "steelblue")

globocan.plt <- ggplot (polygon.dfr,
                        aes (x = x,
                             y = y,
                             fill = colour)) +
                                 geom_polygon ()

globocan.plt <- globocan.plt +
    scale_x_continuous (limits = c (-45, 45),
                        breaks = seq (-40, 40, 10),
                        labels = as.character (abs (seq (-40, 40, 10))))
globocan.plt <- globocan.plt +
    scale_fill_manual (values = key.colours,
                       breaks = c ("incidence-males",
                           "mortality-males",
                           "incidence-females",
                           "mortality-females"),
                       labels = c ("male\nincidence",
                           "male\nmortality",
                           "female\nincidence",
                           "female\nmortality"))
globocan.plt <- globocan.plt +
    scale_y_continuous (breaks = seq (.5, 13.5, 1),
                        labels = top.cancers.srt) +

                            geom_vline (xintercept = 0)

globocan.plt +
    theme_bw () +
        ## theme (text = element_text (family = "OfficinaSansITC", size = 14),
        theme (text = element_text (family = "OfficinaSansITC-Book", size = 20),
               axis.title.x = element_blank (),
               axis.title.y = element_blank (),
               axis.ticks = element_blank (),
               panel.grid.major = element_line (colour = "grey60"),
               panel.grid.major.y = element_blank (),
               panel.grid.minor = element_blank (),
               legend.position = "top",
               legend.title = element_blank (),
               panel.border = element_blank (),
               panel.background = element_blank (),
               plot.margin = unit (c (-0.5, 1, 0, 0.5), "lines"))

ggsave (width = 13,
        height = 7,
        "img/globocan.svg")
