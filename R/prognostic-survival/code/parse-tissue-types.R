## Produces table with tissues of origin for cancers analyzed

## setup work directory
work.dir <- "~/Dropbox/thesis/R/prognostic-survival"
setwd (work.dir)

library (data.table)

## loads studies.csv
studies.dfr <- read.csv2 ("raw-data/studies.csv",
                          stringsAsFactors = FALSE)

## correspondence file
emptyCrs.dfr <- data.frame (cancer = unique (studies.dfr$cancer),
                       tissue = NA)

write.csv2 (emptyCrs.dfr,
            file = "csv/empty-correspondence.csv")

crs.dfr <- read.csv2 ("csv/filled-correspondence.csv",
                      row.names = 1,
                      stringsAsFactors = FALSE)

dfr <- studies.dfr[, "cancer", drop = FALSE]
dfr$tissue <- NA

dtb <- data.table (dfr, key = "cancer")
crs.dtb <- data.table (crs.dfr, key = "cancer")

tab.dfr <- as.data.frame (table (studies.dfr$cancer))
names (tab.dfr)[1] <- "cancer"

tab.dtb <- data.table (tab.dfr, key = "cancer")

mgd.dtb <- crs.dtb[tab.dtb]

out.dtb <- mgd.dtb[, sum (Freq), by = tissue][order (V1, decreasing = TRUE)]

write.csv2 (out.dtb,
            file = "csv/cancer-studies-by-tissue-type.csv")
