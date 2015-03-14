prj.dir <- "~/Dropbox/thesis/R/datasets"
setwd (prj.dir)

dfr <- read.csv2 ("csv/studies.csv",
                  row.names = 1,
                  stringsAsFactors = FALSE)

cols2keep <- c ("cancer",
                "platform",
                "normalization",
                "nSamples",
                "nGenes",
                "os",
                "dfs",
                "dss",
                "dmfs")

tab.dfr <- dfr[cols2keep]

write.csv2 (tab.dfr,
            "csv/datasets-table.csv")
