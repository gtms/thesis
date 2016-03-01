setwd("~/ownCloud/projects/prognostic-survival")

library (ProjectTemplate)
load.project ()

## require(RCurl)
## source(textConnection(getURL("https://gist.github.com/mages/5339689/raw/576263b8f0550125b61f4ddba127f5aa00fa2014/add.alpha.R")))

## Add an alpha value to a colour
add.alpha <- function(col=NULL, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
              rgb(x[1], x[2], x[3], alpha=alpha))
}

## load data
## expression matrix
dset <- loadDset ("GSE9893-breast")
gct <- dset$gct
str(gct)

## survival object
os <- dset$pheno$os
head(os)

## define the two batches from gpr files
gpr.dfr <- read.csv2 ("data/raw-data-reanalysis-GSE9893/GSE9893-gpr.csv")
colnames (gpr.dfr) <- sapply (gpr.dfr[1, ], gsub,
                              pattern = "^(.*)=(.*)",
                              replacement = "\\1")
colnames (gpr.dfr)[1] <- "gsm"
gpr.dfr$gsm <- gsub (".gpr", "", gpr.dfr$gsm)
## removes constant variables
gpr.dfr <- gpr.dfr[sapply (gpr.dfr, function (x) length (unique (x)) > 1)]
## to data table
gpr.dtb <- data.table (gpr.dfr, key = "gsm")
## reorders according to dset$pheno
gpr.dtb <- gpr.dtb[rownames (dset$pheno)]
## whenever appropriate, removes everything before '=' (included) from each column
gpr.dtb <- gpr.dtb[, lapply (.SD, gsub, pattern = "^(.*)=(.*)", replacement = "\\2")]
## transform variables to numeric
numVar <- c ("Temperature", "PMTGain", "ScanPower", "LaserPower")
for (j in numVar) set (gpr.dtb, j = j, value = as.numeric (gpr.dtb[[j]]))
## creates new variable defined by the presence of the string 'bordeaux' in the
## 'Settings' variable
gpr.dtb <- gpr.dtb[, Bordeaux := grepl ("bordeaux", Settings)]
## creates new variable by extracting the series number referenced in the
## variable 'ImageFiles'
gpr.dtb <- gpr.dtb[, Series := gsub (".*rie[[:space:]]*([0-9]+).*", "\\1", ImageFiles)]
gpr.dtb$Series <- as.numeric (gpr.dtb$Series)

gpr.dtb <- gpr.dtb[, Date := as.Date (DateTime, "%Y/%m/%d")]
gpr.dtb <- gpr.dtb[, Time := as.numeric (Date)]
gpr.dtb <- gpr.dtb[, Year := as.numeric (gsub ("^([0-9]+)-(.*)", "\\1", Date))]

os.dtb <- data.table(unclass(os))
os.dtb[, sampleName := rownames(dset[["pheno"]])]
setkey(os.dtb, "sampleName")

## merge gpr.dtb and os.dtb3
setnames(gpr.dtb, "gsm", "sampleName")
setkey(gpr.dtb, "sampleName")

dtb <- os.dtb[gpr.dtb][,
                       list(sampleName,
                            time,
                            status,
                            Date,
                            Year)]

## plot survival events 2005

ggplot(dtb,
       aes(x = time,
           colour = factor(status))) +
    geom_point(aes(y = -.02),
               position = position_jitter(height = .02),
               size = 2) +
    scale_colour_manual(values = c("0" = "black",
                                   "1" = "darkred")) +
    geom_density(data = dtb[status == 0],
                 aes(x = time)) +
    geom_density(data = dtb[status == 1],
                 aes(x = time),
                 colour = "darkred") +
    scale_y_continuous(limits = c(-.05, .05)) +
    ylab("") +
    xlab("time (months)") +
    facet_wrap(~ Year,
               ncol = 1) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

ggsave("~/Dropbox/thesis/graphics/gse9893-density-distribution-survival-events.pdf",
       width = 4,
       height = 3)

## chi-square test
chisq.test(table(dtb[, list(Year, status)]))

## pca
pcaRes <- prcomp(t(gct$data))
plot(pcaRes)

pt <- rep(1, nrow(gct$data))
pt[dset$pheno$status == "RF"] <- 3
col <- rep("black", nrow(gct$data))
col[gpr.dtb$Year == "2006"] <- "red"
colAlpha <- add.alpha(col, alpha = 0.5)

pdf("~/Dropbox/thesis/graphics/gse9893-pca-prior.pdf",
    width = 7,
    height = 7)

plot(pcaRes$x[, 1],
     pcaRes$x[, 2], col = colAlpha,
     pch = 19,
     xlab = paste("PC1(", round(100 * pcaRes$sdev[1] ^ 2 / sum(pcaRes$sdev ^ 2),
                                d = 1), "%)", sep = ""),
     ylab = paste("PC2(", round(100 * pcaRes$sdev[2] ^ 2 / sum(pcaRes$sdev ^ 2),
                                d = 1), "%)", sep = ""),
     cex = 2)

dev.off()

## renormalize data
## reads gpr files (these are one-color arrays scans)
gpr.dir <- "data/raw-data-reanalysis-GSE9893/GSE9893-gpr"
gpr.files <- file.path (gpr.dir, list.files (gpr.dir))

ff <- function (x) as.numeric(x$Flags > -99)
Cy5 <- "F635 Mean" # fools limma to read single clolor genepix file
RG <- read.maimages (gpr.files, source = "genepix.median",
                     wt.fun = ff,
                     columns = list (R = Cy5, G = Cy5),
                     verbose = FALSE)
RG$G <- NULL

## quantile-nomalize arrays
norm <- normalizeBetweenArrays (RG$R, method = "quantile")
norm <- log2 (norm)

## reorder columns
colnames (norm) <- gsub ("(.*)/(.*)$", "\\2", colnames (norm))
norm <- norm[, colnames (gct$data)]

## gene-wise average
norm <- aggregate (norm, by = list (RG$genes$Name), mean)
rownames (norm) <- norm[, 1]
norm <- as.matrix (norm[, -1])

## plot
pcaResNorm <- prcomp (t (norm))
plot (pcaResNorm)

pdf("~/Dropbox/thesis/graphics/gse9893-pca-post.pdf",
    width = 7,
    height = 7)

plot(pcaResNorm$x[, 1],
     pcaResNorm$x[,2],
     col = colAlpha,
     pch = 19,
     xlab = paste("PC1(", round(100 * pcaResNorm$sdev[1] ^ 2 / sum(pcaResNorm$sdev ^ 2),
                                d = 1), "%)", sep = ""),
     ylab = paste("PC2(", round(100 * pcaResNorm$sdev[2] ^ 2 / sum(pcaResNorm$sdev ^ 2),
                                d = 1), "%)", sep = ""),
     cex = 2)

dev.off()

## compute significant fractions post re-normalization
norm.gct <- list (row.descriptions = rownames (norm),
                  data = norm)

(frac.signif.pval (gct = norm.gct,
                   s = os))

prior <- testMSigDBgct(gct = dset$gct, surv = os, nCores = 1)
with(prior, sum(p.val < .05, na.rm = TRUE) / sum(!is.na(p.val)))

post <- testMSigDBgct(gct = norm.gct, surv = os, nCores = 1)
with(post, sum(p.val < .05, na.rm = TRUE) / sum(!is.na(p.val)))

postRnd <- testMSigDBgct(gct = norm.gct, surv = os, nCores = 1, is.rnd = TRUE)
with(postRnd, sum(p.val < .05, na.rm = TRUE) / sum(!is.na(p.val)))
