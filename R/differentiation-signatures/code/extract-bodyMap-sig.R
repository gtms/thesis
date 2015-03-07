## extractBodyMapSig
extractBodyMapSig <- function (tss, ## one of the 16 tissues of the
                               ## bodyMapRnk.mtx (see colnames)
                               top.tss = 1000,
                               top.other = 5000,
                               tol = 3,
                               prjDir = prj.dir) { ## must contain
    ## object bodyMapRnk.Rda in file.path (prj.dir, "raw-data")
    require (data.table)
    ## prj.dir <- "/Users/giltomas/projects/illumina"
    load (file.path (prjDir, "raw-data/bodyMapRnk.Rda"))

    top.keep <- bodyMapRnk.mtx[, tss] <= top.tss
    other.keep <- apply (bodyMapRnk.mtx[, colnames (bodyMapRnk.mtx) != tss], 1, function (gene) {
        sum (gene <= top.other) <= tol
    })

    rownames (bodyMapRnk.mtx)[top.keep & other.keep]
}
