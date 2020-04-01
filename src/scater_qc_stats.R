library("DropletUtils")
library(scater)

qcstat_dir = function(dirx){
    dfx = read10xCounts(dirx)
    per.cell = perCellQCMetrics(dfx, 
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    per.feat = perFeatureQCMetrics(dfx)
    tfx = per.feat$detected * dim(dfx)[2]/100.0
    libsize_mads = isOutlier(per.cell$sum, nmads=3,
                                type="lower", log=FALSE)
    genects_mads = isOutlier(per.cell$detected, nmads=3,
                                type="lower", log=FALSE)
    feat_cell_cts = per.feat$detected * dim(dfx)[2]/100.0
    feat_cell_mean = per.feat$mean >= 1.0

    libsize_95p = per.cell$sum < (quantile(per.cell$sum, .95) * 0.05)
    genects_99p = per.cell$detected < (quantile(per.cell$detected, 0.99) * 0.05)

    mcg_nmads = isOutlier(per.cell$subsets_MCG_percent, nmads=3, typ="higher")
    mcg_5p = per.cell$subsets_MCG_percent > 5
    mcg_4p = per.cell$subsets_MCG_percent > 4
    mcg_3p = per.cell$subsets_MCG_percent > 3
    c(#dirname = dirx,
      total_cells = dim(dfx)[2],
      libsize_mads=sum(libsize_mads),
      genects_mads=sum(genects_mads),
      libsize_95p=sum(libsize_95p),
      genects_99p=sum(genects_99p),
      mcg_nmads=sum(mcg_nmads),
      mcg_5p=sum(mcg_5p),
      mcg_4p=sum(mcg_4p),
      mcg_3p=sum(mcg_3p),
      total_features = dim(dfx)[1], 
      feat_cell_mean=sum(feat_cell_mean)
       )

}

get_qcstats = function(outfile, indirs){
    print(outfile)
    print(indirs)
    qcl = lapply(indirs, qcstat_dir)
    qdf = t(as.data.frame(qcl))
    rownames(qdf) = indirs
    write.table(qdf, outfile, sep="\t", quote=FALSE)
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) >= 2){
    get_qcstats(args[1], args[2:length(args)])
}  else {
    print(args)
    print("Usage: Rscript scater_qc_stats.R out_file indir1 indir2 ...")
}
