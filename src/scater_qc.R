source("seurat_plots.R")
source("scater_plots.R")

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

plot_flip_hist = function(per.cell, per.feat, dfx, dirx, out_dir, out_prefix){
    fname = paste(out_dir, dirx, paste(out_prefix, "-scater-qc-feat-filp-hist.png", sep=""), sep="/"  )
    plot_feat_counts_flip_hist(per.feat, dfx, dirx, fname)
    
    fname = paste(out_dir, dirx, paste(out_prefix, "-scater-qc-counts-flip-hist.png", sep=""), sep="/"  )
    plot_counts_flip_hist(per.cell, dirx, fname)
    
    fname = paste(out_dir, dirx, paste(out_prefix, "-scater-qc-ngenes-flip-hist.png", sep=""), sep="/"  )
    plot_cell_ngenes_flip_hist(per.cell, dirx, fname)

    fname = paste(out_dir, dirx, paste(out_prefix, "-scater-qc-mcgpct-flip-hist.png", sep=""), sep="/"  )
    plot_mcg_flip_hist(per.cell, dirx, fname)

    fname = paste(out_dir, dirx, paste(out_prefix, "-scater-qc-favg-flip-hist.png", sep=""), sep="/"  )
    plot_feat_avg_flip_hist(per.feat, dirx, fname)
}

plot_hist = function(per.cell, per.feat, dfx, dirx, out_dir, out_prefix){
    cat("Dims: ",  dirx, dim(dfx), "\n")

    fname = paste(out_dir, dirx, 
        paste(out_prefix,
            "-scater-qc-mcgpct-hist.png", sep=""),
        sep="/"  )
    mcg_threshold = plot_mcg_hist(per.cell, dirx, fname)

    fname = paste(out_dir, dirx, 
        paste(out_prefix, 
            "-scater-qc-feat-hist.png", sep=""),
        sep="/"  )
    feat_cells_threshold = plot_feat_counts_hist(per.feat, dfx, dirx, fname)

    fname = paste(out_dir, dirx,
        paste(out_prefix, "-scater-qc-feats-avg-hist.png", sep=""),
        sep="/"  )
    plot_feat_avg_hist(per.feat, dirx, fname)

    fname = paste(out_dir, dirx,
        paste(out_prefix, "-scater-qc-counts-avg-hist.png", sep=""),
        sep="/"  )
    plot_counts_avg_hist(per.cell, dfx, dirx, fname)

    fname = paste(out_dir, dirx, 
        paste(out_prefix, "-scater-qc-counts-logavg-hist.png", sep=""),
        sep="/"  )
    lmean_3madsp = plot_counts_logavg_hist(per.cell, dfx, dirx, fname)

    fname = paste(out_dir, dirx,
        paste(out_prefix, "-scater-qc-counts-logavg-hist2.png", sep=""),
        sep="/"  )
    plot_counts_logavg_hist2(per.cell, lmean_3madsp, dfx, dirx, fname)

    fname = paste(out_dir, dirx,
        paste(out_prefix, "-scater-qc-ngenes-hist.png", sep=""),
        sep="/"  )
    plot_ngenes(per.cell, dirx, fname)

    c(mcg=mcg_threshold, feature=feat_cells_threshold, 200)
}

plot_qcstats = function(out_dir, out_prefix, in_base_dir, in_dirs){
    outfile = paste(out_dir, "scater-qc.tsv", sep="/")
    get_qcstats(outfile, paste(in_base_dir, in_dirs, sep="/"))
    for(dirx in in_dirs){
        cx.indir = paste(in_base_dir, dirx, sep="/")
        dfx = read10xCounts(cx.indir)
        per.cell = perCellQCMetrics(dfx, 
            subset=list(MCG=grep("AT[MC]G",
                            rownames(dfx))))
        per.feat = perFeatureQCMetrics(dfx)
        #plot_flip_hist(per.cell, per.feat, dfx, dirx, out_dir, out_prefix)
        all_thrs = plot_hist(per.cell, per.feat, dfx, dirx, out_dir, out_prefix)
        seurat_ba_qcplot(cx.indir, all_thrs, dirx, out_dir, out_prefix, ".png")
    }
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) >= 4){
    plot_qcstats(args[1], args[2], args[3], args[4:length(args)])
}  else {
    print(args)
    print("Usage: Rscript scater_qcplot.R outdir out_prefix in_base_dir indir1 indir2 ...")
}
