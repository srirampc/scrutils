source("seurat_plots.R")
source("scater_plots.R")

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

    c(mcg=mcg_threshold, feature=feat_cells_threshold)
}

plot_qcstats = function(out_dir, out_prefix, in_dirs){
    for(dirx in in_dirs){
        dfx = read10xCounts(dirx)
        per.cell = perCellQCMetrics(dfx, 
            subset=list(MCG=grep("AT[MC]G",
                            rownames(dfx))))
        per.feat = perFeatureQCMetrics(dfx)
        #plot_flip_hist(per.cell, per.feat, dfx, dirx, out_dir, out_prefix)
        all_thrs = plot_hist(per.cell, per.feat, dfx, dirx, out_dir, out_prefix)
        #seurat_ba_qcplot(dirx, all_thrs, out_dir, out_prefix, ".png")
    }
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) >= 3){
    plot_qcstats(args[1], args[2], args[3:length(args)])
}  else {
    print(args)
    print("Usage: Rscript scater_qcplot.R outdir out_prefix indir1 indir2 ...")
}
