library(scater)
library(DropletUtils)
library(Seurat)
library(ggplot2)

mcg_cell_filter = function(dfx){
    per.cell = perCellQCMetrics(dfx, 
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    mcgpct = per.cell$subsets_MCG_percent
    log_mcg = log10(mcgpct)
    log_mcgf = log_mcg
    log_mcg = log_mcg[log_mcg > -2.3]
    mcg_mad = mad(log_mcg)
    mcg_med = median(log_mcg)
    mcg_3madsp = mcg_med + 3 * mcg_mad
    mcg_threshold = min(5, 10^mcg_3madsp)
    cat(" ", mcg_threshold)
    (mcgpct > mcg_threshold)
}

ngenes_cell_filter_lb = function(dfx){
    per.cell = perCellQCMetrics(dfx, 
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    ngenes = per.cell$detected
    log_ngenes = log10(ngenes)
    ngenes_mad = mad(log_ngenes)
    ngenes_med = median(log_ngenes)
    ngenes_1mads = ngenes_med -  ngenes_mad
    cat(" ", 10^ngenes_1mads)
    (log_ngenes < ngenes_1mads)
}

ngenes_cell_filter_ub = function(dfx){
    per.cell = perCellQCMetrics(dfx, 
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    ngenes = per.cell$detected
    log_ngenes = log10(ngenes)
    ngenes_mad = mad(log_ngenes)
    ngenes_med = median(log_ngenes)
    ngenes_3madsp = ngenes_med + (3*ngenes_mad)
    cat(" ", 10^ngenes_3madsp)
    (log_ngenes > ngenes_3madsp)
}

logcounts_cell_filter = function(dfx){
    per.cell = perCellQCMetrics(dfx, 
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    libsize = per.cell$sum 
    log_libsize = log10(libsize)
    lsize_mad = mad(log_libsize)
    lsize_med = median(log_libsize)
    lsize_1mads = lsize_med - lsize_mad
    cat(" ", 10^lsize_1mads)
    (log_libsize < lsize_1mads)
}



avgcounts_cell_filter = function(dfx){
    per.cell = perCellQCMetrics(dfx, 
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    libmean = per.cell$sum / dim(dfx)[2]
    log_libmean = log10(libmean)
    lmean_mad = mad(log_libmean)
    lmean_med = median(log_libmean)
    lmean_1mads = lmean_med - lmean_mad
    cat(" ", 10^lmean_1mads)
    (log_libmean < lmean_1mads)
}

avg_reads_feat_filter = function(dfx){
    per.feat = perFeatureQCMetrics(dfx)
    featavg = per.feat$mean 
    log_featavg = log10(featavg)

    xhist = hist(log_featavg, breaks=120, xlab=expression(Log[10]~"Avg. No. Reads"), ylab="No. cells", 
        main=dirx, prob=FALSE, plot=FALSE)

    xrange = (xhist$mids < -1) & (xhist$mids > -4)  & (xhist$density > 0.05)
    xrngcts =xhist$counts[xrange]
    xrngmids = xhist$mids[xrange]
    xrnmidv = xrngmids[which.min(xrngcts)]
    (log_featavg > xrnmidv)
}


plot_cells_hist = function(dfx, dirx, name_prefix,
                           out_dir, out_prefix){
    per.cell = perCellQCMetrics(dfx,
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    fname = paste(out_dir, dirx,
        paste(out_prefix, name_prefix,
              "-scater-qc-logcounts-hist.png", sep=""),
        sep="/"  )
    # print(fname)
    plot_logcounts_hist(per.cell, dfx, dirx, fname, FALSE)

    fname = paste(out_dir, dirx,
        paste(out_prefix, name_prefix,
              "-scater-qc-counts-logavg-hist.png", sep=""),
        sep="/"  )
    # print(fname)
    plot_counts_logavg_hist(per.cell, dfx, dirx, fname, FALSE)

    fname = paste(out_dir, dirx,
        paste(out_prefix, name_prefix,
          "-scater-qc-ngenes-hist.png", sep=""),
        sep="/"  )
    plot_ngenes(per.cell, dirx, fname, FALSE)
}

apply_filters = function(out_dir, out_prefix, in_base_dir, in_dirs){

    cat("DIR", "NGENES", "NCELLS",
        "MCG_UB", "MCG_CELLS", "MCG_PCT",
        "NGENES_LB", "NGENES_LB_CELLS", "NGENES_LB_PCT",
        "NGENES_UB", "NGENES_UB_CELLS", "NGENES_UB_PCT",
        "AVGCTS_LB", "AVGCTS_CELS", "AVGCTS_PCT",
        # "MCGNG_CELLS", "MCGNG_PCT",
        # "ALLLB_CELLS", "ALLLB_PCT",
        "ALL_CELLS", "ALL_PCT",
        "\n")
    for(dx in in_dirs){
        dirx = paste(in_base_dir, dx, sep="/")
        dfx = read10xCounts(dirx)
        cat(dirx, dim(dfx)[1], dim(dfx)[2])
        nlength = dim(dfx)[2]
        # per.feat = perFeatureQCMetrics(dfx)

        plot_cells_hist(dfx, dx, "-before-drop",
                        out_dir, out_prefix)
        mcg_drop = mcg_cell_filter(dfx)
        cat(" ", sum(mcg_drop), sum(mcg_drop)*100/nlength)

        ngenes_lb_drop = ngenes_cell_filter_lb(dfx)
        cat(" ", sum(ngenes_lb_drop), sum(ngenes_lb_drop)*100/nlength)
        ngenes_ub_drop = ngenes_cell_filter_ub(dfx)
        cat(" ", sum(ngenes_ub_drop), sum(ngenes_ub_drop)*100/nlength)

        #avgcts_drop = avgcounts_cell_filter(dfx)
        #cat(" ", sum(avgcts_drop), sum(avgcts_drop)*100/nlength)
        #mcg_ngenes_drop = mcg_drop | ngenes_lb_drop | ngenes_ub_drop
        #all_lb_drop = mcg_drop | ngenes_lb_drop | avgcts_drop
        #all_drop = mcg_drop | ngenes_lb_drop | ngenes_ub_drop | avgcts_drop

        logcts_drop = logcounts_cell_filter(dfx)
        cat(" ", sum(logcts_drop), sum(logcts_drop)*100/nlength)
        mcg_ngenes_drop = mcg_drop | ngenes_lb_drop | ngenes_ub_drop
        all_lb_drop = mcg_drop | ngenes_lb_drop | logcts_drop
        all_drop = mcg_drop | ngenes_lb_drop | ngenes_ub_drop | logcts_drop

        # cat(" ",
        #     sum(mcg_ngenes_drop),
        #     sum(mcg_ngenes_drop)*100/nlength)
        # cat(" ",
        #     sum(all_lb_drop),
        #     sum(all_lb_drop)*100/nlength)
        cat(" ",
            sum(all_drop),
            sum(all_drop)*100/nlength)
        dfx2 = dfx[, !mcg_ngenes_drop]
        plot_cells_hist(dfx2, dx, "-after-mcg-ng-drop",
                        out_dir, out_prefix)
        dfx3 = dfx[, !all_lb_drop]
        plot_cells_hist(dfx3, dx, "-after-all-lb-drop",
                        out_dir, out_prefix)
        dfx4 = dfx[, !all_drop]
        #print(dim(dfx4))
        dfx4 = dfx4[feat_drip, ]
        cat(" ", sum(feat_drip),sum(feat_drip)*100/nfeatures,
            dim(dfx4)[1], dim(dfx4)[2])
        plot_cells_hist(dfx4, dx, "-after-all-drop-",
                        out_dir, out_prefix)
        cat("\n")
    }
}

