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

    xhist = hist(log_featavg, breaks=120, plot=FALSE)

    xrange = (xhist$mids < -1) & (xhist$mids > -4)  & (xhist$density > 0.05)
    xrngcts =xhist$counts[xrange]
    xrngmids = xhist$mids[xrange]
    xrnmidv = xrngmids[which.min(xrngcts)]
    (log_featavg <= xrnmidv)
}


apply_cell_filters = function(dfx){
    mcg_drop = mcg_cell_filter(dfx)
    ngenes_lb_drop = ngenes_cell_filter_lb(dfx)
    ngenes_ub_drop = ngenes_cell_filter_ub(dfx)
    #avgcts_drop = avgcounts_cell_filter(dfx)
    logcts_drop = logcounts_cell_filter(dfx)
    all_drop = mcg_drop | ngenes_lb_drop | ngenes_ub_drop | logcts_drop
    dfx[,!all_drop]
}

apply_gene_filters = function(dfx){
    feat_drop = avg_reads_feat_filter(dfx)
    dfx[!feat_drop, ]
}
