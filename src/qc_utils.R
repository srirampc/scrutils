library(scater, quietly=T)
library(DropletUtils, quietly=T)
library(Seurat, quietly=T)
library(ggplot2, quietly=T)

mcg_cell_filter = function(dfx, print=T){
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
    if (print) cat(" ", mcg_threshold)
    (mcgpct > mcg_threshold)
}

ngenes_cell_filter_lb = function(dfx, print=T){
    per.cell = perCellQCMetrics(dfx, 
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    ngenes = per.cell$detected
    log_ngenes = log10(ngenes)
    ngenes_mad = mad(log_ngenes)
    ngenes_med = median(log_ngenes)
    ngenes_1mads = ngenes_med -  ngenes_mad
    if(print) cat(" ", 10^ngenes_1mads)
    (log_ngenes < ngenes_1mads)
}

ngenes_cell_filter_ub = function(dfx, print=T){
    per.cell = perCellQCMetrics(dfx, 
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    ngenes = per.cell$detected
    log_ngenes = log10(ngenes)
    ngenes_mad = mad(log_ngenes)
    ngenes_med = median(log_ngenes)
    ngenes_3madsp = ngenes_med + (3*ngenes_mad)
    if(print) cat(" ", 10^ngenes_3madsp)
    (log_ngenes > ngenes_3madsp)
}

logcounts_cell_filter = function(dfx, print=T){
    per.cell = perCellQCMetrics(dfx, 
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    libsize = per.cell$sum 
    log_libsize = log10(libsize)
    lsize_mad = mad(log_libsize)
    lsize_med = median(log_libsize)
    lsize_1mads = lsize_med - lsize_mad
    if(print) cat(" ", 10^lsize_1mads)
    (log_libsize < lsize_1mads)
}

avgcounts_cell_filter = function(dfx, print=T){
    per.cell = perCellQCMetrics(dfx, 
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    libmean = per.cell$sum / dim(dfx)[2]
    log_libmean = log10(libmean)
    lmean_mad = mad(log_libmean)
    lmean_med = median(log_libmean)
    lmean_1mads = lmean_med - lmean_mad
    if (print) cat(" ", 10^lmean_1mads)
    (log_libmean < lmean_1mads)
}

avg_reads_feat_filter = function(dfx, print=T){
    per.feat = perFeatureQCMetrics(dfx)
    featavg = per.feat$mean 
    log_featavg = log10(featavg)

    xhist = hist(log_featavg, breaks=120, plot=FALSE)

    xrange = (xhist$mids < -1) & (xhist$mids > -4)  & (xhist$density > 0.05)
    xrngcts =xhist$counts[xrange]
    xrngmids = xhist$mids[xrange]
    xrnmidv = xrngmids[which.min(xrngcts)]
    if (print) cat(" ", 10^xrnmidv)
    (log_featavg <= xrnmidv)
}


apply_cell_filters = function(dfx, print=T){
    mcg_drop = mcg_cell_filter(dfx, print)
    ngenes_lb_drop = ngenes_cell_filter_lb(dfx, print)
    ngenes_ub_drop = ngenes_cell_filter_ub(dfx, print)
    #avgcts_drop = avgcounts_cell_filter(dfx, print)
    logcts_drop = logcounts_cell_filter(dfx, print)
    all_drop = mcg_drop | ngenes_lb_drop | ngenes_ub_drop | logcts_drop
    dfx[,!all_drop]
}

apply_gene_filters = function(dfx, print=T){
    feat_drop = avg_reads_feat_filter(dfx, print)
    dfx[!feat_drop, ]
}

qc_normalize = function(expt.full.dir.path, 
                        include_genes_file,
                        exclude_genes_file){
    dfx = read10xCounts(expt.full.dir.path)
    dfx = apply_cell_filters(dfx, F)
    dfx = apply_gene_filters(dfx, F)
    #print(dfx)
    if(!is.null(include_genes_file) && !is.na(include_genes_file)) {
        gnames = as.character(rownames(dfx))
        inc_df = read.table(include_genes_file, 
            header=TRUE, stringsAsFactors=FALSE)
        inc_flag = (gnames %in% inc_df[,'ID'])
        dfx = dfx[inc_flag, ]
    }
    if(!is.null(exclude_genes_file) && !is.na(exclude_genes_file)) {
        gnames = as.character(rownames(dfx))
        exc_df = read.table(exclude_genes_file, header=TRUE, stringsAsFactors=FALSE)
        exc_drop = !(gnames %in% exc_df[,'ID'])
        dfx = dfx[exc_drop, ]
    }
    dfx = scran_normalize(dfx)
    dfx
}

qc_normalize_matrix = function(expt.full.dir.path,
                               include_genes_file,
                               exclude_genes_file){
    dfx = qc_normalize(expt.full.dir.path,
                       include_genes_file,
                       exclude_genes_file)
    ctmtx = counts(dfx)
    print(dim(ctmtx))
    rownames(ctmtx) = rownames(dfx)
    colnames(ctmtx) = 1:dim(dfx)[2]
    ctmtx
}
