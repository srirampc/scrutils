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


apply_filters = function(out_dir, out_prefix, in_dirs){

    cat("DIR", "NGENES", "NCELLS", 
        "MCG_UB", "MCG_CELLS", "MCG_PCT",
        "NGENES_LB", "NGENES_LB_CELLS", "NGENES_LB_PCT",
        "NGENES_UB", "NGENES_UB_CELLS", "NGENES_UB_PCT",
        "AVGCTS_LB", "AVGCTS_CELS", "AVGCTS_PCT",
        "MCGNG_CELLS", "MCGNG_PCT", 
        "ALLLB_CELLS", "ALLLB_PCT", 
        "ALL_CELLS", "ALL_PCT", 
        "\n")
    for(dirx in in_dirs){
        dfx = read10xCounts(dirx)
        cat(dirx, dim(dfx)[1], dim(dfx)[2])
        nlength = dim(dfx)[2]
        # per.feat = perFeatureQCMetrics(dfx)
        mcg_drop = mcg_cell_filter(dfx)
        cat(" ", sum(mcg_drop), sum(mcg_drop)*100/nlength)
        ngenes_lb_drop = ngenes_cell_filter_lb(dfx)
        cat(" ", sum(ngenes_lb_drop), sum(ngenes_lb_drop)*100/nlength)
        ngenes_ub_drop = ngenes_cell_filter_ub(dfx)
        cat(" ", sum(ngenes_lb_drop), sum(ngenes_lb_drop)*100/nlength)
        avgcts_drop = avgcounts_cell_filter(dfx)
        cat(" ", sum(avgcts_drop), sum(avgcts_drop)*100/nlength)
        cat(" ",
            sum(mcg_drop | ngenes_lb_drop | ngenes_ub_drop),
            sum(mcg_drop | ngenes_lb_drop | ngenes_ub_drop)*100/nlength)
        cat(" ",
            sum(mcg_drop | ngenes_lb_drop | avgcts_drop),
            sum(mcg_drop | ngenes_lb_drop | avgcts_drop)*100/nlength)
        cat(" ",
            sum(mcg_drop | ngenes_lb_drop | ngenes_ub_drop | avgcts_drop),
            sum(mcg_drop | ngenes_lb_drop | ngenes_ub_drop | avgcts_drop)*100/nlength)
       cat("\n")
    }
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) >= 3){
    apply_filters(args[1], args[2], args[3:length(args)])
}  else {
    print(args)
    print("Usage: Rscript qcfilter_workflow.R outdir out_prefix indir1 indir2 ...")
}
