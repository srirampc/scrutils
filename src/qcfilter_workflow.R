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
    cat(mcg_threshold)
    (mcgpct > mcg_threshold)
}

ngenes_cell_filter = function(dfx){
    per.cell = perCellQCMetrics(dfx, 
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    ngenes = per.cell$detected
    log_ngenes = log10(ngenes)
    ngenes_mad = mad(log_ngenes)
    ngenes_med = median(log_ngenes)
    ngenes_1mads = ngenes_med -  ngenes_mad
    cat(10^ngenes_1mads)
    (log_ngenes < ngenes_1mads)
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
    cat(10^lmean_1mads)
    (log_libmean < lmean_1mads)
}


apply_filters = function(out_dir, out_prefix, in_dirs){

    cat("DIR", "NGENES", "NCELLS", 
        "MCG_UB", "NGENES_LB", "AVGCTS_LB",
        "MCG_CELLS", "MCG_PCT",
        "NGENES_CELLS", "NGENES_PCT",
        "AVGCTS_CELS", "AVGCTS_PCT",
        "ALLF_CELLS", "ALLF_PCT", "\n")
    for(dirx in in_dirs){
        dfx = read10xCounts(dirx)
        cat(dirx, dim(dfx)[1], dim(dfx)[2])
        # per.feat = perFeatureQCMetrics(dfx)
        mcg_drop = mcg_cell_filter(dfx)
        ngenes_drop = ngenes_cell_filter(dfx)
        avgcts_drop = avgcounts_cell_filter(dfx)
        nlength = dim(dfx)[2]
        cat(sum(mcg_drop), sum(mcg_drop)*100/nlength, 
            sum(ngenes_drop), sum(ngenes_drop)*100/nlength, 
            sum(avgcts_drop), sum(avgcts_drop)*100/nlength,
            sum(mcg_drop | ngenes_drop | avgcts_drop),
            sum(mcg_drop | ngenes_drop | avgcts_drop)*100/nlength)
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
