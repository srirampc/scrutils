source("scater_qc.R")
source("scater_plots.R")

apply_filters = function(out_dir, out_prefix, in_base_dir, in_dirs){

    cat("DIR", "NGENES", "NCELLS",
        "MCG_UB", "MCG_CELLS", "MCG_PCT",
        "NGENES_LB", "NGENES_LB_CELLS", "NGENES_LB_PCT",
        "NGENES_UB", "NGENES_UB_CELLS", "NGENES_UB_PCT",
        "AVGCTS_LB", "AVGCTS_CELS", "AVGCTS_PCT",
        # "MCGNG_CELLS", "MCGNG_PCT",
        # "ALLLB_CELLS", "ALLLB_PCT",
        "ALL_CELLS", "ALL_PCT",
        "FEAT_FILT", "FEAT_PCT",
        "NGENES_FINAL", "NCELLS_FINAL",
        "\n")
    for(dx in in_dirs){
        dirx = paste(in_base_dir, dx, sep="/")
        dfx = read10xCounts(dirx)
        cat(dirx, dim(dfx)[1], dim(dfx)[2])
        nlength = dim(dfx)[2]
        nfeatures = dim(dfx)[1]
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
        # mcg_ngenes_drop = mcg_drop | ngenes_lb_drop | ngenes_ub_drop
        # all_lb_drop = mcg_drop | ngenes_lb_drop | logcts_drop
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
        # dfx2 = dfx[, !mcg_ngenes_drop]
        # plot_cells_hist(dfx2, dx, "-after-mcg-ng-drop",
        #                 out_dir, out_prefix)
        # dfx3 = dfx[, !all_lb_drop]
        # plot_cells_hist(dfx3, dx, "-after-all-lb-drop",
        #                 out_dir, out_prefix)
        dfx4 = dfx[, !all_drop]
        #print(dim(dfx4))
        feat_drop = avg_reads_feat_filter(dfx)
        dfx4 = dfx4[!feat_drop, ]
        plot_cells_hist(dfx4, dx, "-after-all-drop-",
                        out_dir, out_prefix)

        cat(" ", sum(feat_drop),sum(feat_drop)*100/nfeatures,
            dim(dfx4)[1], dim(dfx4)[2])
        ncell_list = 1:nlength
        filter_list = list(ncell_list[mcg_drop], 
                ncell_list[ngenes_ub_drop|ngenes_lb_drop],
                ncell_list[logcts_drop])
        filter_names = c("MCG", "Genes filter", "Counts filter")
        venn_fname = paste(out_dir, dx,
           paste(out_prefix, "-filter-venn.png", sep=""),
        sep="/"  )
        filter_venn(filter_list, filter_names, venn_fname)

        cat("\n")
    }
}



args = commandArgs(trailingOnly=TRUE)

if(length(args) >= 4){
    apply_filters(args[1], args[2], args[3], args[4:length(args)])
}  else {
    print(args)
    print("Usage: Rscript qcfilter_workflow.R outdir out_prefix in_base_dir indir1 indir2 ...")
}
