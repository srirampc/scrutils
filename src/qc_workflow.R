source("qc_utils.R")
source("plot_utils.R")


apply_filter_dirs = function(out_dir, out_prefix, glist_file,
                             in_base_dir, in_dirs, image.option){

    cat("DIR", "NGENES", "NCELLS",
        "MCG_UB", "MCG_CELLS", "MCG_PCT",
        "NGENES_LB", "NGENES_LB_CELLS", "NGENES_LB_PCT",
        "NGENES_UB", "NGENES_UB_CELLS", "NGENES_UB_PCT",
        "AVGCTS_LB", "AVGCTS_CELS", "AVGCTS_PCT",
        # "MCGNG_CELLS", "MCGNG_PCT",
        # "ALLLB_CELLS", "ALLLB_PCT",
        "ALL_CELLS", "ALL_PCT",
        "FEAT_FILT", "FEAT_PCT",
        "CODING_FILT", "CODING_PCT",
        "NGENES_FINAL", "NCELLS_FINAL",
        "\n")
    gdf = read.table(glist_file, header=TRUE, stringsAsFactors=FALSE)
    for(dx in in_dirs){
        dirx = paste(in_base_dir, dx, sep="/")
        dfx = read10xCounts(dirx)
        cat(dirx, dim(dfx)[1], dim(dfx)[2])
        nlength = dim(dfx)[2]
        nfeatures = dim(dfx)[1]

        plot_cells_hist(dfx, dx, "-before-drop",
                        out_dir, out_prefix, image.option)

        mcg_drop = mcg_cell_filter(dfx)
        ngenes_lb_drop = ngenes_cell_filter_lb(dfx)
        ngenes_ub_drop = ngenes_cell_filter_ub(dfx)
        #avgcts_drop = avgcounts_cell_filter(dfx)
        logcts_drop = logcounts_cell_filter(dfx)


        #mcg_ngenes_drop = mcg_drop | ngenes_lb_drop | ngenes_ub_drop
        #all_lb_drop = mcg_drop | ngenes_lb_drop | avgcts_drop
        all_drop = mcg_drop | ngenes_lb_drop | ngenes_ub_drop | logcts_drop
        
        cat(" ", sum(mcg_drop), sum(mcg_drop)*100/nlength)
        cat(" ", sum(ngenes_lb_drop), sum(ngenes_lb_drop)*100/nlength)
        cat(" ", sum(ngenes_ub_drop), sum(ngenes_ub_drop)*100/nlength)
        #cat(" ", sum(avgcts_drop), sum(avgcts_drop)*100/nlength)
        cat(" ", sum(logcts_drop), sum(logcts_drop)*100/nlength)

        # cat(" ", sum(mcg_ngenes_drop), sum(mcg_ngenes_drop)*100/nlength)
        # dfx2 = dfx[, !mcg_ngenes_drop]
        # plot_cells_hist(dfx2, dx, "-after-mcg-ng-drop", out_dir, out_prefix)
        # cat(" ", sum(all_lb_drop), sum(all_lb_drop)*100/nlength)
        # dfx3 = dfx[, !all_lb_drop]
        # plot_cells_hist(dfx3, dx, "-after-all-lb-drop",  out_dir, out_prefix)

        cat(" ", sum(all_drop), sum(all_drop)*100/nlength)
        dfx4 = dfx[, !all_drop]
        #print(dim(dfx4))

        feat_drop = avg_reads_feat_filter(dfx)
        dfx4 = dfx4[!feat_drop, ]
        plot_cells_hist(dfx4, dx, "-after-all-drop-", out_dir, out_prefix, image.option)

        cat(" ", sum(feat_drop), sum(feat_drop)*100/nfeatures)
        gnames = as.character(rownames(dfx4))
        #print(gnames[1:4])
        #print(gdf$name[1:4])
        coding_drop = !(gnames %in% gdf$name)
        dfx4 = dfx4[!coding_drop, ]
        cat(" ", sum(coding_drop), sum(coding_drop)*100/nfeatures)
        cat(" ", dim(dfx4)[1], dim(dfx4)[2])
        ncell_list = 1:nlength
        filter_list = list(
           ncell_list[mcg_drop], 
           ncell_list[ngenes_ub_drop|ngenes_lb_drop],
           ncell_list[logcts_drop])
        filter_names = c("MCG filter", "Genes filter", "Counts filter")
        venn_fname = paste(out_dir, dx,
           paste(out_prefix, "-filter-venn.png", sep=""),
        sep="/"  )
        filter_venn(filter_list, filter_names, venn_fname image.option)

        cat("\n")
    }
}

qcwf_main(in_base_dir, data.file, glist_file, 
          out_dir, out_prefix, image.option){
    data.df = read.csv(data.file, header=TRUE, stringsAsFactors=FALSE)
    expt.dir.paths = data.df$dir.paths
    short.names = data.df$short.names
    project.names = data.df$project.names

    apply_filter_dirs(out_dir, out_prefix, glist_file,
                      in_base_dir, expt.dir.paths, image.option)
    
}

args = commandArgs(trailingOnly=TRUE)
cmd_usage = "Usage: Rscript qcfilter_workflow.R in_base_dir data_file gene_list outdir out_prefix png/pdf"
if(length(args) >= 6){
    if((args[6] == "png" || args[6] == "pdf")){
        qcwf_main(args[1], args[2], args[3], args[4], args[5], args[6])
    } else {
        print(args)
        print(cmd_usage)
    }
}  else {
    print(args)
    print(cmd_usage)
}
