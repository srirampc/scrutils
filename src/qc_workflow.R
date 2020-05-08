library(argparser, quietly=TRUE)
source("qc_utils.R")
source("plot_utils.R")


apply_filter_dirs = function(out_dir, out_prefix, in_base_dir,
                             in_dirs, image_option,
                             inc_file, exc_file){

    cat("DIR", "NGENES", "NCELLS",
        "MCG_UB", "MCG_CELLS", "MCG_PCT",
        "NGENES_LB", "NGENES_LB_CELLS", "NGENES_LB_PCT",
        "NGENES_UB", "NGENES_UB_CELLS", "NGENES_UB_PCT",
        "AVGCTS_LB", "AVGCTS_CELS", "AVGCTS_PCT",
        # "MCGNG_CELLS", "MCGNG_PCT",
        # "ALLLB_CELLS", "ALLLB_PCT",
        "ALL_CELLS", "ALL_PCT",
        "FEAT_FILT", "FEAT_PCT")

    if(!is.null(inc_file) && !is.na(inc_file)) 
        cat(" ", "CODING_FILT", "CODING_PCT")

    if(!is.null(exc_file) && !is.na(exc_file)) 
        cat(" ", "PROTOPLAST_FILT", "PROTOPLAST_PCT")

    cat(" ", "NGENES_FINAL", "NCELLS_FINAL", "\n")

    for(dx in in_dirs){
        dirx = paste(in_base_dir, dx, sep="/")
        dfx = read10xCounts(dirx)
        cat(dirx, dim(dfx)[1], dim(dfx)[2])
        nlength = dim(dfx)[2]
        nfeatures = dim(dfx)[1]

        plot_cells_hist(dfx, dx, "-before-drop",
                        out_dir, out_prefix, image_option)

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
        plot_cells_hist(dfx4, dx, "-after-all-drop-", out_dir, out_prefix, image_option)

        cat(" ", sum(feat_drop), sum(feat_drop)*100/nfeatures)
        gnames = as.character(rownames(dfx4))
        #print(gnames[1:4])
        #print(gdf$name[1:4])

        if(!is.null(inc_file) && !is.na(inc_file)) {
          gdf = read.table(inc_file, header=TRUE, stringsAsFactors=FALSE)
          coding_drop = !(gnames %in% gdf$ID)
          dfx4 = dfx4[!coding_drop, ]
          cat(" ", sum(coding_drop), sum(coding_drop)*100/nfeatures)
	}

        if(!is.null(exc_file) && !is.na(exc_file)) {
          pdf = read.table(exc_file, header=TRUE, stringsAsFactors=FALSE)
          protoplast_drop = (gnames %in% pdf$ID)
          dfx4 = dfx4[!protoplast_drop, ]
          cat(" ", sum(protoplast_drop), sum(protoplast_drop)*100/nfeatures)
	}
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
        filter_venn(filter_list, filter_names, venn_fname, image_option)

        cat("\n")
    }
}

qcwf_main = function(in_base_dir, data.file,
                    out_dir, out_prefix, image_option,
                    include_genes_file, exclude_genes_file){
    data.df = read.csv(data.file, header=TRUE, stringsAsFactors=FALSE)
    expt_dir_paths = data.df$dir.paths
    apply_filter_dirs(out_dir, out_prefix, in_base_dir,
                      expt_dir_paths, image_option,
                      include_genes_file, exclude_genes_file)
    
}

# Create a parser
p <- arg_parser("Run quality control workflow to generate qc plots and  qc stats")

# Add command line arguments
p <- add_argument(p, "root_dir", help="Root directory of datasets", type="character")
p <- add_argument(p, "data_file_csv", 
                  help="CSV file with dataset info (See ath.control.csv for example)",
                  type="character")
p <- add_argument(p, "out_dir", help="Output directory", type="character")
p <- add_argument(p, "out_prefix", help="Output Prefix", type="character")
p <- add_argument(p, "--img", 
                  help="Output image option should be one png/pdf (default:png)", short='-g', default='png')
p <- add_argument(p, "--exc", help="File containing list of genes to be included eg. protien coding (default:None)", short='-e', default=NULL)
p <- add_argument(p, "--inc", help="File containing list of genes to be excluded eg. protoplast (default:None)", short='-i', default=NULL)

# Parse the command line arguments
argv <- parse_args(p)

if(argv$img == "png" || argv$img == "pdf") {
	print(argv$inc)
    qcwf_main(argv$root_dir, argv$data_file_csv,
        argv$out_dir, argv$out_prefix,  
        argv$img, argv$inc, argv$exc)    
} else {
    print("Invalid image option.")
    print(p)
}

