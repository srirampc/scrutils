
library(scran)
library(DropletUtils)
source("qc_utils.R")
source("plot_utils.R")
source("data_utils.R")

variance_plot = function(dfx, fname){
    dec = scran::modelGeneVar(dfx)
    if(endsWith(fname, "png")) png(file=fname) else pdf(file=fname)
    plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
    curve(metadata(dec)$trend(x), col="blue", add=TRUE)
    dev.off()
}

scran_normalize_dir = function(out_dir, out_prefix, 
                            in_base_dir, dirx, 
                            plot=FALSE, image.option="png"){
    rdx = paste(in_base_dir, dirx, sep="/")
    dfx = DropletUtils::read10xCounts(rdx)
    cat(dirx)
    if(plot == TRUE){
        dfx2 = scater::logNormCounts(dfx)
        #print(dim(dfx2))
        fname = paste(out_dir, dirx, 
            paste(out_prefix, "-variance-before.", image.option, sep=""),
            sep="/"  )
        variance_plot(dfx2, fname)
        fname = paste(out_dir, dirx,
            paste(out_prefix, "-seurat-scatter-before.", image.option, sep=""),
        sep="/"  )
        ctmtx = counts(dfx2)
        #print(dim(ctmtx))
        rownames(ctmtx) = rownames(dfx2)
        colnames(ctmtx) = 1:dim(dfx2)[2]
        scrj = CreateSeuratObject(counts = ctmtx, project = dirx)
        seurat_fscatter(scrj, fname)
    }
    dfx = apply_cell_filters(dfx)
    dfx = apply_gene_filters(dfx)
    #print(dfx)
    dfx = scran_normalize(dfx)

    if(plot == TRUE){
        fname = paste(out_dir, dirx,
            paste(out_prefix, "-variance-after.", image.option, sep=""),
        sep="/"  )
        variance_plot(dfx, fname)
        fname = paste(out_dir, dirx,
            paste(out_prefix, "-seurat-scatter-after.", image.option, sep=""),
        sep="/"  )
        ctmtx = counts(dfx)
        print(dim(ctmtx))
        rownames(ctmtx) = rownames(dfx)
        colnames(ctmtx) = 1:dim(dfx)[2]
        scrj = CreateSeuratObject(counts = ctmtx, project = dirx)
        seurat_fscatter(scrj, fname)
    }
    cat("\n")
    dfx
}

scran_main = function(in_base_dir, data.file, 
                      out_dir, out_prefix, plot, image.option){
    data.df = read.csv(data.file, header=TRUE, stringsAsFactors=FALSE)
    expt.dir.paths = data.df$dir.paths
    short.names = data.df$short.names
    project.names = data.df$project.names
    plot = as.logical(plot)

    for(dirx in expt.dir.paths){
        scran_normalize_dir(out_dir, out_prefix, 
                            in_base_dir, dirx, plot, image.option)
    }
}

# Create a parser
p <- arg_parser("Normalize with scran and generate variance plots")

# Add command line arguments
p <- add_argument(p, "root_dir", help="Root directory of datasets", type="character")
p <- add_argument(p, "data_file_csv", 
                  help="CSV file with dataset info (See ath.control.csv for example)",
                  type="character")
p <- add_argument(p, "out_dir", help="Output directory", type="character")
p <- add_argument(p, "out_prefix", help="Output Prefix", type="character")
p <- add_argument(p, "--plot", help="Flag if generate plots or not TRUE/FALSE (default: TRUE)", short='-p', default=TRUE)
p <- add_argument(p, "--img", help="Output image option should be one png/pdf (default:png)", short='-g', default='png')

if(!(argv$img == "png" || argv$img == "pdf")) {
    scran_main(argv$in_base_dir, argv$data_file_csv, argv$out_dir, argv$out_prefix, argv$plot, argv$img)
} else {
    print("Invalid image option.")
    print.arg.parser()
}


# args = commandArgs(trailingOnly=TRUE)
# cmd_usage = "Usage:  Rscript scran_normalize.R in_base_dir data.file.csv out_dir out_prefix plot(TRE/FALSE) png/pdf"

# if(length(args) >= 6){
#     if((args[6] == "png" || args[6] == "pdf")){
#         scran_main(args[1], args[2], args[3], args[4], args[5], args[6])
#     } else {
#         print(args)
#         print(cmd_usage)
#     }
# }  else {
#     print(args)
#     print(cmd_usage)
# }
