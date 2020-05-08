source("plot_utils.R")

seurat_qc_main = function(data.dir,plot.prefix, plot.suffix){
   seurat_allqc_plot(data.dir,plot.prefix, plot.suffix)
}

# Create a parser
p <- arg_parser("Generate Seurat qc plots for a given dir.")

# Add command line arguments
p <- add_argument(p, "input_dir", help="Input directory of datasets", type="character")
p <- add_argument(p, "out_dir", help="Output directory", type="character")
p <- add_argument(p, "--img", help="Output image option should be one png/pdf (default:png)", short='-g', default='png')

if(!(argv$img == "png" || argv$img == "pdf")) {
    seurat_qc_main(argv$in_base_dir,  argv$out_dir, argv$img)
} else {
    print("Invalid image option.")
    print.arg.parser()
}

# args = commandArgs(trailingOnly=TRUE)
# if(length(args) == 3){
#    main(args[1], args[2], args[3])
# } else {
#    cat("Need 3 args : data.dir, out.dir, png/pdf")
# }
