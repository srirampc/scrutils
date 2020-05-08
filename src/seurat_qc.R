library(argparser, quietly=T)
source("plot_utils.R")

seurat_qc_main = function(data.dir, out.dir, out.prefix, img.option){
   seurat_allqc_plot(data.dir, out.dir, out.prefix, img.option)
}

# Create a parser
p <- arg_parser("Generate Seurat qc plots for a given dataset.")

# Add command line arguments
p <- add_argument(p, "in_dir", help="Input directory of a dataset", type="character")
p <- add_argument(p, "out_dir", help="Output directory for the dataset", type="character")
p <- add_argument(p, "output_prefix", help="Output preifx of the files", type="character")
p <- add_argument(p, "--img", help="Output image option should be one png/pdf (default:png)", short='-g', default='png')

# Parse the command line arguments
argv <- parse_args(p)

if(argv$img == "png" || argv$img == "pdf") {
    print(argv$img)
    seurat_qc_main(argv$in_dir,  argv$out_dir, argv$out_prefix, argv$img)
} else {
    print("Invalid image option.")
    print(p)
}

