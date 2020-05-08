library(argparser, quietly=T)
source("plot_utils.R")


# Create a parser
p <- arg_parser("Generate qc plots using scater library")

# Add command line arguments
p <- add_argument(p, "root_dir", help="Root directory of datasets", type="character")
p <- add_argument(p, "data_file_csv", 
                  help="CSV file with dataset info (See ath.control.csv for example)",
                  type="character")
p <- add_argument(p, "out_dir", help="Output directory", type="character")
p <- add_argument(p, "out_prefix", help="Output Prefix", type="character")
p <- add_argument(p, "--img", help="Output image option should be one png/pdf (default:png)", short='-g', default='png')

if(!(argv$img == "png" || argv$img == "pdf")) {
    data_df = read.csv(argv$data_file_csv, header=TRUE, stringsAsFactors=FALSE)
    print(head(data_df))
    in_dirs = data_df[,'dir.paths']
    plot_qcstats(argv$in_base_dir, in_dirs, argv$out_dir, argv$out_prefix, argv$img)
} else {
    print("Invalid image option.")
    print.arg.parser()
}
# args = commandArgs(trailingOnly=TRUE)

# if((length(args) >= 5) && (args[5] == 'png' || args[5] == 'pdf')){
#     in_base_dir=args[1]
#     data.file=args[2]
#     out_dir=args[3]
#     out_prefix=args[4]
#     image.option=args[5]

#     data.df = read.csv(data.file, header=TRUE, stringsAsFactors=FALSE)
#     print(head(data.df))
#     in_dirs = data.df[,'dir.paths']
#     plot_qcstats(in_base_dir, in_dirs, out_dir, out_prefix, image.option)
# }  else {
#     print(args)
#     print("Usage: Rscript scater_plot_dir.R in_base_dir _data_file outdir out_prefix png/pdf")
# }
