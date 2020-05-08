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

# Parse the command line arguments
argv <- parse_args(p)

if(argv$img == "png" || argv$img == "pdf") {
    data_df = read.csv(argv$data_file_csv, header=TRUE, stringsAsFactors=FALSE)
    print(head(data_df))
    in_dirs = data_df[,'dir.paths']
    plot_qcstats(argv$root_dir, in_dirs, argv$out_dir, argv$out_prefix, argv$img)
} else {
    print("Invalid image option.")
    print(p)
}

