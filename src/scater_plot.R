source("plot_utils.R")


args = commandArgs(trailingOnly=TRUE)

if((length(args) >= 5) && (args[5] == 'png' || args[5] == 'pdf')){
    in_base_dir=args[1]
    data.file=args[2]
    out_dir=args[3]
    out_prefix=args[4]
    image.option=args[5]

    data.df = read.csv(data.file, header=TRUE, stringsAsFactors=FALSE)
    print(head(data.df))
    in_dirs = data.df[,'dir.paths']
    plot_qcstats(in_base_dir, in_dirs, out_dir, out_prefix, image.option)
}  else {
    print(args)
    print("Usage: Rscript scater_plot_dir.R in_base_dir _data_file outdir out_prefix png/pdf")
}
