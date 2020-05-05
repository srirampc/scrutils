source("plot_utils.R")


args = commandArgs(trailingOnly=TRUE)

if(length(args) >= 3){
    plot_qcstats(args[1], args[2], args[3:length(args)])
}  else {
    print(args)
    print("Usage: Rscript scater_plot_dir.R in_data_file outdir out_prefix ")
}
