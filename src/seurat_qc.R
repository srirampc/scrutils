source("seurat_plots.R")

main = function(data.dir,plot.prefix, plot.suffix){
   seurat_allqc_plot(data.dir,plot.prefix, plot.suffix)
}


args = commandArgs(trailingOnly=TRUE)
if(length(args) == 3){
   main(args[1], args[2], args[3])
} else {
   cat("Need 3 args : data.dir, out.dir, plot.suffix")
}
