source("scater_qc.R")

args = commandArgs(trailingOnly=TRUE)

if(length(args) >= 3){
    apply_filters(args[1], args[2], args[3:length(args)])
}  else {
    print(args)
    print("Usage: Rscript qcfilter_workflow.R outdir out_prefix indir1 indir2 ...")
}
