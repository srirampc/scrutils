
library(scran)
library("DropletUtils")
source("scater_qc.R")
source("seurat_plots.R")

scran_normalize = function(dfx){
    clusters <- quickCluster(dfx)
    dfx <- computeSumFactors(dfx, clusters=clusters)
    summary(sizeFactors(dfx))
    dfx <- logNormCounts(dfx)
    dfx
}
variance_plot = function(dfx, fname){
    dec = modelGeneVar(dfx)
    png(file=fname)
    plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
    curve(metadata(dec)$trend(x), col="blue", add=TRUE)
    dev.off()

}

scran_normalize_dir = function(out_dir, out_prefix, in_base_dir, dirx){
    dirx = paste(in_base_dir, dx, sep="/")
    dfx = read10xCounts(dirx)
    fname = paste(out_dir, dirx, 
        paste(out_prefix, "-variance-before.png", sep=""),
        sep="/"  )
    variance_plot(dfx, fname)
    fname = paste(out_dir, dirx,
        paste(out_prefix, "-seurat-scatter-before.png", sep=""),
    sep="/"  )
    scrj = CreateSeuratObject(counts = dfx, project = data.dir)
    seurat_fscatter(scrj, fname)

    dfx = apply_cell_filters(dfx)
    dfx = apply_gene_filters(dfx)
    dfx = scran_normalize(dfx)

    fname = paste(out_dir, dirx,
        paste(out_prefix, "-variance-after.png", sep=""),
    sep="/"  )
    variance_plot(dfx, fname)
    fname = paste(out_dir, dirx,
        paste(out_prefix, "-seurat-scatter-after.png", sep=""),
    sep="/"  )
    dev.off()
    scrj = CreateSeuratObject(counts = dfx, project = data.dir)
    seurat_fscatter(scrj, fname)
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) >= 4){
    for(dirx in args[4:length(args)]){
        normalize_dir(args[1], args[2], args[3], dirx)
    }
}  else {
    print(args)
    print("Usage: Rscript scran_normalize.R outdir out_prefix in_base_dir indir1 indir2 ...")
}
