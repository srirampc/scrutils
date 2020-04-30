library(Seurat)
library(cowplot)
library(ggplot2)
library(harmony)
source("data_utils.R")

qc_normalize = function(root.dir, expt.dir.path){
    rdx = paste(root.dir, expt.dir.path, sep="/")
    dfx = DropletUtils::read10xCounts(rdx)
    dfx = apply_cell_filters(dfx)
    dfx = apply_gene_filters(dfx)
    print(dfx)
    dfx = scran_normalize(dfx)
    dfx
}

dim_plot = function(sobj, reduce_by, group = "stim"){
    DimPlot(object = sobj, reduction = reduce_by, pt.size = .1, group.by = group)
}
violin_plot = function(sobj, feats, group) {
    VlnPlot(object = sobj, features = feats, group.by = group, pt.size = .1)
}

dim_violin_plot = function(sobj, reduce_by, feats, group, out_file = NULL){
    options(repr.plot.height = 5, repr.plot.width = 12)
    p1 = dim_plot(sobj, reduce_by, group)
    p2 = violin_plot(sobj, feats, group)
    p3 = plot_grid(p1, p2)
    if(is.null(out_file)){
        ggsave(out_file, p2, width=12, height=5)
    }
    p2
}

harmony_cluster = function(root.dir, data.file){
    data.df = read.csv(data.file, header=TRUE, stringsAsFactors=FALSE)
    print(head(data.df))
    expt.dir.paths = data.df[,'dir.paths']
    short.names = data.df[,'short.names']
    project.names = data.df[, 'project.names']
    print(expt.dir.paths)
    athaliana.mlist = load_10X_matrices(root.dir, expt.dir.paths,
                                        project.names)

    athaliana = combined_seurat_object(athaliana.mlist, short.names)

    dim_violin_plot(athaliana, "pca", "PC_1", "stim", "p3-sath1.png")
    ath.list = integrate_data_harmony(athaliana)
    athaliana = ath.list[[1]]
    harmony_embeddings = ath.list[[2]]
    #
    cluster_umap_seurat(athaliana, "harmony", 0.5, 1:20)
    dim_violin_plot(athaliana, "harmony", "harmony_1", "stim", "p3-sath2.png")

    #
    options(repr.plot.height = 4, repr.plot.width = 10)
    p3 = DimPlot(athaliana, reduction = "umap", group.by = "stim", pt.size = .1, split.by = 'stim')
    ggsave("p3-aths3.png", p3, width=10, height=4)
    #
    options(repr.plot.height = 4, repr.plot.width = 6)
    p3 = DimPlot(athaliana, reduction = "umap", label = TRUE, pt.size = .1)
    ggsave("p3-aths4.png", p3, width=6, height=4)
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) >= 2){
    harmony_cluster(args[1], args[2])
}  else {
    print(args)
    print("Usage: Rscript harmony_cluster.R root_dir data.file.csv")
}

