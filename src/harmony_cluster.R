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

dim_plot = function(sobj, reduce_by, group = "dataset",
		    split=NULL, label=FALSE, 
		    width=6, height=4,
		    out_file=NULL){
    if(!is.null(out_file)) {
       options(repr.plot.height = height, repr.plot.width = width)
    }
    px = DimPlot(object = sobj, reduction = reduce_by, pt.size = .1, 
	    group.by = group, split.by=split, label=label)
    if(!is.null(out_file)) {
        ggsave(out_file, px, width=width, height=height)
    }
    px
}
violin_plot = function(sobj, feats, group = "dataset", 
		       width=6, height=4,
		       out_file=NULL) {
    if(!is.null(out_file)) {
       options(repr.plot.height = height, repr.plot.width = width)
    }
    px = VlnPlot(object = sobj, features = feats, group.by = group, pt.size = .1)
    if(!is.null(out_file)) {
        ggsave(out_file, px, width=width, height=height)
    }
    px
}

dim_violin_plot = function(sobj, reduce_by, feats, group, out_file = NULL){
    options(repr.plot.height = 5, repr.plot.width = 12)
    p1 = dim_plot(sobj, reduce_by, group)
    p2 = violin_plot(sobj, feats, group)
    p3 = plot_grid(p1, p2)
    if(!is.null(out_file)){
        ggsave(out_file, p2, width=12, height=5)
    }
    p2
}

harmony_umap = function(athaliana, out.dir){
    #
    athaliana = cluster_umap_seurat(athaliana, "harmony", 0.5, 1:20)
    #
    print(athaliana@reductions)
    dim_plot(athaliana, reduce_by="umap",group="dataset", split="dataset",
	     width=10, height=4,
	     out_file=paste(out.dir, "dim-umap-grouped.png", sep="/"))
    #
    dim_plot(athaliana, reduce_by="umap",group=NULL, label=TRUE,
	     width=6, height=4,
	     out_file=paste(out.dir, "dim-umap-integrated.png", sep="/"))
}

harmony_tsne = function(athaliana, out.dir){
    #
    athaliana = cluster_tsne_seurat(athaliana, "harmony", 0.5, 1:20)
    #
    print(athaliana@reductions)
    dim_plot(athaliana, reduce_by="tsne",group="dataset", split="dataset",
	     width=10, height=4,
	     out_file=paste(out.dir, "dim-tsne-grouped.png", sep="/"))
    #
    dim_plot(athaliana, reduce_by="tsne",group=NULL, label=TRUE,
	     width=6, height=4,
	     out_file=paste(out.dir, "dim-tsne-integrated.png", sep="/"))
}

harmony_cluster = function(root.dir, data.file, out.dir){
    data.df = read.csv(data.file, header=TRUE, stringsAsFactors=FALSE)
    print(head(data.df))
    expt.dir.paths = data.df[,'dir.paths']
    short.names = data.df[,'short.names']
    project.names = data.df[, 'project.names']
    
    athaliana.mlist = load_10X_matrices(root.dir, expt.dir.paths,
                                        project.names)

    athaliana = combined_seurat_object(athaliana.mlist, short.names)

    dim_violin_plot(athaliana, "pca", "PC_1", "dataset", 
		    paste(out.dir, "dim-pca-grouped.png",sep="/"))
    ath.list = integrate_data_harmony(athaliana)
    athaliana = ath.list[[1]]
    harmony_embeddings = ath.list[[2]]
    dim_violin_plot(athaliana, "harmony", "harmony_1", "dataset",
                    paste(out.dir, "dim-violin-harmony.png", sep="/"))
    #harmony_umap(athaliana, out.dir)
    harmony_tsne(athaliana, out.dir)
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) >= 3){
    harmony_cluster(args[1], args[2], args[3])
}  else {
    print(args)
    print("Usage: Rscript harmony_cluster.R root_dir data.file.csv out_dir")
}

