library(Seurat)
library(cowplot)
library(ggplot2)
library(harmony)
source("data_utils.R")
source("qc_utils.R")
source("plot_utils.R")

combined_umap = function(athaliana, out.dir, reduce_by="harmony"){
    #
    athaliana = cluster_umap_seurat(athaliana, reduce_by, 0.5, 1:20)
    #
    print(athaliana@reductions)
    dim_plot(athaliana, reduce_by="umap",group="dataset", split="dataset",
	     width=10, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-dim-umap-grouped.png", sep=""), 
			    sep="/"))
    #
    dim_plot(athaliana, reduce_by="umap",group=NULL, label=TRUE,
	     width=6, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-dim-umap-integrated.png", sep=""),
			    sep="/"))
    athaliana
}

combined_tsne = function(athaliana, out.dir, reduce_by="harmony"){
    #
    athaliana = cluster_tsne_seurat(athaliana, reduce_by, 0.5, 1:20)
    #
    print(athaliana@reductions)
    dim_plot(athaliana, reduce_by="tsne",group="dataset", split="dataset",
	     width=10, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-dim-tsne-grouped.png", sep=""), 
			    sep="/"))
    #
    dim_plot(athaliana, reduce_by="tsne",group=NULL, label=TRUE,
	     width=6, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-dim-tsne-integrated.png", sep=""), 
			    sep="/"))
    athaliana
}

combined_pca = function(athalina, out.dir){
    print(athaliana@reductions)
    dim_plot(athaliana, reduce_by="pca",group="dataset", split="dataset",
	     width=10, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-dim-pca-grouped.png", sep=""),
			    sep="/"))
    #
    dim_plot(athaliana, reduce_by="pca",group=NULL, label=TRUE,
	     width=6, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "dim-pca-combined.png", sep=""),
			    sep="/"))
    athaliana
}

harmony_cluster = function(root.dir, data.file, out.dir, qc.flag, vis.option){
    data.df = read.csv(data.file, header=TRUE, stringsAsFactors=FALSE)
    print(head(data.df))
    expt.dir.paths = data.df[,'dir.paths']
    short.names = data.df[,'short.names']
    project.names = data.df[, 'project.names']
    
    athaliana.mlist = if(as.logical(qc.flag)){
        qcload_10X_matrices(root.dir, expt.dir.paths,
                            project.names, qc_normalize_matrix)
    } else {
        load_10X_matrices(root.dir, expt.dir.paths,
                           project.names)
    }

    athaliana = combined_seurat_object(athaliana.mlist, short.names)
    if(vis.option == "umap"){
       combined_umap(athaliana, out.dir, "pca")
    }
    if(vis.option == "tsne"){
       combined_tsne(athaliana, out.dir, "pca")
    }
#
#    dim_violin_plot(athaliana, "pca", "PC_1", "dataset", 
#		    paste(out.dir, "dim-pca-grouped.png",sep="/"))
    ath.list = integrate_data_harmony(athaliana)
    athaliana = ath.list[[1]]
    harmony_embeddings = ath.list[[2]]
#    dim_violin_plot(athaliana, "harmony", "harmony_1", "dataset",
#                    paste(out.dir, "dim-violin-harmony.png", sep="/"))
    athaliana = if(vis.option == "umap"){
     combined_umap(athaliana, out.dir, "harmony")
    }
    else {
        if(vis.option == "tsne"){
           combined_tsne(athaliana, out.dir, "harmony")
        } else {
		athaliana
	}
    }
    mkdf = FindAllMarkers(athaliana)
    write.table(mkdf, paste(out.dir, "harmony-markers.tsv", sep="/"), 
		row.names=FALSE, sep="\t")
    athaliana
}

args = commandArgs(trailingOnly=TRUE)
cmd_usage = "Usage: Rscript harmony_cluster.R root_dir data.file.csv out_dir qc_flag tnse/umap"
if(length(args) >= 5){
    if(args[5] == "tsne" || args[5] == "umap"){
        harmony_cluster(args[1], args[2], args[3], args[4], args[5])
    } else {
        print(args)
        print(cmd_usage)
    }
}  else {
    print(args)
    print(cmd_usage)
}

