library(Seurat)
source("qc_utils.R")
source("plot_utils.R")
source("data_utils.R")


combined_umap = function(athaliana, out.dir, reduce_by="pca",
			 dims=1:30, img.option="png"){
    #
    athaliana = cluster_umap_seurat(athaliana, reduce_by, 0.5, dims)
    #
    print(athaliana@reductions)
    dim_plot(athaliana, reduce_by="umap",group="dataset", split="dataset",
	     width=10, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-seurat-umap-grouped.", 
				  img.option, sep=""), 
			    sep="/"))
    #
    dim_plot(athaliana, reduce_by="umap",group=NULL, label=TRUE,
	     width=6, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-seurat-umap-integrated.", 
				  img.option, sep=""), 
			    sep="/"))
    athaliana
}


combined_tsne = function(athaliana, out.dir, reduce_by="pca", 
			 dims=1:30, img.option="png"){
    #
    athaliana = cluster_tsne_seurat(athaliana, reduce_by, 0.5, dims)
    #
    print(athaliana@reductions)
    dim_plot(athaliana, reduce_by="tsne",group="dataset", split="dataset",
	     width=10, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-seurat-tsne-grouped.",
				  img.option, sep=""), 
			    sep="/"))
    #
    dim_plot(athaliana, reduce_by="tsne",group=NULL, label=TRUE,
	     width=6, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-seurat-tsne-integrated.", 
				  img.option, sep=""), 
			    sep="/"))
    athaliana
}

seurat_cluster = function(root.dir, data.file, out.dir, qc.flag,
			  vis.option, img.option){
    data.df = read.csv(data.file, header=TRUE, stringsAsFactors=FALSE)
    expt.dir.paths = data.df$dir.paths
    short.names = data.df$short.names
    project.names = data.df$project.names

    athaliana.sobj = if(as.logical(qc.flag)){
        qcload_10X_seurat_objects(root.dir, expt.dir.paths, project.names, 
 				  short.names, qc_normalize_matrix)
    } else {
        load_10X_seurat_objects(root.dir, expt.dir.paths,
                                project.names, short.names)
    }

    athaliana.integrated = integrate_seurat_objects(athaliana.sobj)
    DefaultAssay(athaliana.integrated) <- "integrated"

    # Run the standard workflow for visualization and clustering
    athaliana.integrated <- ScaleData(athaliana.integrated, verbose = FALSE)
    athaliana.integrated <- RunPCA(athaliana.integrated, npcs = 30, verbose = FALSE)
    if(vis.option == "umap"){
        athaliana.integrated = combined_umap(athaliana.integrated, out.dir,
					     "pca", 1:30, img.option)
    }
    if(vis.option == "tsne"){
        athaliana.integrated = combined_tsne(athaliana.integrated, out.dir, 
					     "pca", 1:30, img.option)
    }
    mkdf = FindAllMarkers(athaliana.integrated)
    write.table(mkdf, paste(out.dir, "seurat-markers.tsv", sep="/"), 
		row.names=FALSE, sep="\t")
    athaliana.integrated
}

args = commandArgs(trailingOnly=TRUE)
cmd_usage = "Usage: Rscript seurat_cluster.R root_dir data.file.csv out_dir qc_flag tnse/umap png/pdf"
if(length(args) >= 6){
    if((args[5] == "tsne" || args[5] == "umap") && 
       (args[6] == "png" || args[6] == "pdf")){
        seurat_cluster(args[1], args[2], args[3], args[4], args[5], args[6])
    } else {
        print(args)
        print(cmd_usage)
    }
}  else {
    print(args)
    print(cmd_usage)
}
