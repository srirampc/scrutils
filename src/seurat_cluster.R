library(argparser, quietly=TRUE)
library(Seurat, quietly=TRUE)
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
			              vis.option, img.option, include_genes_file,
                          exec_list_file){
    data.df = read.csv(data.file, header=TRUE, stringsAsFactors=FALSE)
    expt.dir.paths = data.df$dir.paths
    short.names = data.df$short.names
    project.names = data.df$project.names

    athaliana.sobj = if(as.logical(qc.flag)){
        qcload_10X_seurat_objects(root.dir, expt.dir.paths, project.names, 
                    short.names, qc_normalize_matrix,
                    include_genes_file, exclude_genes_file)
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
    write.table(mkdf, paste(out.dir, paste(vis.option, "-seurat-markers.tsv", 
					   sep=""),
			    sep="/"), 
		row.names=FALSE, sep="\t")
    athaliana.integrated
}

# Create a parser
p <- arg_parser("Pre-process, Normalize, Integrate w. Seurat and cluster")

# Add command line arguments
p <- add_argument(p, "root_dir", help="Root directory of datasets", type="character")
p <- add_argument(p, "data_file_csv", 
                  help="CSV file with dataset info (See ath.control.csv for example)",
                  type="character")
p <- add_argument(p, "out_dir", help="Output directory", type="character")
p <- add_argument(p, "--qc", help="Flag to indicate to preform qc", short='-q', default=TRUE)
p <- add_argument(p, "--vis", help="Visualization option should be tsne/umap (default:umap)", short='-v', default='umap')
p <- add_argument(p, "--img", help="Output image option should be one png/pdf (default:png)", short='-g', default='png')
p <- add_argument(p, "--exc", help="File containing list of inc. genes (default:None)", short='-e', default=NULL)
p <- add_argument(p, "--inc", help="File containing list of exc. genes (default:None)", short='-i', default=NULL)

# Parse the command line arguments
argv <- parse_args(p)

if(!((argv$vis == "tsne" || argv$vis == "umap") && 
     (argv$img == "png" || argv$img == "pdf"))) {
         seurat_cluster(argv$root_dir, argv$data_file_csv,
                argv$out_dir, argv$qc, argv$vis, argv$img, 
                argv$inc, argv$exec)
    
} else {
    print("Invalid image/visualization option.")
    print.arg.parser()
}



# args = commandArgs(trailingOnly=TRUE)
# cmd_usage = "Usage: Rscript seurat_cluster.R root_dir data.file.csv out_dir qc_flag tnse/umap png/pdf"
# if(length(args) >= 6){
#     if((args[5] == "tsne" || args[5] == "umap") && 
#        (args[6] == "png" || args[6] == "pdf")){
#         seurat_cluster(args[1], args[2], args[3], args[4], args[5], args[6])
#     } else {
#         print(args)
#         print(cmd_usage)
#     }
# }  else {
#     print(args)
#     print(cmd_usage)
# }
