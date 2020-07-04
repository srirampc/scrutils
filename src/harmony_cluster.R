library(argparser, quietly=TRUE)
library(Seurat, quietly=TRUE)
library(cowplot, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(harmony, quietly=TRUE)
source("data_utils.R")
source("qc_utils.R")
source("plot_utils.R")

combined_umap = function(athaliana, out.dir, reduce_by="harmony",
                         image.option="png"){
    #
    athaliana = cluster_umap_seurat(athaliana, reduce_by, 0.5, 1:20)
    #
    print(athaliana@reductions)
    dim_plot(athaliana, reduce_by="umap",group="dataset", split="dataset",
	     width=10, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-umap-grouped.",
                      image.option, sep=""), 
			    sep="/"))
    #
    dim_plot(athaliana, reduce_by="umap",group=NULL, label=TRUE,
	     width=6, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-umap-integrated.",
                      image.option, sep=""),
			    sep="/"))
    athaliana
}

combined_tsne = function(athaliana, out.dir,
                         reduce_by="harmony", image.option="png"){
    #
    athaliana = cluster_tsne_seurat(athaliana, reduce_by, 0.5, 1:20)
    #
    #print(athaliana@reductions)
    dim_plot(athaliana, reduce_by="tsne",group="dataset", split="dataset",
	     width=10, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-tsne-grouped.",
                      image.option, sep=""), 
			    sep="/"))
    #
    dim_plot(athaliana, reduce_by="tsne",group=NULL, label=TRUE,
	     width=6, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-tsne-integrated.",
                      image.option, sep=""), 
			    sep="/"))
    athaliana
}

combined_pca = function(athalina, out.dir, image.option="png"){
    #print(athaliana@reductions)
    dim_plot(athaliana, reduce_by="pca",group="dataset", split="dataset",
	     width=10, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-pca-grouped.",
                      image.option, sep=""),
			    sep="/"))
    #
    dim_plot(athaliana, reduce_by="pca",group=NULL, label=TRUE,
	     width=6, height=4,
	     out_file=paste(out.dir, 
			    paste(reduce_by, "-pca-combined.",
                      image.option, sep=""),
			    sep="/"))
    athaliana
}

harmony_cluster = function(root.dir, data.file, out.dir, 
                            qc.flag, vis.option, img.option,
                            inc_list_file, exec_list_file,
                            gen_markers, dot_markers_file){
    data.df = read.csv(data.file, header=TRUE, stringsAsFactors=FALSE)
    print(head(data.df))
    expt.dir.paths = data.df[,'dir.paths']
    short.names = data.df[,'short.names']
    project.names = data.df[, 'project.names']
    
    athaliana.mlist = if(as.logical(qc.flag)){
        qcload_10X_matrices_union(root.dir, expt.dir.paths,
                            project.names, qc_normalize_matrix,
                            inc_list_file, exec_list_file)
    } else {
        load_10X_matrices(root.dir, expt.dir.paths,
                           project.names)
    }

    athaliana = combined_seurat_object(athaliana.mlist, short.names)
    if(vis.option == "umap"){
        combined_umap(athaliana, out.dir, "pca", img.option)
    }
    if(vis.option == "tsne"){
        combined_tsne(athaliana, out.dir, "pca", img.option)
    }
    #
    dim_violin_plot(athaliana, "pca", "PC_1", "dataset", 
            paste(out.dir, 
                paste("violin-pca-grouped.", img.option, sep=""),
                sep="/"))
    ath.list = integrate_data_harmony(athaliana)
    athaliana = ath.list[[1]]
    harmony_embeddings = ath.list[[2]]
    dim_violin_plot(athaliana, "harmony", "harmony_1", "dataset",
                    paste(out.dir, 
                    paste("violin-harmony.", img.option, sep=""), 
                    sep="/"))
    athaliana = if(vis.option == "umap"){
        combined_umap(athaliana, out.dir, "harmony", img.option)
    } else {
        if(vis.option == "tsne"){
           combined_tsne(athaliana, out.dir, "harmony", img.option)
        } else {
            athaliana
        }
    }
    print(athaliana@reductions)
    #print(athaliana@seurat_clusters)

    if(gen_markers) {
        mkdf = FindAllMarkers(athaliana)
        write.table(mkdf, paste(out.dir, "harmony-markers.tsv", sep="/"), 
                    row.names=FALSE, sep="\t")
    }

    if(!is.na(dot_markers_file)){
        mdf = read.table(dot_markers_file, stringsAsFactors=F, sep="\t")
        for(dxf in mdf$V1){
          gdf = read.table(dxf, header=T, 
                           sep="\t", stringsAsFactors=F)
          genes.plot = gdf$ID
          pfx = gsub("\\.tsv", "" , basename(dxf))
          npresent = genes.plot %in% rownames(athaliana)
          cat(pfx, length(genes.plot), sum(npresent), sum(!npresent), "\n")
          if(sum(npresent) > 0) {
              dot_fname = paste(out.dir, 
                          paste(pfx, "-dot-markers-harmony.", 
                              img.option, sep=""), 
                       sep="/")
              px = DotPlot(athaliana, features=genes.plot) + 
                  ggtitle(dxf)
               ggsave(dot_fname, px, width=10, height=4)
           }
        }
    }
    athaliana
}

# Create a parser
p <- arg_parser("Pre-process, Normalize, Integrate w. Harmony and Cluster")

# Add command line arguments
p <- add_argument(p, "root_dir", help="Root directory of datasets", type="character")
p <- add_argument(p, "data_file_csv", 
                  help="CSV file with dataset info (See ath.control.csv for example)",
                  type="character")
p <- add_argument(p, "out_dir", help="Output directory", type="character")
p <- add_argument(p, "--qc", help="Flag to indicate to preform qc", short='-q', default=TRUE)
p <- add_argument(p, "--vis", help="Visualization option should be tsne/umap ", short='-v', default='umap')
p <- add_argument(p, "--img", help="Output image option should be one png/pdf ", short='-g', default='png')
p <- add_argument(p, "--exc [glist]", help="File containing list of inc. genes (default:None)", short='-e', default=NULL)
p <- add_argument(p, "--inc [glist]", help="File containing list of exc. genes (default:None)", short='-i', default=NULL)
p <- add_argument(p, "--gen_markers", help="Flag to indicate to genreate all markers", short='-m', default=FALSE)
p <- add_argument(p, "--dot_markers [mfile]", help="Generate marker dot plots for each list of markers list in mfile. Each line in mfile is location for a list of markers", short='-d', default=NULL)

# Parse the command line arguments
argv <- parse_args(p)

if((argv$vis == "tsne" || argv$vis == "umap") && 
   (argv$img == "png" || argv$img == "pdf")) {
         harmony_cluster(argv$root_dir, argv$data_file_csv,
                argv$out_dir, argv$qc, argv$vis, argv$img, 
                argv$inc, argv$exec, 
                argv$gen_markers, argv$dot_markers)
    
} else {
    print("Invalid image/visualization option.")
    print(p)
}

